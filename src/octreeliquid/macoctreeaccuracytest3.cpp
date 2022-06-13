/*
**	macoctreeaccuracytest3.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on November 17, 2019. All rights reserved.
**	Permission is hereby granted, free of charge, to any person obtaining a copy of
**	this software and associated documentation files (the "Software"), to deal in
**	the Software without restriction, including without limitation the rights to use,
**	copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the
**	Software, and to permit persons to whom the Software is furnished to do so,
**	subject to the following conditions:
**
**	The above copyright notice and this permission notice shall be included in all copies
**	or substantial portions of the Software.
**
**	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
**	INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
**	PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
**	HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
**	CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
**	OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
//
#ifndef SHKZ_OCTREEACCURACYTEST3_H
#define SHKZ_OCTREEACCURACYTEST3_H
//
#include "macoctreegrid3.h"
#include "macoctreeproject3.h"
#include "macoctreemesher3.h"
#include "mitsuba_xml.h"
#include <shiokaze/ui/drawable.h>
#include <shiokaze/image/color.h>
#include <shiokaze/core/filesystem.h>
#include <shiokaze/core/console.h>
#include <shiokaze/meshexporter/meshexporter3_interface.h>
#include <cmath>
#include <mutex>
#include <typeinfo>
//
SHKZ_BEGIN_NAMESPACE
//
namespace macotreeliquid3_namespace {
	//
	class macoctreeaccuracytest3 : public drawable {
	protected:
		//
		LONG_NAME("MAC Octree Accuracy Test 3D")
		ARGUMENT_NAME("macoctreeaccuracytest")
		//
		virtual void configure( configuration &config ) override {
			//
			m_dx = m_shape.dx();
			//
			set_environment("shape",&m_shape);
			set_environment("dx",&m_dx);
			//
			config.get_integer("Case",m_case,"Case number");
			config.set_default_double("Residual",1e-18);
			config.set_default_double("EpsFluid",0.0);
			config.set_default_double("EpsSolid",0.0);
			config.set_default_double("AccuracyEps",0.0);
			config.set_default_bool("SteepAdaptivity",false);
			config.set_default_bool("DegradeWarning",true);
			config.set_default_bool("VolumeCorrection",false);
			config.set_default_unsigned("DilateCount",0);
			config.set_default_unsigned("MaxIterations",300000);
			config.set_default_bool("RemoveOneDegreesOfFreedom",true);
			//
			config.get_double("Radius",m_r,"Circle radius");
			config.get_unsigned("SampleCount",m_sample_count,"Sample count for rendering");
			config.get_bool("RenderMovie",m_render_movie,"Render movie");
			config.get_bool("RenderImage",m_render_image,"Render image");
			config.get_unsigned("ImageWidth",m_canvas_shape[0],"Image width");
			config.get_unsigned("ImageHeight",m_canvas_shape[1],"Image height");
			config.get_integer("TrialCount",m_max_trial_count,"Number of trial count");
			config.get_integer("SubdivisionCount",m_max_subdivion_count,"Grid subdivision count");
			config.get_vec3d("Shift",m_shift.v,"Shift amount");
			config.get_bool("UniformAdaptivity",m_uniform_adaptivity,"Use uniform grid");
		}
		//
		double analytical_function( const vec3d &center, const vec3d &p, double r ) {
			//
			if( m_case == 0 ) {
				double d2 = (p-center).norm2();
				return d2-r*r;
			} else if( m_case == 1 ) {
				const vec3d q = p-center;
				const double x = 2.0*q[0];
				const double y = 2.0*q[1];
				const double z = 2.0*q[2];
				return x*(y*y*y)+(x*x)*(z*z);
			}
			return 0.0;
		}
		//
		vec3d derivative_analytical_function( const vec3d &center, const vec3d &p ) {
			//
			vec3d r = p-center;
			if( m_case == 0 ) {
				return 2.0*r;
			} else if( m_case == 1 ) {
				const double x = 2.0*r[0];
				const double y = 2.0*r[1];
				const double z = 2.0*r[2];
				return vec3d(y+(y*y*y)+2.0*x*(z*z),-x+3.0*x*(y*y),2.0*(x*x)*z);
			}
			return vec3d();
		}
		//
		double measure_surface_norm ( const vec3d &center, double r, unsigned depth ) {
			//
			double inf_norm (0.0);
			if( m_case == 0 ) {
				m_grid.serial_iterate_active_cells([&]( const cell_id3 &cell_id ) {
					if( cell_id.depth == depth ) {
						if( m_grid.levelset[cell_id.index] < 0.0 ) {
							double analytical_value = analytical_function(center,m_grid.get_cell_position(cell_id),r);
							double error = std::abs(analytical_value-m_pressure[cell_id.index]);
							if( error > inf_norm ) {
								inf_norm = error;
							}
						}
					}
				});
			} else if( m_case == 1 ) {
				//
				auto check_valid_cell = [&]( const cell_id3 &cell_id ) {
					bool valid (false);
					for( char dim : DIMS3 ) {
						m_grid.iterate_face_neighbors(cell_id,dim,[&]( const face_id3 &face_id ) {
							if( m_grid.area[face_id.index] ) valid = true;
						});
					}
					return valid;
				};
				//
				double analytical_average (0.0);
				double discrete_average (0.0);
				unsigned sample_sum (0);
				//
				m_grid.serial_iterate_active_cells([&]( const cell_id3 &cell_id ) {
					if( cell_id.depth == depth && check_valid_cell(cell_id)) {
						const vec3d &p = m_grid.get_cell_position(cell_id);
						analytical_average += analytical_function(center,p,r);
						discrete_average += m_pressure[cell_id.index];
						sample_sum ++;
					}
				});
				analytical_average /= (double)sample_sum;
				discrete_average /= (double)sample_sum;
				//
				m_grid.serial_iterate_active_cells([&]( const cell_id3 &cell_id ) {
					if( cell_id.depth == depth && check_valid_cell(cell_id)) {
						const vec3d &p = m_grid.get_cell_position(cell_id);
						const double v0 = 0.5*(analytical_function(center,p,r)-analytical_average);
						const double v1 = m_pressure[cell_id.index]-discrete_average;
						const double error = std::abs(v0-v1);
						if( error > inf_norm ) {
							inf_norm = error;
						}
					}
				});
			}
			return inf_norm;
		}
		//
		void export_cutaway_mesh () {
			//
			unsigned prefix_num = m_shape[0];
			//
			std::string export_path = console::get_root_path() + "/mesh";
			if( ! filesystem::is_exist(export_path)) {
				filesystem::create_directory(export_path);
			}
			//
			auto remove_half = [&]( std::vector<vec3d> &vertices, std::vector<std::vector<size_t> > &faces, int mode) {
				//
				for( auto it=faces.begin(); it!=faces.end(); ) {
					vec3d center;
					int sum (0);
					double x_max (0.0);
					for( const auto &n : *it ) {
						center += vertices[n];
						x_max = std::max(x_max,vertices[n][0]);
						sum ++;
					}
					center /= sum;
					const double x_bound (0.5-m_dx);
					if( mode == 0 ) {
						if( center[2] > 0.5 && x_max > x_bound ) {
							it = faces.erase(it);
						} else {
							for( const auto &n : *it ) {
								if( center[2] > 0.5 && vertices[n][0] > x_bound-m_dx ) vertices[n][0] = x_bound-m_dx;
							}
							++it;
						}
					} else if( mode == 1 ) {
						if( center[2] > 0.5 ) {
							it = faces.erase(it);
						} else {
							++it;
						}
					}
				}
			};
			//
			// XML file export
			std::vector<FILE *> files(3);
			for( int k=0; k<3; ++k ) {
				files[k] = fopen((export_path+"/"+std::to_string(prefix_num)+"_cutaway_mesh"+std::to_string(k)+".xml").c_str(),"w");
			}
			//
			for( int k=0; k<2; ++k ) {
				assert(files[k]);
				fprintf(files[k],mitsuba_xml::header.c_str());
				fprintf(files[k],mitsuba_xml::sensor.c_str(),m_sample_count);
				fprintf(files[k],mitsuba_xml::sunlight.c_str());
				if( k == 0 && ! m_render_movie ) fprintf(files[k],mitsuba_xml::inlight.c_str());
				fprintf(files[k],mitsuba_xml::floor.c_str());
			}
			//
			// Export fluid mesh
			std::vector<vec3d> vertices;
			std::vector<std::vector<size_t> > faces;
			m_macoctreemesher.generate_mesh(m_grid,m_grid.levelset,0.0,nullptr,vertices,faces,false);
			remove_half(vertices,faces,0);
			//
			m_mesh_exporter->set_mesh(vertices,faces);
			m_mesh_exporter->export_ply(console::format_str("%s/%d_cutaway_mesh.ply",export_path.c_str(),prefix_num));
			m_mesh_exporter->export_mitsuba(console::format_str("%s/%d_cutaway_mesh.serialized",export_path.c_str(),prefix_num));
			//
			for( int k=0; k<1; ++k ) {
				fprintf(files[k],mitsuba_xml::wireframe_mesh.c_str(),console::format_str("%d_cutaway_mesh.serialized",prefix_num).c_str(),0.5,0.6,1.0);
			}
			//
			int max_count (6);
			double theta = 1.0 / max_count;
			double rgb_colors[][3] = {{0.0,0.3,1.0},{0.2,0.4,0.8},{0.5,0.3,0.6},{0.7,0.4,0.7},{0.8,0.4,0.4},{1.0,0.2,0.2}};
			//
			for( int n=0; n<max_count; ++n ) {
				//
				m_macoctreemesher.generate_mesh(m_grid,m_pressure,0.025*(n+1),nullptr,vertices,faces,false);
				if( n < max_count-1 ) remove_half(vertices,faces,1);
				//
				m_mesh_exporter->set_mesh(vertices,faces);
				m_mesh_exporter->export_ply(console::format_str("%s/%d_cutaway_mesh_level%d.ply",export_path.c_str(),prefix_num,n));
				m_mesh_exporter->export_mitsuba(console::format_str("%s/%d_cutaway_mesh_level%d.serialized",export_path.c_str(),prefix_num,n));
				//
				for( int k=0; k<2; ++k ) {
					fprintf(files[k],mitsuba_xml::mesh.c_str(),console::format_str("%d_cutaway_mesh_level%d.serialized",prefix_num,n).c_str(),
					rgb_colors[n][0],rgb_colors[n][1],rgb_colors[n][2]);
				}
			}
			//
			for( int k=0; k<2; ++k ) {
				fprintf(files[k],mitsuba_xml::footer.c_str());
				fclose(files[k]);
			}
		}
		//
		virtual void post_initialize( bool initialized_from_file ) override {
			//
			if( typeid(Real) == typeid(double)) {
				console::dump( "Real = double\n");
			} else if( typeid(Real) == typeid(float)) {
				console::dump( "Real = float\n");
			}
			//
			m_grid.clear();
			console::set_time(m_shape[0]);
			//
			unsigned n (1);
			for( unsigned k=0; k<4; ++k ) {
				shape3 shape = m_shape/n;
				m_grid.add_layer(shape,n*m_dx);
				n *= 2;
			}
			//
			const vec3d center = vec3d(0.5,0.5,0.5)+m_shift;
			auto setup = [&]( double r ) {
				//
				auto fluid_func = [&](const vec3d &p) { 
					if( m_case == 0 ) {
						return (p-center).len()-r;
					} else {
						return -1.0;
					}
				};
				auto solid_func = [&](const vec3d &p) {
					if( m_case == 0 ) {
						return 1.0;
					} else if( m_case == 1 ) {
						return r-(p-center).len();
					}
					return 1.0;
				};
				auto adaptivity_func = [&](char depth, const vec3d &p) {
					if( m_uniform_adaptivity ) {
						return true;
					} else {
						if( fluid_func(p)-3.0*m_grid.layers[depth]->dx < 0.0 ) {
							return 4*p[0] < depth+1;
						} else {
							return false;
						}
					}
				};
				//
				m_grid.activate_cells(adaptivity_func);
				m_grid.balance_layers();
				m_grid.assign_indices();
				m_grid.assign_levelset(fluid_func,solid_func);
				m_grid.set_velocity([&]( const vec3d &p, char dim ) { return derivative_analytical_function(center,p)[dim]; });
				m_grid.serial_iterate_active_faces([&]( const face_id3 &face_id ) {
					if( ! m_grid.area[face_id.index] ) m_grid.velocity[face_id.index] = 0.0;
				});
				m_macoctreeproject.assemble_matrix(m_grid);
				m_macoctreeproject.project(m_grid,1.0,nullptr,&m_pressure);
			};
			//
			unsigned size = m_grid.layers.size();
			std::vector<double> max_norm (size,0.0);
			size_t compromised_count (0), total_count (0);
			std::vector<int> q_list;
			for( int q=-m_max_trial_count; q<=m_max_trial_count; ++q ) {
				if( q ) q_list.push_back(q);
			}
			q_list.push_back(0);
			for( int q : q_list ) {
				//
				const double t = 2.0*sqrt(3.0)*m_dx/m_max_trial_count;
				const double r = m_r+t*q;
				setup(r);
				//
				for( unsigned depth=0; depth<size; ++depth ) {
					double inf_norm = measure_surface_norm(center,r,depth);
					console::dump( "R=%d, q=%d, depth=%u, r=%e, inf_norm = %.2e\n", (int)m_shape[0], q, depth, r, inf_norm );
					if( max_norm[depth] < inf_norm ) {
						max_norm[depth] = inf_norm;
						compromised_count = m_grid.get_compromised_gradient_count();
						total_count = m_grid.get_surface_T_junction_count();
					}
				}
			}
			for( unsigned depth=0; depth<size; ++depth ) {
				console::write(std::string("max_norm")+std::to_string(depth),max_norm[depth]);
			}
			//
			static std::vector<double> prev_norm;
			console::dump( "========== Result ===========\n" );
			console::dump( "compromised_count = %u\n", compromised_count );
			console::dump( "total_count = %u\n", total_count );
			console::write("num_compromised",compromised_count);
			console::write("num_surface_t_junction",total_count);
			for( unsigned depth=0; depth<size; ++depth ) {
				if( prev_norm.empty()) {
					console::dump( "R=%d, depth=%u, max_norm = %.2e\n", (int)m_shape[0], depth, max_norm[depth] );
				} else {
					double order = log(prev_norm[depth]/max_norm[depth])/log(2);
					console::dump( "R=%d, depth=%u, max_norm = %.2e (order=%.2f)\n", (int)m_shape[0], depth, max_norm[depth], order );
					console::write(std::string("order")+std::to_string(depth),order);
				}
			}
			prev_norm = max_norm;
			//
			if( console::get_root_path().size()) {
				export_cutaway_mesh();
				if( console::system("mitsuba > /dev/null 2>&1") == 0 ) {
					//
					unsigned prefix_num = m_shape[0];
					std::string export_path = console::get_root_path() + "/mesh";
					//
					if( m_render_movie ) {
						std::string mov_path = console::get_root_path() + "/mov";
						if( ! filesystem::is_exist(mov_path)) {
							filesystem::create_directory(mov_path);
						}
						filesystem::create_directory(mov_path+"/"+std::to_string(prefix_num));
						const double dt (0.02);
						const double r (2.7+m_r);
						for( unsigned count=0; count<2.0*M_PI/dt; ++count ) {
							//
							const double t = count * dt;
							const double x = 0.5+r*sin(t);
							const double y = 1.5;
							const double z = 0.5+r*cos(t);
							//
							for( int k=0; k<2; ++k ) {
								std::string output_dir = console::format_str(mov_path+"/%d/%d",prefix_num,k);
								if( ! filesystem::is_exist(output_dir)) {
									filesystem::create_directory(output_dir);
								}
								for( unsigned count0=0; count0<=count; ++count0 ) {
									std::string dest_exr_path = console::format_str("%s/%d_%d_cutaway_mesh%d.exr", output_dir.c_str(), count0, prefix_num, k );
									std::string dest_exr_png = console::format_str("%s/%d_%d_cutaway_mesh%d.png", output_dir.c_str(), count0, prefix_num, k );
									if( ! filesystem::is_exist(dest_exr_path)) {
										console::system("mitsuba -Dx=%g -Dy=%g -Dz=%g, -Dw=%d -Dh=%d -o %s %s/%d_cutaway_mesh%d.xml;",
														x, y, z, 1280, 720, dest_exr_path.c_str(), export_path.c_str(), prefix_num, k);
									}
									if( ! filesystem::is_exist(dest_exr_png)) {
										console::system("mtsutil tonemap -g 3 %s", dest_exr_path.c_str());
									}
								}
								//
								if( k == 1 ) {
									for( unsigned count0=0; count0<=count; ++count0 ) {
										std::string output_dir1 = console::format_str(mov_path+"/%d/%d",prefix_num,0);
										std::string output_dir2 = console::format_str(mov_path+"/%d/%d",prefix_num,1);
										std::string dest_png = console::format_str("%s/%d_%d_composite.png",output_dir.c_str(), count0, prefix_num);
										if( ! filesystem::is_exist(dest_png)) {
											console::system("composite -blend 50 %s/%d_%d_cutaway_mesh0.png %s/%d_%d_cutaway_mesh1.png %s",
														output_dir1.c_str(), count0, prefix_num, output_dir2.c_str(), count0, prefix_num, dest_png.c_str() );
										}
									}
								}
							}
						}
						//
						for( int k=0; k<2; ++k ) {
							std::string output_dir = console::format_str(mov_path+"/%d/%d",prefix_num,k);
							console::system("avconv -r 60 -i %s/%%d_%d_cutaway_mesh%d.png -b:v 12000k -pix_fmt yuv420p %s/%d_cutaway_mesh%d.mp4",
								output_dir.c_str(), prefix_num, k, mov_path.c_str(), prefix_num, k );
							if( k == 1 ) {
								console::system("avconv -r 60 -i %s/%%d_%d_composite.png -b:v 12000k -pix_fmt yuv420p %s/%d_composite.mp4",
												output_dir.c_str(), prefix_num, mov_path.c_str(), prefix_num );
							}
						}
						//
					} else if( m_render_image ) {
						//
						std::string img_path = console::get_root_path() + "/img";
						if( ! filesystem::is_exist(img_path)) {
							filesystem::create_directory(img_path);
						}
						std::string export_path = console::get_root_path() + "/mesh";
						for( int k=0; k<2; ++k ) {
							console::system("mitsuba -Dx=0.5 -Dy=1.0 -Dz=2.3, -Dw=%d -Dh=%d -o %s/%d_cutaway_mesh%d.exr %s/%d_cutaway_mesh%d.xml;"
											"mtsutil tonemap -g 3 %s/%d_cutaway_mesh%d.exr;",
							m_canvas_shape[0], m_canvas_shape[1], img_path.c_str(), prefix_num, k, export_path.c_str(), prefix_num, k,
							img_path.c_str(), prefix_num, k, prefix_num, k );
						}
						console::system("composite -blend 50 %s/%d_cutaway_mesh0.png %s/%d_cutaway_mesh1.png %s/%d_composite.png",
							img_path.c_str(), prefix_num, img_path.c_str(), prefix_num, img_path.c_str(), prefix_num );
					}
				}
			}
			//
			m_camera->set_bounding_box(vec3d().v,m_shape.box(m_dx).v);
		}
		//
		virtual void idle () override {
			if( ! UI_interface::has_graphical_interface() && console::get_root_path().size()) {
				if( m_step < m_max_subdivion_count ) {
					m_dx *= 0.5;
					m_shape *= 2.0;
					m_step ++;
					reinitialize();
				}
			}
		}
		//
		virtual bool should_quit() const override { return m_step >= m_max_subdivion_count; }
		//
		virtual bool keyboard( int key, int action, int mods ) override {
			//
			if( action == UI_interface::PRESS && key == KEY_C ) {
				//
				m_dx *= 0.5;
				m_shape *= 2.0;
				m_step ++;
				reinitialize();
				return true;
			}
			return drawable::keyboard(key,action,mods);
		}
		//
		void draw_fluid( graphics_engine &g ) const {
			//
			// Draw surface
			std::vector<vec3d> vertices;
			std::vector<std::vector<size_t> > faces;
			m_macoctreemesher.generate_mesh(m_grid,m_grid.levelset,0.0,nullptr,vertices,faces,false);
			//
			g.color4(1.0,1.0,1.0,0.5);
			for( unsigned i=0; i<faces.size(); i++ ) {
				g.begin(graphics_engine::MODE::LINE_LOOP);
				for( unsigned j=0; j<faces[i].size(); j++ ) g.vertex3v(vertices[faces[i][j]].v);
				g.end();
			}
		}
		//
		void draw_pressure( graphics_engine &g ) const {
			//
			int max_count (6);
			double theta = 1.0 / max_count;
			//
			for( int n=0; n<max_count; ++n ) {
				//
				double rgb[4];
				color::heatcolor((n+1)*theta,rgb); rgb[3] = 1.0;
				//
				std::vector<vec3d> vertices;
				std::vector<std::vector<size_t> > faces;
				m_macoctreemesher.generate_mesh(m_grid,m_pressure,0.025*(n+1),nullptr,vertices,faces,false);
				//
				g.color4v(rgb);
				for( unsigned i=0; i<faces.size(); i++ ) {
					g.begin(graphics_engine::MODE::LINE_LOOP);
					for( unsigned j=0; j<faces[i].size(); j++ ) g.vertex3v(vertices[faces[i][j]].v);
					g.end();
				}
			}
		}
		//
		virtual void draw( graphics_engine &g ) const override {
			//
			if( m_case == 0 ) {
				draw_fluid(g);
				draw_pressure(g);
			}
		}
		//
		virtual void setup_window( std::string &name, int &width, int &height ) const override {
			double ratio = m_shape[1] / (double) m_shape[0];
			height = ratio * width;
		}
		//
		struct NormPoint3 {
			vec3d p;
			double dx;
		};
		//
		shape3 m_shape {64,64,64};
		shape2 m_canvas_shape {1800,1800};
		double m_dx;
		double m_r {0.4};
		vec3d m_shift;
		unsigned m_step {0};
		unsigned m_sample_count {64};
		int m_max_trial_count {4};
		int m_max_subdivion_count {4};
		int m_case {0};
		bool m_uniform_adaptivity {false};
		bool m_render_movie {false};
		bool m_render_image {true};
		std::vector<Real> m_pressure;
		//
		grid3 m_grid{this};
		macoctreeproject3 m_macoctreeproject{this};
		macoctreemesher3 m_macoctreemesher{this};
		meshexporter3_driver m_mesh_exporter{this,"meshexporter3"};
		//
	};
	//
	extern "C" module * create_instance() {
		return new macoctreeaccuracytest3;
	}
};
//
SHKZ_END_NAMESPACE
//
#endif