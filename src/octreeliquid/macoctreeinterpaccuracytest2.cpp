/*
**	macoctreeinterpaccuracytest2.cpp
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
#ifndef SHKZ_OCTREEINTERPACCURACY2_H
#define SHKZ_OCTREEINTERPACCURACY2_H
//
#include "macoctreegrid2.h"
#include "mitsuba_xml.h"
#include <shiokaze/ui/drawable.h>
#include <shiokaze/graphics/graphics_interface.h>
#include <shiokaze/core/filesystem.h>
#include <shiokaze/core/console.h>
#include <shiokaze/meshexporter/meshexporter3_interface.h>
#include <shiokaze/visualizer/gridvisualizer2_interface.h>
#include <shiokaze/array/shared_array2.h>
#include <shiokaze/image/color.h>
#include <cmath>
#include <mutex>
#include <typeinfo>
//
SHKZ_BEGIN_NAMESPACE
//
namespace macotreeliquid2_namespace {
	//
	class macoctreeinterpaccuracytest2 : public drawable {
	protected:
		//
		LONG_NAME("MAC Octree Interpolation Accuracy 2D")
		ARGUMENT_NAME("macinterpoctreeaccuracy")
		//
		virtual void configure( configuration &config ) override {
			//
			config.get_double("Speed",m_param.speed,"Speed");
			config.get_double("Wavelength",m_param.wavelen,"Wave length");
			config.get_double("Height",m_param.height,"Height");
			config.get_unsigned("SampleCount",m_param.sample_count,"SampleCount");
			config.get_integer("Cycle",m_param.cycle,"Function example number");
			config.get_unsigned("MaxSteps",m_param.max_steps,"Max steps to record");
			config.get_unsigned("UpsampleScale",m_param.upsampling_scale,"Upsample scale");
			config.get_bool("RenderMesh",m_param.render_mesh,"Render mesh");
			config.get_integer("ContourCount",m_param.contour_count,"Contour count");
			config.get_bool("WriteSVG",m_param.write_svg,"Write SVG");
			config.get_integer("Mode",m_mode,"0 (cell) or 1 (face)");
			config.get_integer("Dimension",m_dim,"0 (x) or 1 (y)");
			//
			double resolution_scale (1.0);
			config.get_double("ResolutionScale",resolution_scale,"Resolution doubling scale");
			//
			m_shape *= resolution_scale;
			m_high_shape = m_param.upsampling_scale * m_shape;
			//
			m_dx = m_shape.dx();
			m_high_dx = m_high_shape.dx();
			//
			set_environment("shape",&m_shape);
			set_environment("dx",&m_dx);
		}
		//
		struct Parameters {
			double speed {5.0};
			double wavelen {0.1};
			double height {0.5};
			unsigned sample_count {8};
			unsigned max_steps {420};
			unsigned upsampling_scale {4};
			bool render_mesh {true};
			bool write_svg {true};
			int contour_count {10};
			int cycle {0};
		};
		//
		static double center_sin( const vec2d &p, double time, const Parameters &param ) {
			const double d = (p-vec2d(0.5,0.5)).len();
			return 0.5*(1+cos(d/param.wavelen+param.speed*time));
		}
		//
		static double linear_sin_x( const vec2d &p, double time, const Parameters &param ) {
			const double d = p[0];
			return 0.5*(1+cos(d/param.wavelen+param.speed*time));
		}
		//
		static double linear_curve_x( const vec2d &p, double time, const Parameters &param ) {
			return 0.5+(p[0]-0.5)*cos(param.speed*time);
		}
		//
		static double linear_curve_xy( const vec2d &p, double time, const Parameters &param ) {
			return 0.5+2.0*(p[1]-0.5)*(p[0]-0.5)*cos(param.speed*time);
		}
		//
		virtual void post_initialize( bool initialized_from_file ) override {
			//
			m_step = 0;
			m_time = 0.0;
			select_function();
			//
			m_grid.clear();
			unsigned n (1);
			for( unsigned k=0; k<4; ++k ) {
				shape2 shape = m_shape/n;
				m_grid.add_layer(shape,n*m_dx);
				n *= 2;
			}
			//
			auto fluid_func = [&](const vec2d &p) { return -1.0; };
			auto adaptivity_func = [&](char depth, const vec2d &p) {
				return 10*(p-vec2d(0.5,0.5)).len() < depth+1;
			};
			//
			m_data.initialize(m_high_shape);
			m_low_data.initialize(m_shape);
			m_grid.activate_cells(adaptivity_func);
			m_grid.balance_layers();
			m_grid.assign_indices();
			m_grid.assign_levelset(fluid_func,nullptr);
			m_cell_values.resize(m_grid.cell_count);
			m_face_values.resize(m_grid.face_count);
			//
			m_camera->set_bounding_box(vec2d().v,m_shape.box(m_dx).v);
			//
			if( console::get_root_path().size()) {
				export_xml();
				m_svg_export_path = console::get_root_path() + "/SVG";
				if( ! filesystem::is_exist(m_svg_export_path)) {
					filesystem::create_directory(m_svg_export_path);
				}
			}
			//
			update();
		}
		//
		void select_function() {
			int selection = m_param.cycle % 4;
			if( selection == 0 ) {
				m_func = center_sin;
			} else if( selection == 1 ) {
				m_func = linear_sin_x;
			} else if( selection == 2 ) {
				m_func = linear_curve_x;
			} else if( selection == 3 ) {
				m_func = linear_curve_xy;
			}
		}
		//
		virtual bool keyboard( int key, int action, int mods ) override {
			//
			if( action == UI_interface::PRESS ) {
				if( key == KEY_C ) {
					m_param.cycle ++;
					select_function();
					update();
					return true;
				} else if( key == KEY_M ) {
					m_mode ++;
					update();
					return true;
				} else if( key == KEY_D ) {
					m_dim ++;
					update();
					return true;
				}
			}
			return drawable::keyboard(key,action,mods);
		}
		//
		virtual void cursor( double x, double y, double z ) override {
			m_mouse_pos = vec2d(x,y);
		}
		//
		void update() {
			//
			if( m_mode % 2 == 0 ) {
				m_grid.iterate_active_cells([&]( const cell_id2 &cell_id, int tid ) {
					m_cell_values[cell_id.index] = m_func(m_grid.get_cell_position(cell_id),m_time,m_param);
				});
				m_data.parallel_all([&]( int i, int j, auto &it ) {
					const vec2d p = m_high_dx*vec2i(i,j).cell();
					it.set(m_grid.sample_cell(p,m_cell_values));
				});
				m_low_data.parallel_all([&]( int i, int j, auto &it ) {
					const vec2d p = m_dx*vec2i(i,j).cell();
					it.set(m_grid.sample_cell(p,m_cell_values));
				});
			} else {
				m_grid.iterate_active_faces([&]( const face_id2 &face_id, int tid ) {
					m_face_values[face_id.index] = m_func(m_grid.get_face_position(face_id),m_time,m_param);
				});
				m_data.parallel_all([&]( int i, int j, auto &it ) {
					const vec2d p = m_high_dx*vec2i(i,j).cell();
					it.set(m_grid.sample_face(p,m_dim%2,m_face_values));
				});
				m_low_data.parallel_all([&]( int i, int j, auto &it ) {
					const vec2d p = m_dx*vec2i(i,j).cell();
					it.set(m_grid.sample_face(p,m_dim%2,m_face_values));
				});
			}
		}
		virtual void idle () override {
			//
			update();
			if( console::get_root_path().size()) {
				export_mesh(m_step);
			}
			//
			if( m_param.write_svg && ! m_svg_export_path.empty()) {
				m_svg_writer->setup_graphics();
				m_svg_writer->clear();
				draw(*m_svg_writer.get());
				std::string final_path = m_svg_export_path + "/" + "output_" + std::to_string(m_step) + ".svg";
				m_svg_writer->const_send_message("write",(char *)final_path.c_str());
				std::string jpg_path = m_svg_export_path + "/" + "output_" + std::to_string(m_step) + ".jpg";
				console::system("convert -flip %s %s",final_path.c_str(),jpg_path.c_str());
				//
				if( m_step && m_step % 10 == 0 ) {
					console::system("rm -rf %s/contour.mp4; avconv -r 60 -i %s/output_%%d.jpg -b:v 12000k -pix_fmt yuv420p %s/contour.mp4",
								m_svg_export_path.c_str(), m_svg_export_path.c_str(), m_svg_export_path.c_str() );
				}
			}
			m_step ++;
			m_time += 0.01;
		}
		//
		virtual bool should_quit() const override { return m_step > m_param.max_steps; }
		//
		void export_xml () {
			//
			std::string export_path = console::get_root_path() + "/mesh";
			if( ! filesystem::is_exist(export_path)) {
				filesystem::create_directory(export_path);
			}
			//
			FILE * file = fopen((export_path+"/mesh.xml").c_str(),"w");
			//
			std::string mesh_xml { 
			"<shape type=\"serialized\">\n"
			"	<string name=\"filename\" value=\"$meshfile\"/>\n"
			"	<boolean name=\"faceNormals\" value=\"true\"/>\n"
			"	<bsdf type=\"twosided\">\n"
			"		<bsdf type=\"diffuse\">\n"
			"			<texture type=\"checkerboard\" name=\"reflectance\">\n"
			"				<float name=\"uvscale\" value=\"6\"/>\n"
			"			</texture>\n"
			"		</bsdf>\n"
			"	</bsdf>\n"
			"</shape>\n"
			};
			fprintf(file,mitsuba_xml::header.c_str());
			fprintf(file,mitsuba_xml::sensor.c_str(),m_param.sample_count);
			fprintf(file,mitsuba_xml::sunlight.c_str());
			fprintf(file,mitsuba_xml::floor.c_str());
			fprintf(file,mesh_xml.c_str());
			fprintf(file,mitsuba_xml::footer.c_str());
			//
			fclose(file);
		}
		//
		void export_mesh( unsigned step ) {
			//
			std::string export_path = console::get_root_path() + "/mesh";
			if( ! filesystem::is_exist(export_path)) {
				filesystem::create_directory(export_path);
			}
			//
			std::vector<vec3d> vertices;
			std::vector<vec2d> uv_coord;
			std::vector<std::vector<size_t> > faces;
			//
			m_high_shape.for_each([&]( int i, int j ) {
				//
				const vec2d p = m_high_dx*vec2i(i,j).cell();
				const double v = m_data(i,j);
				const double value = 0.5*(1.0-m_param.height)+m_param.height*v;
				vertices.push_back(vec3d(p[0],value,p[1]));
				uv_coord.push_back(p);
			});
			(m_high_shape-shape2(1,1)).for_each([&]( int i, int j ) {
				std::vector<size_t> face = {
					m_high_shape.encode(i,j),m_high_shape.encode(i+1,j),m_high_shape.encode(i+1,j+1),m_high_shape.encode(i,j+1)
				};
				faces.push_back(face);
			});
			//
			m_mesh_exporter->set_mesh(vertices,faces);
			m_mesh_exporter->set_texture_coordinates(uv_coord);
			m_mesh_exporter->export_ply(console::format_str("%s/%u_mesh.ply",export_path.c_str(),step));
			m_mesh_exporter->export_mitsuba(console::format_str("%s/%u_mesh.serialized",export_path.c_str(),step));
			//
			if( m_param.render_mesh && console::system("mitsuba > /dev/null 2>&1") == 0 ) {
				//
				std::string img_path = console::get_root_path() + "/img";
				if( ! filesystem::is_exist(img_path)) {
					filesystem::create_directory(img_path);
				}
				for( unsigned step1=0; step1<=step; ++step1 ) {
					if( ! filesystem::is_exist(console::format_str("%s/%u_mesh.exr",img_path.c_str(),step1))) {
						std::string command = console::format_str("mitsuba -Dx=0.5 -Dy=3.0 -Dz=3.5 -Dw=1280 -Dh=720 -Dmeshfile=%s/%u_mesh.serialized -o %s/%u_mesh.exr %s/mesh.xml",
							export_path.c_str(),step1,img_path.c_str(),step1,export_path.c_str());
						//
						console::system(command);
					}
					if( ! filesystem::is_exist(console::format_str("%s/%u_mesh.png",img_path.c_str(),step1))) {
						console::system("mtsutil tonemap -g 3 %s/%u_mesh.exr",img_path.c_str(),step1);
					}
				}
				//
				if( step && step % 10 == 0 ) {
					console::system("rm -f %s/mesh.mp4; avconv -r 60 -i %s/%%d_mesh.png -b:v 12000k -pix_fmt yuv420p %s/mesh.mp4",
								img_path.c_str(), img_path.c_str(), img_path.c_str() );
				}
			}
		}
		//
		virtual void draw( graphics_engine &g ) const override {
			//
			g.color4(0.3,0.3,0.6,1.0);
			g.begin(graphics_engine::MODE::TRIANGLE_FAN);
			g.vertex2(0.0,0.0);
			g.vertex2(1.0,0.0);
			g.vertex2(1.0,1.0);
			g.vertex2(0.0,1.0);
			g.end();
			//
			double theta = 1.0 / m_param.contour_count;
			for( int n=0; n<m_param.contour_count; ++n ) {
				double rgb[4];
				color::heatcolor((n+1)*theta,rgb); rgb[3] = 1.0;
				g.color4v(rgb);
				//
				shared_array2<Real> contour_array(m_low_data);
				contour_array() *= -1.0;
				contour_array() += (n+1)*theta;
				m_gridvisualizer->draw_levelset(g,contour_array(),false);
			}
			//
			if( g.get_graphics_engine_name() == "SVG" ) g.line_width(4.0);
			m_grid.draw_grid(g);
			g.line_width(1.0);
			//
			if( ! m_mouse_pos.empty()) {
				//
				if( m_mode % 2 == 0 ) {
					std::vector<grid2::point_info2> points;
					m_grid.gather_neighbor_cells(m_mouse_pos,points,m_cell_values);
					m_grid.draw_connections(g,m_mouse_pos,points);
				} else {
					std::vector<grid2::point_info2> points;
					m_grid.gather_neighbor_faces(m_mouse_pos,m_dim%2,points,m_face_values);
					m_grid.draw_connections(g,m_mouse_pos,points);
				}
			}
			//
			if( g.get_graphics_engine_name() == "OpenGL" ) {
				if( m_mode % 2 == 0 ) {
					m_grid.draw_ghost_cells(g);
					g.color4(1.0,1.0,1.0,1.0);
					g.draw_string(vec2d(0.01,0.02).v,console::format_str("Error = %e\n",m_func(m_mouse_pos,m_time,m_param)-m_grid.sample_cell(m_mouse_pos,m_cell_values)).c_str());
				} else {
					m_grid.draw_ghost_faces(g,m_dim%2);
					g.color4(1.0,1.0,1.0,1.0);
					g.draw_string(vec2d(0.01,0.02).v,console::format_str("Error = %e\n",m_func(m_mouse_pos,m_time,m_param)-m_grid.sample_face(m_mouse_pos,m_dim%2,m_face_values)).c_str());
				}
			}
		}
		//
		virtual void setup_window( std::string &name, int &width, int &height ) const override {
			double ratio = m_shape[1] / (double) m_shape[0];
			height = ratio * width;
		}
		//
		shape2 m_shape {64,64};
		shape2 m_high_shape;
		unsigned m_step;
		int m_mode {0};
		int m_dim {0};
		double m_dx, m_high_dx;
		double m_time;
		array2<Real> m_data{this,"lineararray2"};
		array2<Real> m_low_data{this,"lineararray2"};
		std::string m_svg_export_path;
		//
		vec2d m_mouse_pos;
		std::vector<Real> m_cell_values;
		std::vector<Real> m_face_values;
		std::function<double( const vec2d &p, double time, const Parameters &param )> m_func {center_sin};
		//
		Parameters m_param;
		//
		grid2 m_grid{this};
		graphics_interface_driver m_svg_writer{this,"graphics_svg"};
		meshexporter3_driver m_mesh_exporter{this,"meshexporter3"};
		gridvisualizer2_driver m_gridvisualizer{this,"gridvisualizer2"};
	};
	//
	extern "C" module * create_instance() {
		return new macoctreeinterpaccuracytest2;
	}
};
//
SHKZ_END_NAMESPACE
//
#endif