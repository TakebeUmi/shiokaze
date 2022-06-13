/*
**	macoctreeaccuracytest2.cpp
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
#ifndef SHKZ_OCTREEACCURACYTEST2_H
#define SHKZ_OCTREEACCURACYTEST2_H
//
#include "macoctreegrid2.h"
#include "macoctreeproject2.h"
#include <shiokaze/ui/drawable.h>
#include <shiokaze/graphics/graphics_interface.h>
#include <shiokaze/core/filesystem.h>
#include <shiokaze/core/console.h>
#include <shiokaze/image/color.h>
#include <cmath>
#include <mutex>
#include <typeinfo>
//
SHKZ_BEGIN_NAMESPACE
//
namespace macotreeliquid2_namespace {
	//
	class macoctreeaccuracytest2 : public drawable {
	protected:
		//
		LONG_NAME("MAC Octree Accuracy Test 2D")
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
			config.set_default_bool("DrawVelocity",true);
			config.set_default_unsigned("DilateCount",0);
			config.set_default_unsigned("MaxIterations",300000);
			config.set_default_bool("RemoveOneDegreesOfFreedom",true);
			//
			config.get_double("Radius",m_r,"Circle radius");
			config.get_string("SVGPath",m_svg_export_path,"SVG export path");
			config.get_integer("TrialCount",m_max_trial_count,"Number of trial count");
			config.get_integer("SubdivisionCount",m_max_subdivion_count,"Grid subdivision count");
			config.get_vec2d("Shift",m_shift.v,"Shift amount");
			config.get_bool("UniformAdaptivity",m_uniform_adaptivity,"Use uniform grid");
		}
		//
		double analytical_function( const vec2d &center, const vec2d &p, double r ) {
			//
			if( m_case == 0 ) {
				double d2 = (p-center).norm2();
				return d2-r*r;
			} else if( m_case == 1 ) {
				const vec2d q = p-center;
				const double x = 2.0*q[0];
				const double y = 2.0*q[1];
				return x*y*y*y;
			}
			return 0.0;
		}
		//
		vec2d derivative_analytical_function( const vec2d &center, const vec2d &p ) {
			//
			const vec2d r = p-center;
			if( m_case == 0 ) {
				return 2.0*r;
			} else if( m_case == 1 ) {
				const double x = 2.0*r[0];
				const double y = 2.0*r[1];
				return vec2d(y+y*y*y,-x+3.0*x*y*y);
			}
			return vec2d();
		}
		//
		double measure_surface_norm ( const vec2d &center, double r, unsigned depth ) {
			//
			double inf_norm (0.0);
			if( m_case == 0 ) {
				m_grid.serial_iterate_active_cells([&]( const cell_id2 &cell_id ) {
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
				auto check_valid_cell = [&]( const cell_id2 &cell_id ) {
					bool valid (false);
					for( char dim : DIMS2 ) {
						m_grid.iterate_face_neighbors(cell_id,dim,[&]( const face_id2 &face_id ) {
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
				m_grid.serial_iterate_active_cells([&]( const cell_id2 &cell_id ) {
					if( cell_id.depth == depth && check_valid_cell(cell_id)) {
						const vec2d &p = m_grid.get_cell_position(cell_id);
						analytical_average += analytical_function(center,p,r);
						discrete_average += m_pressure[cell_id.index];
						sample_sum ++;
					}
				});
				analytical_average /= (double)sample_sum;
				discrete_average /= (double)sample_sum;
				//
				m_grid.serial_iterate_active_cells([&]( const cell_id2 &cell_id ) {
					if( cell_id.depth == depth && check_valid_cell(cell_id)) {
						const vec2d &p = m_grid.get_cell_position(cell_id);
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
		virtual void post_initialize( bool initialized_from_file ) override {
			//
			if( typeid(Real) == typeid(double)) {
				console::dump( "Real = double\n");
			} else if( typeid(Real) == typeid(float)) {
				console::dump( "Real = float\n");
			}
			//
			if( console::get_root_path().size()) {
				if( m_svg_export_path.empty()) {
					m_svg_export_path = console::get_root_path() + "/SVG";
					if( ! filesystem::is_exist(m_svg_export_path)) {
						filesystem::create_directory(m_svg_export_path);
					}
				}
			}
			//
			m_grid.clear();
			console::set_time(m_shape[0]);
			//
			unsigned n (1);
			for( unsigned k=0; k<4; ++k ) {
				shape2 shape = m_shape/n;
				m_grid.add_layer(shape,n*m_dx);
				n *= 2;
			}
			//
			const vec2d center = vec2d(0.5,0.5)+m_shift;
			auto setup = [&]( double r ) {
				//
				auto fluid_func = [&](const vec2d &p) { 
					if( m_case == 0 ) {
						return (p-center).len()-r;
					} else {
						return -1.0;
					}
				};
				auto solid_func = [&](const vec2d &p) {
					if( m_case == 0 ) {
						return 1.0;
					} else if( m_case == 1 ) {
						return r-(p-center).len();
					}
					return 1.0;
				};
				auto adaptivity_func = [&](char depth, const vec2d &p) {
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
				m_grid.set_velocity([&]( const vec2d &p, char dim ) { return derivative_analytical_function(center,p)[dim]; });
				m_grid.serial_iterate_active_faces([&]( const face_id2 &face_id ) {
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
				const double t = 2.0*sqrt(2.0)*m_dx/m_max_trial_count;
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
			if( ! m_svg_export_path.empty()) {
				unsigned prefix_num = m_shape[0];
				m_svg_writer->setup_graphics();
				m_svg_writer->clear();
				draw(*m_svg_writer.get());
				std::string final_path = m_svg_export_path + "/" + "output_" + std::to_string(prefix_num) + ".svg";
				m_svg_writer->const_send_message("write",(char *)final_path.c_str());
			}
			m_camera->set_bounding_box(vec2d().v,m_shape.box(m_dx).v);
		}
		//
		virtual bool keyboard( int key, int action, int mods ) override {
			//
			if( action == UI_interface::PRESS && key == KEY_C ) {
				//
				m_dx *= 0.5;
				m_shape *= 2.0;
				m_doubled_count ++;
				reinitialize();
				return true;
			}
			return drawable::keyboard(key,action,mods);
		}
		//
		virtual void idle () override {
			//
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
		virtual void draw( graphics_engine &g ) const override {
			//
			m_grid.draw_fluid(g);
			//
			if( m_case == 0 ) {
				int max_count (6);
				double theta = 1.0 / max_count;
				for( int n=0; n<max_count; ++n ) {
					double rgb[4];
					color::heatcolor((n+1)*theta,rgb); rgb[3] = 1.0;
					m_grid.draw_contour(g,m_pressure,0.025*(n+1),rgb,false);
				}
			}
			m_grid.draw_grid(g);
		}
		//
		virtual void setup_window( std::string &name, int &width, int &height ) const override {
			double ratio = m_shape[1] / (double) m_shape[0];
			height = ratio * width;
		}
		//
		struct NormPoint2 {
			vec2d p;
			double dx;
		};
		//
		shape2 m_shape {64,64};
		double m_dx;
		double m_r {0.4};
		vec2d m_shift;
		int m_doubled_count {0};
		unsigned m_step {0};
		int m_max_trial_count {4};
		int m_max_subdivion_count {4};
		int m_case {0};
		bool m_uniform_adaptivity {false};
		std::vector<Real> m_pressure;
		std::string m_svg_export_path;
		//
		grid2 m_grid{this};
		macoctreeproject2 m_macoctreeproject{this};
		graphics_interface_driver m_svg_writer{this,"graphics_svg"};
	};
	//
	extern "C" module * create_instance() {
		return new macoctreeaccuracytest2;
	}
};
//
SHKZ_END_NAMESPACE
//
#endif