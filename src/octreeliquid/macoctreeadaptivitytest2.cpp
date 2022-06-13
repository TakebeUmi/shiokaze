/*
**	macoctreeadaptivitytest2.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on December 31, 2019. All rights reserved.
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
#ifndef SHKZ_OCTREEADAPTIVITYTEST2_H
#define SHKZ_OCTREEADAPTIVITYTEST2_H
//
#include <shiokaze/ui/drawable.h>
#include <shiokaze/graphics/graphics_interface.h>
#include <shiokaze/core/filesystem.h>
#include <shiokaze/core/console.h>
//
#include "macoctreesizingfunc2.h"
#include "macoctreegrid2.h"
//
SHKZ_BEGIN_NAMESPACE
//
namespace macotreeliquid2_namespace {
	//
	class macoctreeadaptivitytest2 : public drawable {
	protected:
		//
		LONG_NAME("MAC Octree Adaptivity Test 2D")
		ARGUMENT_NAME("macadaptivitytest")
		//
		virtual void configure( configuration &config ) override {
			//
			double scale (1.0);
			config.get_double("ResolutionScale",scale,"Resolution doubling scale");
			//
			m_shape *= scale;
			m_dx = m_shape.dx();
			//
			set_environment("shape",&m_shape);
			set_environment("dx",&m_dx);
			//
			config.set_default_unsigned("CanvasX",1280);
			config.set_default_unsigned("CanvasY",640);
			config.set_default_bool("SteepAdaptivity",false);
			//
			config.get_integer("Scene",m_param.scene,"Scene");
			config.get_string("SVGPath",m_param.svg_export_path,"SVG export path");
			config.get_unsigned("MinResolution",m_param.min_resolution,"Minimal resolution");
			config.get_bool("ReuseGrid",m_param.reuse_grid,"Reuse grid");
		}
		//
		virtual void post_initialize( bool initialized_from_file ) override {
			//
			m_step = 0;
			refine([&]() {
				m_grid->activate_cells([&](char depth, const vec2d &p) {
					return depth > 3;
				});
			});
			update(false);
			m_camera->set_bounding_box(vec2d().v,m_shape.box(m_dx).v);
		}
		//
		void refine( std::function<void()> activate_func ) {
			//
			std::swap(m_grid,m_grid_prev);
			m_grid->clear();
			shape2 shape(m_shape);
			double dx (m_dx);
			while( shape.max() >= m_param.min_resolution ) {
				m_grid->add_layer(shape,dx);
				shape = shape/2;
				dx *= 2.0;
			}
			activate_func();
			m_grid->balance_layers();
			m_grid->assign_indices();
		}
		//
		void update( bool use_prev_grid ) {
			//
			auto solid_func = [&](const vec2d &p) { return 1.0; };
			//
			if( use_prev_grid ) {
				m_grid->assign_levelset([&](const vec2d &p) {
					return m_grid_prev->sample_levelset(p);
				},solid_func);
			} else {
				if( m_param.scene == 0 ) {
					m_grid->assign_levelset([&](const vec2d &p) {
						double value = (p-vec2d(0.5,0.5)).len()-0.2;
						value = std::min(value,p[1]-0.5);
						return value;
					},solid_func);
				} else if( m_param.scene == 1 ) {
					m_grid->assign_levelset([&](const vec2d &p) {
						double value = (p-vec2d(0.5,0.7)).len()-0.2;
						value = std::min(value,p[1]-0.5);
						return value;
					},solid_func);
				}
			}
			//
			if( m_param.scene == 0 ) {
				m_grid->set_velocity([&]( const vec2d &p, char dim ) { return 0.0; });
			} else if( m_param.scene == 1 ) {
				m_grid->set_velocity([&]( const vec2d &p, char dim ) { return p[1]<0.5 ? 0.0 : vec2d(0.0,-1.0)[dim]; });
			}
			//
			m_macoctreesizingfunc.compute_sizing_function(*m_grid_prev,*m_grid,1.0);
			m_step ++;
		}
		//
		virtual bool keyboard( int key, int action, int mods ) override {
			//
			if( ! m_macoctreesizingfunc.keyboard(key,action,mods)) {
				if( action == UI_interface::PRESS && key == UI_interface::KEY_Q ) {
					//
					refine([&]() {
						m_macoctreesizingfunc.activate_cells(*m_grid_prev,*m_grid);
					});
					update(m_param.reuse_grid);
					//
					if( ! m_param.svg_export_path.empty()) {
						m_svg_writer->setup_graphics();
						m_svg_writer->clear();
						draw(*m_svg_writer.get());
						std::string final_path = m_param.svg_export_path + "/" + "output_" + std::to_string(m_step) + ".svg";
						m_svg_writer->const_send_message("write",(char *)final_path.c_str());
					}
					return true;
				}
			}
			return false;
		}
		//
		virtual void draw( graphics_engine &g ) const override {
			//
			m_grid->draw_fluid(g);
			m_grid->draw_grid(g,false);
			m_macoctreesizingfunc.draw(g,*m_grid);
		}
		//
		virtual void setup_window( std::string &name, int &width, int &height ) const override {
			double ratio = m_shape[1] / (double) m_shape[0];
			height = ratio * width;
		}
		//
		virtual void cursor( double x, double y, double z ) override {
			m_macoctreesizingfunc.cursor(x,y,z);
		}
		//
		shape2 m_shape {64,64};
		double m_dx;
		unsigned m_step;
		//
		struct Parameters {
			std::string svg_export_path;
			int scene {0};
			bool reuse_grid {false};
			unsigned min_resolution {8};
		};
		Parameters m_param;
		//
		grid2 m_grid0{this};
		grid2 m_grid1{this};
		//
		grid2 *m_grid {&m_grid0};
		grid2 *m_grid_prev {&m_grid1};
		//
		macoctreesizingfunc2 m_macoctreesizingfunc{this};
		graphics_interface_driver m_svg_writer{this,"graphics_svg"};
		//
		std::vector<Real> lap_levelset;
	};
	//
	extern "C" module * create_instance() {
		return new macoctreeadaptivitytest2;
	}
};
//
SHKZ_END_NAMESPACE
//
#endif