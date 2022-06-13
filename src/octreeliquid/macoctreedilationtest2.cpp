/*
**	macoctreedilationtest2.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on December 2, 2019. All rights reserved.
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
#ifndef SHKZ_OCTREESTANDINGWATERTEST2_H
#define SHKZ_OCTREESTANDINGWATERTEST2_H
//
#include <shiokaze/ui/drawable.h>
#include <shiokaze/graphics/graphics_interface.h>
#include <shiokaze/core/filesystem.h>
#include <shiokaze/core/console.h>
#include "macoctreegrid2.h"
#include "macoctreeproject2.h"
#include "macoctreehelper2.h"
//
SHKZ_BEGIN_NAMESPACE
//
namespace macotreeliquid2_namespace {
	//
	class macoctreedilationtest2 : public drawable {
	protected:
		//
		LONG_NAME("MAC Octree Dilation Test 2D")
		ARGUMENT_NAME("macoctreeaccuracytest")
		//
		virtual void configure( configuration &config ) override {
			//
			m_dx = m_shape.dx();
			//
			set_environment("shape",&m_shape);
			set_environment("dx",&m_dx);
			//
			config.set_default_unsigned("CanvasX",640);
			config.set_default_unsigned("CanvasY",640);
			config.set_default_integer("Threads",1);
			config.set_default_bool("SteepAdaptivity",false);
			config.set_default_unsigned("DilateCount",0);
			config.get_string("SVGPath",m_param.svg_export_path,"SVG export path");
		}
		//
		void update() {
			//
			m_grid.param.dilate_count = m_step;
			m_grid.clear();
			unsigned n (1);
			for( unsigned k=0; k<4; ++k ) {
				shape2 shape = m_shape/n;
				m_grid.add_layer(shape,n*m_dx);
				n *= 2;
			}
			auto adaptivity_func = [&](char depth, const vec2d &p) {
				return (p-m_center).len()-m_r < m_dx ? true : false;
			};
			m_grid.activate_cells(adaptivity_func);
			m_grid.balance_layers();
			m_grid.assign_indices();
			//
			auto fluid_func = [&](const vec2d &p) { return (p-m_center).len()-m_r; };
			auto solid_func = [&](const vec2d &p) { return 1.0; };
			//
			m_grid.assign_levelset(fluid_func,solid_func);
			//
			if( ! m_param.svg_export_path.empty()) {
				m_svg_writer->setup_graphics();
				m_svg_writer->clear();
				draw(*m_svg_writer.get());
				std::string final_path = m_param.svg_export_path + "/" + "output_" + std::to_string(m_step) + ".svg";
				m_svg_writer->const_send_message("write",(char *)final_path.c_str());
			}
		}
		//
		virtual void post_initialize( bool initialized_from_file ) override {
			//
			m_step = 0;
			update();
			m_camera->set_bounding_box(vec2d().v,m_shape.box(m_dx).v);
		}
		//
		virtual bool keyboard( int key, int action, int mods ) override {
			if( action == UI_interface::PRESS ) {
				if( key == UI_interface::KEY_C ) {
					m_step ++;
					update();
					return true;
				}
			}
			return false;
		}
		//
		virtual void draw( graphics_engine &g ) const override {
			//
			m_grid.draw_fluid(g);
			m_grid.draw_grid(g);
		}
		//
		virtual void setup_window( std::string &name, int &width, int &height ) const override {
			double ratio = m_shape[1] / (double) m_shape[0];
			height = ratio * width;
		}
		//
		shape2 m_shape {64,64};
		double m_dx;
		unsigned m_step;
		vec2d m_center {0.5,0.5};
		double m_r {0.15};
		//
		struct Parameters {
			std::string svg_export_path;
		};
		Parameters m_param;
		//
		grid2 m_grid{this};
		graphics_interface_driver m_svg_writer{this,"graphics_svg"};
	};
	//
	extern "C" module * create_instance() {
		return new macoctreedilationtest2;
	}
};
//
SHKZ_END_NAMESPACE
//
#endif