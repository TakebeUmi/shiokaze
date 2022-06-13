/*
**	macoctreeaccuracytest2.cpp
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
	class macoctreeaccuracytest2 : public drawable {
	protected:
		//
		LONG_NAME("MAC Octree Standing Water Test 2D")
		ARGUMENT_NAME("macoctreeaccuracytest")
		//
		virtual void configure( configuration &config ) override {
			//
			m_dx = m_shape.dx();
			//
			set_environment("shape",&m_shape);
			set_environment("dx",&m_dx);
			//
			config.set_default_unsigned("CanvasX",1280);
			config.set_default_unsigned("CanvasY",640);
			config.set_default_integer("Threads",1);
			config.set_default_double("Residual",1e-18);
			config.set_default_double("EpsFluid",0.0);
			config.set_default_double("AccuracyEps",0.0);
			config.set_default_bool("SteepAdaptivity",false);
			config.set_default_bool("DegradeWarning",true);
			config.set_default_bool("VolumeCorrection",false);
			config.set_default_unsigned("DilateCount",0);
			config.set_default_unsigned("MaxIterations",300000);
			config.get_string("SVGPath",m_param.svg_export_path,"SVG export path");
			//
			config.get_bool("Debug",m_param.debug,"Debug mode");
			config.get_bool("OnlyVertical",m_param.only_vertical,"Allow only vertical movement");
			m_param.debug = m_param.debug || config.is_parameter_set("LinSolver");
			config.get_double("Level",m_param.water_level,"Water level");
			config.get_double("Scale",m_param.scale,"Force scale");
		}
		//
		virtual void post_initialize( bool initialized_from_file ) override {
			//
			m_step = 0;
			m_time = 0.0;
			m_grid.clear();
			m_macoctreehelper.initialize();
			unsigned n (1);
			for( unsigned k=0; k<4; ++k ) {
				shape2 shape = m_shape/n;
				m_grid.add_layer(shape,n*m_dx);
				n *= 2;
			}
			auto adaptivity_func = [&](char depth, const vec2d &p) {
				return 4*p[0] < depth+1;
			};
			m_grid.activate_cells(adaptivity_func);
			m_grid.balance_layers();
			m_grid.assign_indices();
			//
			idle();
			m_camera->set_bounding_box(vec2d().v,m_shape.box(m_dx).v);
		}
		//
		virtual void idle() override {
			//
			m_time += 0.01;
			const vec2d &tilt_vec = vec2d(m_param.only_vertical ? 0.0 : 0.1*cos(3.0*m_time),1.0).normal();
			//
			auto fluid_func = [&](const vec2d &p) { return (tilt_vec*p)-m_param.water_level+0.05*sin(m_time); };
			auto solid_func = [&](const vec2d &p) { return 1.0; };
			//
			m_grid.assign_levelset(fluid_func,solid_func);
			m_grid.set_velocity([&]( const vec2d &p, char dim ) { return -m_param.scale*tilt_vec[dim]; });
			if( m_param.debug ) {
				m_macoctreeproject.assemble_matrix(m_grid);
				m_macoctreehelper.reset_focus();
			}
			m_macoctreeproject.project(m_grid,1.0);
			//
			if( ! m_param.svg_export_path.empty()) {
				m_svg_writer->setup_graphics();
				m_svg_writer->clear();
				draw(*m_svg_writer.get());
				std::string final_path = m_param.svg_export_path + "/" + "output_" + std::to_string(m_step) + ".svg";
				m_svg_writer->const_send_message("write",(char *)final_path.c_str());
				//
				if( m_step == 1000 ) exit(0);
			}
			//
			m_step ++;
		}
		//
		virtual void cursor( double x, double y, double z ) override {
			m_macoctreehelper.cursor(m_grid,x,y,z);
		}
		//
		virtual bool keyboard( int key, int action, int mods ) override {
			return m_macoctreehelper.keyboard(key,action,mods);
		}
		//
		virtual void draw( graphics_engine &g ) const override {
			//
			m_grid.draw_fluid(g);
			m_grid.draw_grid(g);
			//
			if( m_param.debug ) {
				m_macoctreehelper.draw_debug(m_grid,m_macoctreeproject.m_matrix,g);
			}
		}
		//
		virtual void setup_window( std::string &name, int &width, int &height ) const override {
			double ratio = m_shape[1] / (double) m_shape[0];
			height = ratio * width;
		}
		//
		shape2 m_shape {64,32};
		double m_dx;
		double m_time;
		unsigned m_step;
		//
		struct Parameters {
			bool debug {true};
			bool only_vertical {false};
			double water_level {0.25};
			double scale {5e2};
			std::string svg_export_path;
		};
		Parameters m_param;
		//
		grid2 m_grid{this};
		macoctreeproject2 m_macoctreeproject{this};
		macoctreehelper2 m_macoctreehelper;
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