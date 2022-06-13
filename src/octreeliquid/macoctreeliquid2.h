/*
**	macoctreeliquid2.h
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on October 9, 2018.All rights reserved.
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
#ifndef SHKZ_OCTREELIQUID2_H
#define SHKZ_OCTREELIQUID2_H
//
#include <shiokaze/ui/drawable.h>
#include <shiokaze/timestepper/timestepper_interface.h>
#include <shiokaze/visualizer/gridvisualizer2_interface.h>
#include <shiokaze/utility/macutility2_interface.h>
#include <shiokaze/utility/gridutility2_interface.h>
#include <shiokaze/utility/meshutility2_interface.h>
#include <shiokaze/flip/macflip2_interface.h>
#include <shiokaze/array/array2.h>
#include <shiokaze/array/macarray2.h>
#include <shiokaze/core/dylibloader.h>
#include <vector>
#include <memory>
#include <utility>
//
#include "macoctreegrid2.h"
#include "macoctreeproject2.h"
#include "macoctreehelper2.h"
#include "macoctreesegregator2.h"
#include "macoctreesizingfunc2.h"
//
SHKZ_BEGIN_NAMESPACE
//
namespace macotreeliquid2_namespace {
	//
	class macoctreeliquid2 : public drawable {
	public:
		//
		LONG_NAME("MAC Octree Liquid 2D")
		ARGUMENT_NAME("macoctreeliquid")
		//
	protected:
		//
		shape2 m_shape;
		double m_dx;
		double m_narrowband_depth;
		//
		gridvisualizer2_driver m_solid_gridvisualizer{this,"gridvisualizer2"};
		array2<Real> m_solid;
		//
		grid2 m_grid_0{this};
		grid2 m_grid_1{this};
		//
		grid2 *m_grid {&m_grid_0};
		grid2 *m_grid_prev {&m_grid_1};
		//
		// -------- Reginal Volume Correction Members --------
		std::vector<Real> m_y_list;
		uint_type m_region_count;
		std::vector<uint_type> m_regions;
		std::vector<Real> m_current_volumes;
		std::vector<Real> m_volumes;
		// ---------------------------------------------------
		//
		macflip2_driver m_flip{this,"macexnbflip2"};
		//
		virtual void load( configuration &config ) override;
		virtual void configure( configuration &config ) override;
		virtual void post_initialize( bool initialized_from_file ) override;
		//
		virtual void cursor( double x, double y, double z ) override;
		virtual bool keyboard( int key, int action, int mods ) override;
		virtual void drag( double x, double y, double z, double u, double v, double w ) override;
		virtual void idle() override;
		virtual void setup_window( std::string &name, int &width, int &height ) const override;
		virtual void draw( graphics_engine &g ) const override;
		virtual bool should_quit() const override { return m_timestepper->should_quit(); }
		virtual bool should_screenshot() const override { return m_timestepper->should_export_frame(); }
		//
		timestepper_driver m_timestepper{this,"timestepper"};
		parallel_driver m_parallel{this};
		dylibloader m_dylib;
		//
		std::function<double(const vec2d &)> m_solid_func, m_combined_solid_func;
		std::function<void(double,Real [DIM2][2])> m_set_boundary_flux;
		std::function<vec2d(double)> m_gravity_func;
		std::function<void(graphics_engine &,double)> m_draw_func;
		std::function<std::pair<double,vec2d>( double, const vec2d &)> m_moving_solid_func;
		//
		macoctreeproject2 m_macoctreeproject{this};
		macoctreesegregator2 m_macoctreesegregator{this};
		macoctreehelper2 m_macoctreehelper;
		meshutility2_driver m_meshutility{this,"meshutility2"};
		macoctreesizingfunc2 m_macoctreesizingfunc{this};
		//
		environment_setter arg_shape{this,"shape",&m_shape};
		environment_setter arg_dx{this,"dx",&m_dx};
		//
		struct Parameters {
			//
			unsigned min_resolution {16};
			bool use_FLIP {true};
			vec2d gravity{0.0,-9.8};
			double PICFLIP {0.98};
			bool volume_correction {true};
			bool regional_volume_correction {false};
			bool maccormack {false};
			unsigned erode_width {0};
			double surftens_k {0.0};
			bool use_sizing_func {true};
			unsigned initial_refinement {3};
			double maximal_CFL_accumulation {1.0};
		};
		//
		Parameters m_param;
		//
		std::function<bool( double dx, double dt, double time, unsigned step )> m_check_inject_func;
		std::function<bool( const vec2d &p, double dx, double dt, double time, unsigned step, double &fluid, vec2d &velocity )> m_inject_func;
		std::function<void( double dx, double dt, double time, unsigned step, double &volume_change )> m_post_inject_func;
		//
		bool m_do_inject;
		double m_injected_volume;
		double m_accumulated_CFL;
		//
		virtual void begin_inject_external_fluid( double dt, double time, unsigned step );
		virtual void do_inject_external_fluid( double dt, double time, unsigned step );
		virtual void end_inject_external_fluid( double dt, double time, unsigned step );
	};
//
};
//
SHKZ_END_NAMESPACE
//
#endif
//