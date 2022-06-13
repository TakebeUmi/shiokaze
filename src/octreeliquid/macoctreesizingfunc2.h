/*
**	macoctreesizingfunc2.h
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
#ifndef SHKZ_OCTREESIZINGFUNC2_H
#define SHKZ_OCTREESIZINGFUNC2_H
//
#include "macoctreegrid2.h"
#include <shiokaze/graphics/graphics_engine.h>
//
SHKZ_BEGIN_NAMESPACE
//
namespace macotreeliquid2_namespace {
	//
	class macoctreesizingfunc2 : public recursive_configurable, public credit {
	public:
		//
		LONG_NAME("MAC Octee Sizing Func 2D")
		ARGUMENT_NAME("macoctreesizingfunc")
		//
		macoctreesizingfunc2 ( recursive_configurable *parent ) {
			if( parent ) parent->add_child(this);
			else setup_now();
		}
		//
		virtual void configure( configuration &config ) override;
		//
		void compute_sizing_function( const grid2 &grid0, const grid2 &grid1, double dt,
							 std::function<double( const vec2d &p )> solid_func=nullptr,
							 std::function<vec2d( const vec2d &p )> additional_velocity_func=nullptr );
		//
		void activate_cells( const grid2 &grid0, grid2 &grid1,
							 std::function<double( const vec2d &p )> solid_func=nullptr,
							 std::function<double( const vec2d &p )> additional_fluid_func=nullptr ) const;
		//
		bool keyboard( int key, int action, int mods );
		void cursor( double x, double y, double z );
		void draw( graphics_engine &g, const grid2 &grid ) const;
		//
	protected:
		//
		struct Parameters {
			double curvature_base {8.0};
			double velocity_base {3.0};
			double halfband_width {1.0};
			double sizing_strength {1.0};
			double decay_rate {0.9};
			double decay_time {0.01};
			double velocity_scale {1.0};
			int diffuse_count {5};
			int injection_depth {0};
		};
		Parameters m_param;
		std::vector<Real> m_laplacian;
		vec2d m_mouse_pos;
	};
};
//
SHKZ_END_NAMESPACE
//
#endif
//