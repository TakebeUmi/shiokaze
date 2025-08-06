/*
**	macoctreesegregtator2.h
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on December 8, 2019. All rights reserved.
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
#ifndef SHKZ_OCTREESEGREGATOR2_H
#define SHKZ_OCTREESEGREGATOR2_H
//
#include "macoctreegrid2.h"
#include <shiokaze/graphics/graphics_engine.h>
//
SHKZ_BEGIN_NAMESPACE
//
namespace macotreeliquid2_namespace {
	//
	class macoctreesegregator2 : public recursive_configurable, public credit {
	public:
		//
		LONG_NAME("MAC Octee Segregator 2D")
		ARGUMENT_NAME("macoctreesegregator")
		//
		macoctreesegregator2 ( recursive_configurable *parent ) {
			if( parent ) parent->add_child(this);
			else setup_now();
		}
		//
		virtual void configure( configuration &config ) override;
		//
		size_t segregate( const grid2 &grid, std::vector<uint_type> &regions ) const;
		void extrapolate( const grid2 &grid, std::vector<uint_type> &regions ) const;
		void extrapolate_jacobi( const grid2 &grid, std::vector<uint_type> &regions ) const;
		void backtrace( const grid2 &grid0, const grid2 &grid1,
						double dt, std::function<vec2d( const vec2d &p )> velocity,
						const std::vector<uint_type> &regions0, std::vector<uint_type> &regions1 ) const;
		void prune( const grid2 &grid, std::vector<uint_type> &regions ) const;
		void compute( const grid2 &grid, const std::vector<uint_type> &regions, std::vector<Real> &volumes ) const;
		double compute_target_volume( const grid2 &grid1,
									const std::vector<uint_type> &regions01, const std::vector<Real> &volumes0,
									const std::vector<uint_type> &regions11, const std::vector<Real> &volumes1,
									const std::vector<Real> &y_lists0,
									std::vector<Real> &current_volumes1,
									std::vector<Real> &target_volumes1,
									std::vector<Real> &new_y_list1 ) const;
		//
		void draw_region( graphics_engine &g, const grid2 &grid, const std::vector<uint_type> &regions ) const;
		//
	protected:
		//
		struct Parameters {
			double ratio_limit {0.0};
			bool form_matrix {true};
		};
		Parameters m_param;
	};
};
//
SHKZ_END_NAMESPACE
//
#endif
//