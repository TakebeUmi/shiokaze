/*
**	macoctreemesher3.h
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on November 11, 2019. All rights reserved.
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
#ifndef SHKZ_OCTREEMESHER3_H
#define SHKZ_OCTREEMESHER3_H
//
#include "macoctreegrid3.h"
//
SHKZ_BEGIN_NAMESPACE
//
namespace macotreeliquid3_namespace {
	//
	class macoctreemesher3 : public recursive_configurable, public credit {
	public:
		//
		LONG_NAME("MAC Octee Mesher 3D")
		ARGUMENT_NAME("macoctreemesher")
		//
		macoctreemesher3 ( recursive_configurable *parent ) {
			if( parent ) parent->add_child(this);
			else setup_now();
		}
		//
		virtual void configure( configuration &config ) override;
		//
		struct Parameters {
			bool snap_vertices {true};
			bool displace_solid_embedded_vertices {true};
			double displace_solid_embedded_rate {0.8};
			bool project_vertices {true};
			double mesh_proximity {0.25};
			bool remove_face {true};
		};
		//
		Parameters m_param;
		parallel_driver m_parallel{this};
		//
		void generate_mesh( const grid3 &grid, const std::vector<Real> &cell_values, double offset, std::function<double(const vec3d &p)> solid_func, std::vector<vec3d> &vertices, std::vector<std::vector<size_t> > &faces, bool enclose ) const;
	};
};
//
SHKZ_END_NAMESPACE
//
#endif
//