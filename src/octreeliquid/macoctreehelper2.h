/*
**	macoctreehelper2.h
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on November 12, 2019. All rights reserved.
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
#ifndef SHKZ_OCTREEHELPER2_H
#define SHKZ_OCTREEHELPER2_H
//
#include <shiokaze/graphics/graphics_engine.h>
#include "macoctreegrid2.h"
#include "macoctreeproject2.h"
//
SHKZ_BEGIN_NAMESPACE
//
namespace macotreeliquid2_namespace {
	//
	class macoctreehelper2 {
	public:
		//
		void initialize();
		void reset_focus();
		void cursor( const grid2 &grid, double x, double y, double z );
		bool keyboard( int key, int action, int mods );
		void draw_debug( const grid2 &grid, const macoctreeproject2::matrix2 &matrix, graphics_engine &g ) const;
		//
	protected:
		//
		vec2d m_mouse_pos;
		cell_id2 m_focus_cell;
		face_id2 m_focus_face;
		bool m_draw_points;
		std::string m_mode;
	};
};
//
SHKZ_END_NAMESPACE
//
#endif
//