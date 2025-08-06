/*
**	macoctreehelper2.cpp
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
#include "macoctreehelper2.h"
#include <shiokaze/ui/UI_interface.h>
#include <shiokaze/graphics/graphics_utility.h>
#include <mutex>
#include <thread>
//
SHKZ_USING_NAMESPACE
using namespace macotreeliquid2_namespace;
//
void macoctreehelper2::initialize() {
	//
	m_mouse_pos = vec2d();
	reset_focus();
	//
	m_mode = "gradient";
	m_draw_points = false;
}
//
void macoctreehelper2::reset_focus() {
	m_focus_cell.depth = -1;
	m_focus_face.depth = -1;
}
//
void macoctreehelper2::cursor( const grid2 &grid, double x, double y, double z ) {
	//
	// Record the mouse position
	m_mouse_pos = vec2d(x,y);
	//
	// Pick the one with the closet degrees of freedom
	double min_dist (1e9);
	//
	// Reset focus
	reset_focus();
	//
	std::mutex mutex;
	if( m_mode == "gradient" || m_mode == "face" ) {
		grid.serial_iterate_active_faces([&]( const face_id2 &face_id ) {
			std::lock_guard<std::mutex> guard(mutex);
			const vec2d &pos = grid.get_face_position(face_id);
			double dist = (pos-m_mouse_pos).len();
			if( dist < min_dist ) {
				m_focus_face = face_id;
				min_dist = dist;
			}
		});
	} else if( m_mode == "divergence" || m_mode == "laplacian" || m_mode == "cell" ) {
		grid.serial_iterate_active_cells([&]( const cell_id2 &cell_id ) {
			std::lock_guard<std::mutex> guard(mutex);
			const vec2d &pos = grid.get_cell_position(cell_id);
			double dist = (pos-m_mouse_pos).len();
			if( dist < min_dist ) {
				m_focus_cell = cell_id;
				min_dist = dist;
			}
		});
	}
}
//
bool macoctreehelper2::keyboard( int key, int action, int mods ) {
	//
	auto reset_focus = [&]() {
		m_focus_cell.depth = -1;
		m_focus_face.depth = -1;
	};
	//
	if( action == UI_interface::PRESS ) {
		if( key == UI_interface::KEY_G ) {
			reset_focus();
			m_mode = "gradient";
			return true;
		} else if( key == UI_interface::KEY_F ) {
			reset_focus();
			m_mode = "face";
			return true;
		} else if( key == UI_interface::KEY_C ) {
			reset_focus();
			m_mode = "cell";
			return true;
		} else if( key == UI_interface::KEY_D) {
			reset_focus();
			m_mode = "divergence";
			return true;
		} else if( key == UI_interface::KEY_L ) {
			reset_focus();
			m_mode = "laplacian";
			return true;
		} else if( key == UI_interface::KEY_P ) {
			m_draw_points = ! m_draw_points;
			return true;
		}
	}
	return false;
}
//
void macoctreehelper2::draw_debug( const grid2 &grid, const macoctreeproject2::matrix2 &matrix, graphics_engine &g ) const {
	//
	if( m_draw_points ) {
		//
		double r = 0.15 * grid.layers.front()->dx;
		//
		std::vector<vec2d> cell_positions(grid.cell_count);
		grid.serial_iterate_active_cells([&]( const cell_id2 &cell_id ) {
			if( grid.cell_map[cell_id.index] ) {
				uint_type row = grid.cell_map[cell_id.index]-1;
				cell_positions[row] = grid.get_cell_position(cell_id);
			}
		});
		//
		std::vector<vec2d> face_positions(grid.face_count);
		grid.serial_iterate_active_faces([&]( const face_id2 &face_id ) {
			if( grid.face_map[face_id.index] ) {
				uint_type row = grid.face_map[face_id.index]-1;
				face_positions[row] = grid.get_face_position(face_id);
			}
		});
		//
		g.color4(1.0,1.0,0.0,0.5);
		grid.serial_iterate_active_cells([&]( const cell_id2 &cell_id ) {
			if( grid.cell_map[cell_id.index] ) {
				graphics_utility::draw_circle(g,grid.get_cell_position(cell_id).v,r,graphics_engine::MODE::TRIANGLE_FAN);
			}
		});
		//
		g.color4(1.0,1.0,1.0,0.5);
		grid.serial_iterate_active_faces([&]( const face_id2 &face_id ) {
			if( grid.face_map[face_id.index] ) {
				graphics_utility::draw_circle(g,grid.get_face_position(face_id).v,r,graphics_engine::MODE::TRIANGLE_FAN);
			}
		});
		//
		if( m_mode == "cell" ) {
			//
			if( m_focus_cell.depth >= 0 ) {
				vec2d cell_position = grid.get_cell_position(m_focus_cell);
				graphics_utility::draw_circle(g,cell_position.v,r,graphics_engine::MODE::TRIANGLE_FAN);
				grid.iterate_cell_neighbors(m_focus_cell,[&]( char dim, const cell_id2 &cell_id ) {
					//
					const vec2d &cell_p = grid.get_cell_position(cell_id);
					//
					g.begin(graphics_engine::MODE::LINES);
					g.vertex2v(cell_position.v);
					g.vertex2v(cell_p.v);
					g.end();
					//
					g.draw_string((cell_p+vec2d(r,r)).v,console::format_str("%u",cell_id.index).c_str());
					graphics_utility::draw_circle(g,cell_p.v,r,graphics_engine::MODE::TRIANGLE_FAN);
				});
			}
		} else if( m_mode == "face" ) {
			//
			if( m_focus_face.depth >= 0 ) {
				vec2d face_position = grid.get_face_position(m_focus_face);
				graphics_utility::draw_circle(g,face_position.v,r,graphics_engine::MODE::TRIANGLE_FAN);
				grid.iterate_face_neighbors(m_focus_face,[&]( const face_id2 &face_id ) {
					//
					const vec2d &face_p = grid.get_face_position(face_id);
					//
					g.begin(graphics_engine::MODE::LINES);
					g.vertex2v(face_position.v);
					g.vertex2v(face_p.v);
					g.end();
					//
					g.draw_string((face_p+vec2d(r,r)).v,console::format_str("%u",face_id.index).c_str());
					graphics_utility::draw_circle(g,face_p.v,r,graphics_engine::MODE::TRIANGLE_FAN);
				});
			}
		} else if( m_mode == "gradient" ) {
			//
			if( m_focus_face.depth >= 0 && ! matrix.G->empty()) {
				vec2d face_position = grid.get_face_position(m_focus_face);
				if( grid.face_map[m_focus_face.index] ) {
					uint_type row = grid.face_map[m_focus_face.index]-1;
					matrix.G->const_for_each(row,[&]( size_t column, double value ) {
						//
						const vec2d &cell_p = cell_positions[column];
						//
						g.begin(graphics_engine::MODE::LINES);
						g.vertex2v(face_position.v);
						g.vertex2v(cell_p.v);
						g.end();
						//
						g.draw_string((cell_p+vec2d(r,r)).v,console::format_str("%.3f",value).c_str());
						graphics_utility::draw_circle(g,cell_p.v,r,graphics_engine::MODE::TRIANGLE_FAN);
					});
				}
			}
		} else if( m_mode == "divergence" ) {
				//
				if( m_focus_cell.depth >= 0 && ! matrix.D->empty()) {
					if( grid.cell_map[m_focus_cell.index] ) {
						uint_type row = grid.cell_map[m_focus_cell.index]-1;
						vec2d pressure_position = grid.get_cell_position(m_focus_cell);
						graphics_utility::draw_circle(g,pressure_position.v,r,graphics_engine::MODE::TRIANGLE_FAN);
						double sum (0.0);
						matrix.D->const_for_each(row,[&]( size_t column, double value ) {
							//
							sum += value;
							const vec2d &face_p = face_positions[column];
							//
							g.begin(graphics_engine::MODE::LINES);
							g.vertex2v(pressure_position.v);
							g.vertex2v(face_p.v);
							g.end();
							//
							g.draw_string((face_p+vec2d(r,r)).v,console::format_str("%.3f",value).c_str());
							graphics_utility::draw_circle(g,face_p.v,r,graphics_engine::MODE::TRIANGLE_FAN);
						});
						//
						g.draw_string(vec2d(0.05,0.05).v,console::format_str("Sum = %.3e",sum).c_str());
					}
				}
		} else if( m_mode == "laplacian" ) {
			//
			if( m_focus_cell.depth >= 0 ) {
				if( grid.cell_map[m_focus_cell.index] ) {
					uint_type row = grid.cell_map[m_focus_cell.index]-1;
					vec2d pressure_position = grid.get_cell_position(m_focus_cell);
					graphics_utility::draw_circle(g,pressure_position.v,r,graphics_engine::MODE::TRIANGLE_FAN);
					double sum (0.0);
					matrix.Lhs->const_for_each(row,[&]( size_t column, double value ) {
						//
						sum += value;
						const vec2d &cell_p = cell_positions[column];
						//
						g.begin(graphics_engine::MODE::LINES);
						g.vertex2v(pressure_position.v);
						g.vertex2v(cell_p.v);
						g.end();
						//
						g.draw_string((cell_p+vec2d(r,r)).v,console::format_str("%.3f",value).c_str());
						graphics_utility::draw_circle(g,cell_p.v,r,graphics_engine::MODE::TRIANGLE_FAN);
					});
					//
					g.draw_string(vec2d(0.05,0.05).v,console::format_str("Sum = %.3e",sum).c_str());
				}
			}
		}
	}
}
//