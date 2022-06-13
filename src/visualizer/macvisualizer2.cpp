/*
**	macvisualizer2.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on July 21, 2017.
**
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
#include <shiokaze/array/shared_array2.h>
#include <shiokaze/visualizer/macvisualizer2_interface.h>
#include <shiokaze/graphics/graphics_utility.h>
#include <algorithm>
#include <limits>
//
SHKZ_USING_NAMESPACE
//
class macvisualizer2 : public macvisualizer2_interface {
protected:
	//
	virtual void draw_velocity( graphics_engine &g, const macarray2<Real> &velocity ) const override {
		//
		if( m_param.draw_velocity ) {
			//
			if( m_param.draw_mac ) {
				//
				velocity.const_serial_actives([&]( int dim, int i, int j, const auto &it ) {
					vec2d p0 = m_dx*vec2i(i,j).face(dim);
					vec2d p1 = p0+m_dx*it()*vec2d(dim==0,dim==1);
					if( dim == 0 ) {
						g.color4(0.75,0.75,1.0,0.5);
					} else {
						g.color4(1.0,0.75,0.75,0.5);
					}
					graphics_utility::draw_arrow(g,p0.v,p1.v);
				});
			} else {
				//
				shared_array2<vec2r> cell_velocity(velocity.shape());
				velocity.convert_to_full(*cell_velocity.get());
				//
				g.color4(1.0,1.0,1.0,0.5);
				velocity.shape().for_each([&](int i, int j) {
					vec2d p0 = m_dx*vec2i(i,j).cell();
					vec2d p1 = p0+m_dx*cell_velocity()(i,j);
					graphics_utility::draw_arrow(g,p0.v,p1.v);
				});
			}
		}
	}
	//
	virtual void visualize_scalar( graphics_engine &g, const macarray2<Real> &array ) const override {
		//
		Real maxv = std::numeric_limits<Real>::min();
		Real minv = std::numeric_limits<Real>::max();
		array.const_serial_actives([&](int dim, int i, int j, const auto &it) {
			Real value = it();
			maxv = std::max(maxv,value);
			minv = std::min(minv,value);
		});
		double det = maxv-minv;
		if( std::abs(det) > 1e-2 ) {
			//
			g.line_width(2.0);
			g.begin(graphics_engine::MODE::LINES);
			array.const_serial_all([&]( int dim, int i, int j, const auto &it ) {
				//
				auto set_color = [&](unsigned i, unsigned j) {
					Real v = array[dim](i,j);
					double normp = v ? 2.0*(v-minv)/det-1.0 : 0.0;
					g.color4(normp>0,0.3,normp<=0,std::abs(normp));
				};
				//
				set_color(i,j);
				g.vertex2v((m_dx*vec2d(i,j)).v);
				g.vertex2v((m_dx*vec2d(i+(dim!=0),j+(dim!=1))).v);
			});
			g.end();
			g.line_width(1.0);
		}
	}
	//
	virtual void initialize( const shape2 &shape, double dx ) override {
		m_dx = dx;
	}
	virtual void initialize( const filestream &file ) override {
		file.r(m_dx);
	}
	virtual void serialize( const filestream &file ) const override {
		file.w(m_dx);
	}
	virtual void configure( configuration &config ) override {
		config.get_bool("DrawVelocity",m_param.draw_velocity,"Should draw velocity");
		config.get_bool("DrawMAC",m_param.draw_mac,"Draw staggered velocity");
	}
	//
	struct Parameters {
		bool draw_velocity {true};
		bool draw_mac {true};
	};
	//
	Parameters m_param;
	double m_dx;
};
//
extern "C" module * create_instance() {
	return new macvisualizer2();
}
//
extern "C" const char *license() {
	return "MIT";
}
//