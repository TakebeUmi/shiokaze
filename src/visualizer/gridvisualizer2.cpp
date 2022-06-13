/*
**	gridvisualizer2.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on July 28, 2017.
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
#include <shiokaze/utility/gridutility2_interface.h>
#include <shiokaze/visualizer/gridvisualizer2_interface.h>
#include <shiokaze/utility/meshutility2_interface.h>
#include <shiokaze/cellmesher/cellmesher2_interface.h>
#include <shiokaze/array/shared_array2.h>
#include <shiokaze/graphics/graphics_utility.h>
#include <shiokaze/core/console.h>
#include <algorithm>
//
SHKZ_USING_NAMESPACE
//
class gridvisualizer2 : public gridvisualizer2_interface {
protected:
	//
	virtual void draw_active( graphics_engine &g, const array2<Real> &q ) const override {
		if( m_param.draw_active ) {
			g.color4(1.0,0.0,0.0,0.25);
			q.const_serial_actives([&](int i, int j, const auto &it) {
				g.begin(graphics_engine::MODE::TRIANGLE_FAN);
				g.vertex2(i*m_dx,j*m_dx);
				g.vertex2((i+1)*m_dx,j*m_dx);
				g.vertex2((i+1)*m_dx,(j+1)*m_dx);
				g.vertex2(i*m_dx,(j+1)*m_dx);
				g.end();
			});
		}
	}
	virtual void draw_inside( graphics_engine &g, const array2<Real> &q ) const override {
		if( m_param.draw_inside ) {
			g.color4(1.0,0.0,0.0,0.25);
			q.const_serial_inside([&](int i, int j, const auto &it) {
				g.begin(graphics_engine::MODE::TRIANGLE_FAN);
				g.vertex2(i*m_dx,j*m_dx);
				g.vertex2((i+1)*m_dx,j*m_dx);
				g.vertex2((i+1)*m_dx,(j+1)*m_dx);
				g.vertex2(i*m_dx,(j+1)*m_dx);
				g.end();
			});
		}
	}
	virtual void draw_grid( graphics_engine &g ) const override {
		if( m_param.draw_grid ) {
			double fgc[4]; g.get_foreground_color(fgc);
			g.color4(fgc[0],fgc[1],fgc[2],0.4);
			g.begin(graphics_engine::MODE::LINES);
			double width (m_dx * m_shape.w);
			double height (m_dx * m_shape.h);
			for( unsigned i=0; i<m_shape.w+1; i++ ) {
				g.vertex2(i*m_dx,0.0);
				g.vertex2(i*m_dx,height);
			}
			for( unsigned j=0; j<m_shape.h+1; j++ ) {
				g.vertex2(0.0,j*m_dx);
				g.vertex2(width,j*m_dx);
			}
			g.end();
		}
	}
	virtual void draw_density( graphics_engine &g, const array2<Real> &density ) const override {
		if( m_param.draw_density ) {
			density.const_serial_actives([&](int i, int j, const auto &it) {
				g.color4(1.0,1.0,1.0,it());
				g.begin(graphics_engine::MODE::TRIANGLE_FAN);
				g.vertex2(i*m_dx,j*m_dx);
				g.vertex2((i+1)*m_dx,j*m_dx);
				g.vertex2((i+1)*m_dx,(j+1)*m_dx);
				g.vertex2(i*m_dx,(j+1)*m_dx);
				g.end();
			});
		}
	}
	virtual void draw_velocity( graphics_engine &g, const array2<vec2r> &velocity ) const override {
		if( m_param.draw_velocity ) {
			//
			double fgc[4]; g.get_foreground_color(fgc);
			g.color4(fgc[0],fgc[1],fgc[2],0.5);
			velocity.const_serial_actives([&](int i, int j, const auto &it) {
				vec2d p0 = m_dx*vec2i(i,j).cell();
				vec2d p1 = p0+m_dx*it();
				graphics_utility::draw_arrow(g,p0.v,p1.v);
			});
		}
	}
	virtual void draw_levelset( graphics_engine &g, const array2<Real> &levelset, bool draw_contour=true ) const override {
		//
		// Cell Centered
		vec2d origin;
		if( levelset.shape() == m_shape.cell()) origin = m_dx*vec2d(0.5,0.5);
		//
		// Paint fluid
		(levelset.shape()-shape2(1,1)).for_each([&](int i, int j) {
			vec2d p[8];	int pnum; double v[2][2]; vec2d vertices[2][2];
			int count (0);
			for( int ni=0; ni<2; ni++ ) for( int nj=0; nj<2; nj++ ) {
				if(levelset.filled(i+ni,j+nj) || levelset.active(i+ni,j+nj)) count ++;
			}
			if( count == 4 ) {
				for( int ni=0; ni<2; ni++ ) for( int nj=0; nj<2; nj++ ) {
					v[ni][nj] = levelset(i+ni,j+nj);
					vertices[ni][nj] = m_dx*vec2d(i+ni,j+nj);
				}
				m_meshutility->march_points(v,vertices,p,pnum,true);
				g.begin(graphics_engine::MODE::TRIANGLE_FAN);
				for( int m=0; m<pnum; m++ ) {
					g.vertex2v((p[m]+origin).v);
				}
				g.end();
			}
		});
		//
		if( draw_contour ) {
#if 0
			//
			// Draw contour
			(levelset.shape()-shape2(1,1)).for_each([&](int i, int j) {
				//
				int count (0);
				for( int ni=0; ni<2; ni++ ) for( int nj=0; nj<2; nj++ ) {
					if(levelset.filled(i+ni,j+nj) || levelset.active(i+ni,j+nj)) count ++;
				}
				if( count == 4 ) {
					//
					// Check the nearst possible surfaces
					vec2d lines[8];	int pnum; double v[2][2]; vec2d vertices[2][2];
					for( int ni=0; ni<2; ni++ ) for( int nj=0; nj<2; nj++ ) {
						v[ni][nj] = levelset(i+ni,j+nj);
						vertices[ni][nj] = m_dx*vec2d(i+ni,j+nj);
					}
					m_meshutility->march_points(v,vertices,lines,pnum,false);
					g.color4(1.0,1.0,1.0,1.0);
					g.begin(graphics_engine::MODE::LINES);
					for( int m=0; m<pnum; m++ ) {
						g.vertex2v((lines[m]+origin).v);
					}
					g.end();
				}
			});
#else
			std::vector<vec2d> vertices;
			std::vector<std::vector<size_t> > faces;
			m_mesher->generate_contour(levelset,vertices,faces);
			//
			double fgc[4]; g.get_foreground_color(fgc);
			g.color4(fgc[0],fgc[1],fgc[2],1.0);
			//
			for( size_t i=0; i<faces.size(); i++ ) {
				g.begin(graphics_engine::MODE::LINES);
				for( unsigned j=0; j<faces[i].size(); j++ ) g.vertex2v(vertices[faces[i][j]].v);
				g.end();
			}
#endif
		}
	}
	virtual void draw_solid( graphics_engine &g, const array2<Real> &solid ) const override {
		if( m_param.draw_solid ) {
			g.color4(0.9,0.6,0.3,0.5);
			draw_levelset(g,solid);
		}
	}
	virtual void draw_fluid( graphics_engine &g, const array2<Real> &solid, const array2<Real> &fluid ) const override {
		if( m_param.draw_fluid ) {
			//
			g.color4(0.5,0.6,1.0,0.5);
			if( m_param.enclose_fluid ) {
				shared_array2<Real> combined(fluid.type());
				m_gridutility->combine_levelset(solid,fluid,*combined.get());
				draw_levelset(g,*combined.get());
			} else {
				draw_levelset(g,fluid);
			}
			//
			if( m_param.draw_active_levelset ) {
				draw_active(g,fluid);
			}
		}
	}
	virtual void visualize_cell_scalar( graphics_engine &g, const array2<Real> &q, double min_value, double max_value ) const override {
		//
		double maxv = -1e18;
		double minv = 1e18;
		double alpha = 0.5;
		//
		if( ! min_value && ! max_value ) {
			q.const_serial_actives([&](int i, int j, const auto &it) {
				double value = it();
				maxv = std::max(maxv,value);
				minv = std::min(minv,value);
			});
		} else {
			maxv = max_value;
			minv = min_value;
		}
		//
		double det = maxv-minv;
		if( std::abs(det) > 1e-2 ) {
			(q.shape()-shape2(1,1)).for_each([&](int i, int j) {
				auto set_color = [&](unsigned i, unsigned j) {
					if( q.active(i,j)) {
						double v = q(i,j);
						double normp = v ? 2.0*(v-minv)/det-1.0 : 0.0;
						g.color4(normp>0,0.3,normp<=0,alpha*std::abs(normp));
					} else {
						g.color4(0.0,0.0,0.0,0.0);
					}
				};
				g.begin(graphics_engine::MODE::TRIANGLE_FAN);
				set_color(i,j);
				g.vertex2((i+0.5)*m_dx,(j+0.5)*m_dx);
				set_color(i+1,j);
				g.vertex2((i+1.5)*m_dx,(j+0.5)*m_dx);
				set_color(i+1,j+1);
				g.vertex2((i+1.5)*m_dx,(j+1.5)*m_dx);
				set_color(i,j+1);
				g.vertex2((i+0.5)*m_dx,(j+1.5)*m_dx);
				g.end();
			});
		}
	}
	virtual void visualize_nodal_scalar( graphics_engine &g, const array2<Real> &q, double min_value, double max_value ) const override {
		//
		double maxv = -1e18;
		double minv = 1e18;
		double alpha = 0.5;
		//
		if( ! min_value && ! max_value ) {
			q.const_serial_actives([&](int i, int j, const auto &it) {
				double value = it();
				maxv = std::max(maxv,value);
				minv = std::min(minv,value);
			});
		} else {
			maxv = max_value;
			minv = min_value;
		}
		double det = maxv-minv;
		if( std::abs(det) > 1e-2 ) {
			(q.shape()-shape2(1,1)).for_each([&](int i, int j) {
				auto set_color = [&](unsigned i, unsigned j) {
					if( q.active(i,j)) {
						double normp = q(i,j) ? 2.0*(q(i,j)-minv)/det-1.0 : 0.0;
						g.color4(normp>0,0.3,normp<=0,alpha*std::abs(normp));
					} else {
						g.color4(0.0,0.0,0.0,0.0);
					}
				};
				g.begin(graphics_engine::MODE::TRIANGLE_FAN);
				set_color(i,j);
				g.vertex2(i*m_dx,j*m_dx);
				set_color(i+1,j);
				g.vertex2((i+1)*m_dx,j*m_dx);
				set_color(i+1,j+1);
				g.vertex2((i+1)*m_dx,(j+1)*m_dx);
				set_color(i,j+1);
				g.vertex2(i*m_dx,(j+1)*m_dx);
				g.end();
			});
		}
	}
	//
	virtual void configure( configuration &config ) override {
		config.get_bool("DrawActive",m_param.draw_active,"Should draw active");
		config.get_bool("DrawInside",m_param.draw_inside,"Should draw inside");
		config.get_bool("DrawGrid",m_param.draw_grid,"Should draw grid");
		config.get_bool("DrawSolid",m_param.draw_solid,"Should draw solid");
		config.get_bool("DrawDensity",m_param.draw_density,"Should draw density");
		config.get_bool("DrawVelocity",m_param.draw_velocity,"Should draw velocity");
		config.get_bool("DrawActiveFluid",m_param.draw_active_levelset,"Should draw active fluid");
		config.get_bool("EncloseFluid",m_param.enclose_fluid,"Should enclose liquid");
	}
	virtual void initialize( const shape2 &shape, double dx ) override {
		m_shape = shape;
		m_dx = dx;
	}
	virtual void initialize( const filestream &file ) override {
		file.r(m_shape);
		file.r(m_dx);
	}
	virtual void serialize( const filestream &file ) const override {
		file.w(m_shape);
		file.w(m_dx);
	}
	//
	shape2 m_shape;
	double m_dx;
	//
	struct Parameters {
		bool draw_active {true};
		bool draw_inside {true};
		bool draw_grid {true};
		bool draw_solid {true};
		bool draw_fluid {true};
		bool draw_density {true};
		bool draw_velocity {true};
		bool draw_active_levelset {true};
		bool enclose_fluid {false};
	};
	//
	Parameters m_param;
	//
	gridutility2_driver m_gridutility{this,"gridutility2"};
	meshutility2_driver m_meshutility{this,"meshutility2"};
	cellmesher2_driver m_mesher{this,"marchingsquare"};
};
//
extern "C" module * create_instance() {
	return new gridvisualizer2();
}
//
extern "C" const char *license() {
	return "MIT";
}
//