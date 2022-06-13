/*
**	flatrasterizer3.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on July 27, 2017.
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
// Yongning Zhu and Robert Bridson. 2005. Animating sand as a fluid. In ACM SIGGRAPH 2005 Papers (SIGGRAPH '05), Markus Gross (Ed.). ACM, New York, NY, USA, 965-972. DOI: https://doi.org/10.1145/1186822.1073298
//
#include <shiokaze/pointgridhash/pointgridhash3_interface.h>
#include <shiokaze/particlerasterizer/particlerasterizer3_interface.h>
#include <cmath>
//
SHKZ_USING_NAMESPACE
//
class flatrasterizer3 : public particlerasterizer3_interface {
protected:
	//
	LONG_NAME("Flat Rasterizer 3D")
	//
	virtual void build_levelset( array3<Real> &fluid, const bitarray3 &mask, const std::vector<Particle3> &particles ) const override {
		//
		std::vector<vec3r> points(particles.size());
		for( size_t n=0; n<points.size(); ++n ) points[n] = particles[n].p;
		const_cast<pointgridhash3_driver &>(m_pointgridhash)->sort_points(points);
		//
		auto cube = []( double x ) { return x*x*x; };
		auto kernel = [cube]( double s2 ) {
			return std::max(0.0,cube(1.0-s2));
		};
		//
		fluid.clear();
		fluid.activate_as_bit(mask);
		fluid.parallel_actives( [&](int i, int j, int k, auto &it ) {
			if( mask(i,j,k)) {
				double result (1.0);
				vec3d x (m_dx*vec3i(i,j,k).cell());
				std::vector<size_t> neighbors = m_pointgridhash->get_cell_neighbors(vec3i(i,j,k),pointgridhash3_interface::USE_NODAL);
				if( neighbors.size()) {
					double R (0.0), w_sum (0.0);
					for( const auto n : neighbors ) R += m_param.r_factor * particles[n].r;
					R /= neighbors.size(); R *= 2.0 * m_param.w_factor;
					if( R*R ) {
						for( const auto n : neighbors ) w_sum += kernel((particles[n].p-x).norm2() / (R*R));
						if( w_sum ) {
							vec3d avg_x; double avg_r (0.0);
							for( const auto n : neighbors ) {
								double w = kernel((particles[n].p-x).norm2() / (R*R)) / w_sum;
								if( w ) {
									avg_x += w * particles[n].p;
									avg_r += w * m_param.r_factor * particles[n].r;
								}
							}
							result = (avg_x-x).len() - avg_r;
						}
					}
				}
				it.set(result);
			}
		});
	}
	//
	virtual void configure( configuration &config ) override {
		config.get_double("RadiusFactor",m_param.r_factor,"Radius exaggeration factor");
		config.get_double("WeightFactor",m_param.w_factor,"Weight factor");
	}
	virtual void initialize( const shape3 &shape, double dx ) override {
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
	typedef struct {
		double r_factor {1.5};
		double w_factor {1.0};
	} Parameters;
	Parameters m_param;
	//
	pointgridhash3_driver m_pointgridhash{this,"pointgridhash3"};
	shape3 m_shape;
	double m_dx;
};
//
extern "C" module * create_instance() {
	return new flatrasterizer3();
}
//
extern "C" const char *license() {
	return "MIT";
}
//