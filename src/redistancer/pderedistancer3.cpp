/*
**	pderedistancer3.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on Sep 5, 2017.
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
#include <shiokaze/redistancer/redistancer3_interface.h>
#include <shiokaze/utility/gridutility3_interface.h>
#include <shiokaze/array/shared_array3.h>
#include <shiokaze/utility/utility.h>
#include <shiokaze/core/console.h>
//
SHKZ_USING_NAMESPACE
//
class pderedistancer3 : public redistancer3_interface {
protected:
	//
	LONG_NAME("PDE Redistancer 3D")
	ARGUMENT_NAME("PDERedist")
	//
	virtual void redistance( array3<Real> &phi_array, unsigned width ) override {
		//
		double half_bandwidth = width * m_dx;
		auto smoothed_sgn = []( double value, double m_dx ) {
			return value / sqrt(value*value+m_dx*m_dx);
		};
		//
		m_gridutility->trim_narrowband(phi_array);
		phi_array.flood_fill();
		//
		for( int count=0; count<width; ++count) phi_array.dilate([&](int i, int j, int k, auto &it, int tn) {
			vec3i query[] = {vec3i(i-1,j,k),vec3i(i+1,j,k),vec3i(i,j-1,k),vec3i(i,j+1,k),vec3i(i,j,k-1),vec3i(i,j,k+1)};
			double extrapolated_value (0.0);
			for( int nq=0; nq<6; nq++ ) {
				const vec3i &qi = query[nq];
				if( ! phi_array.shape().out_of_bounds(qi) ) {
					if( phi_array.active(qi) ) {
						const double &value = phi_array(qi);
						if( value < 0.0 ) {
							extrapolated_value = std::min(extrapolated_value,extrapolated_value-m_dx);
							break;
						} else {
							extrapolated_value = std::max(extrapolated_value,extrapolated_value+m_dx);
							break;
						}
					}
				}
			}
			it.set(extrapolated_value);
		});
		//
		shared_array3<Real> phi_array0 (phi_array);
		//
		shared_array3<Real> smoothed_sgns (phi_array.type());
		smoothed_sgns->activate_as(phi_array);
		smoothed_sgns->parallel_actives([&](int i, int j, int k, auto &it, int tn) {
			it.set(smoothed_sgn(phi_array0()(i,j,k),m_dx));
		});
		//
		auto derivative = [&]( const array3<Real> &phi_array, array3<Real> &phi_array_derivative ) {
			//
			phi_array_derivative.clear();
			phi_array_derivative.activate_as(phi_array);
			//
			phi_array_derivative.parallel_actives([&](int i, int j, int k, auto &it, int tn) {
				//
				const double &sgn0 = smoothed_sgns()(i,j,k);
				const double &phi0 = phi_array0()(i,j,k);
				const double &phi = phi_array(i,j,k);
				//
				// Evaluate upwind gradient
				vec3d gradient;
				const vec3i ijk (i,j,k);
				for( int dim : DIMS3 ) {
					int select_direction (0);
					double tmp_phi (sgn0*phi), phi_backward, phi_backward0, phi_forward, phi_forward0, g (0.0);
					if( ijk[dim] > 0 && phi_array.active(i-(dim==0),j-(dim==1),k-(dim==2))) {
						phi_backward = phi_array(i-(dim==0),j-(dim==1),k-(dim==2));
						if( sgn0*phi_backward < tmp_phi ) select_direction = -1;
					}
					if( ijk[dim] < phi_array.shape()[dim]-1 && phi_array.active(i+(dim==0),j+(dim==1),k+(dim==2))) {
						phi_forward = phi_array(i+(dim==0),j+(dim==1),k+(dim==2));
						if( sgn0*phi_forward < tmp_phi ) {
							if( select_direction == 0 ) select_direction = 1;
							else if( sgn0*phi_forward < sgn0*phi_backward ) {
								select_direction = 1;
							}
						}
					}
					if( select_direction == -1 ) {
						phi_backward0 = phi_array0()(i-(dim==0),j-(dim==1),k-(dim==2));
						double frac = utility::fraction(phi0,phi_backward0);
						if( frac == 1.0 || frac == 0.0 ) {
							g = (phi-phi_backward) / m_dx;
						} else {
							if( sgn0 < 0.0 ) g = phi / (m_dx * frac);
							else g = phi / (m_dx * (1.0-frac));
						}
					} else if( select_direction == 1 ) {
						phi_forward0 = phi_array0()(i+(dim==0),j+(dim==1),k+(dim==2));
						double frac = utility::fraction(phi0,phi_forward0);
						if( frac == 1.0 || frac == 0.0 ) {
							g = (phi_forward-phi) / m_dx;
						} else {
							if( sgn0 < 0.0 ) g = -phi / (m_dx * frac);
							else g = -phi / (m_dx * (1.0-frac));
						}
					}
					gradient[dim] = g;
				}
				//
				it.set(sgn0 * (1.0 - gradient.len()));
			});
		};
		//
		// Evolve by PDE
		double dt = m_param.rate * m_dx;
		for( int itcount=0; itcount < std::ceil(((double)width)/m_param.rate); ++itcount ) {
			//
			if( m_param.temporal_scheme == "Euler" ) {
				//
				shared_array3<Real>	phi_array_derivative0(phi_array.type()),
										phi_array_save(phi_array);
				//
				phi_array.set_touch_only_actives(true);
				phi_array_derivative0->set_touch_only_actives(true);
				//
				derivative(phi_array,phi_array_derivative0());
				phi_array = phi_array_derivative0();
				phi_array *= dt;
				phi_array += phi_array_save();
				//
			} else if( m_param.temporal_scheme == "RK2" ) {
				//
				shared_array3<Real>	phi_array_derivative0(phi_array.type()),
										phi_array_derivative1(phi_array.type()),
										phi_array_derivative_tmp(phi_array.type());
				//
				phi_array.set_touch_only_actives(true);
				phi_array_derivative0->set_touch_only_actives(true);
				phi_array_derivative1->set_touch_only_actives(true);
				phi_array_derivative_tmp->set_touch_only_actives(true);
				//
				derivative(phi_array,phi_array_derivative0());
				phi_array_derivative_tmp() = phi_array_derivative0();
				phi_array_derivative_tmp() *= dt;
				phi_array_derivative_tmp() += phi_array;
				//
				derivative(phi_array_derivative_tmp(),phi_array_derivative1());
				phi_array_derivative_tmp() = phi_array_derivative0();
				phi_array_derivative_tmp() += phi_array_derivative1();
				phi_array_derivative_tmp() *= 0.5 * dt;
				//
				phi_array += phi_array_derivative_tmp();
				//
			} else {
				console::dump( "Unknown scheme %s\n", m_param.temporal_scheme.c_str());
				exit(0);
			}
		}
		//
		phi_array.parallel_actives([&](auto &it) {
			if( std::abs(it()) > m_dx*(double)width ) it.set_off();
		});
		//
		phi_array.set_as_levelset(m_dx*(double)width);
		phi_array.flood_fill();
	}
	//
protected:
	//
	virtual void configure( configuration &config ) override {
		config.get_double("IntegrationRate",m_param.rate,"Levelset advance rate");
		config.get_string("RedistTemporalScheme",m_param.temporal_scheme,"Temporal integration scheme");
	}
	//
	virtual void initialize( const shape3 &shape, double dx ) override {
		m_dx = dx;
	}
	virtual void initialize( const filestream &file ) override {
		file.r(m_dx);
	}
	virtual void serialize( const filestream &file ) const override {
		file.w(m_dx);
	}
	//
	struct Parameters {
		double rate {0.75};
		std::string temporal_scheme {"Euler"};
	};
	Parameters m_param;
	gridutility3_driver m_gridutility{this,"gridutility3"};
	double m_dx;
};
//
extern "C" module * create_instance() {
	return new pderedistancer3;
}
//
extern "C" const char *license() {
	return "MIT";
}
//