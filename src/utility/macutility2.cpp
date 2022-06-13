/*
**	macutility2.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on April 19, 2017.
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
#include <shiokaze/utility/macutility2_interface.h>
#include <shiokaze/array/array_interpolator2.h>
#include <shiokaze/array/array_extrapolator2.h>
#include <shiokaze/array/array_derivative2.h>
#include <shiokaze/array/array_utility2.h>
#include <shiokaze/array/macarray_interpolator2.h>
#include <shiokaze/array/macarray_extrapolator2.h>
#include <shiokaze/array/shared_bitarray2.h>
#include <shiokaze/cellmesher/cellmesher2_interface.h>
#include <shiokaze/math/WENO2.h>
#include <shiokaze/utility/utility.h>
#include <algorithm>
//
SHKZ_USING_NAMESPACE
using namespace array_utility2;
using namespace array_interpolator2;
//
class macutility2 : public macutility2_interface {
protected:
	//
	virtual double compute_max_u ( const macarray2<Real> &velocity ) const override {
		//
		shared_array2<vec2r> cell_velocity(m_shape);
		velocity.convert_to_full(cell_velocity());
		//
		std::vector<double> max_u_t(cell_velocity->get_thread_num(),0.0);
		cell_velocity->parallel_actives([&]( int i, int j, auto &it, int tn ) {
			max_u_t[tn] = std::max(max_u_t[tn],it().len());
		});
		double max_u (0.0);
		for( double u : max_u_t ) max_u = std::max(max_u,u);
		return max_u;
	}
	virtual void constrain_velocity( const array2<Real> &solid, macarray2<Real> &velocity ) const override {
		//
		shared_macarray2<Real> velocity_save(velocity);
		if( levelset_exist(solid) ) {
			//
			velocity.parallel_actives([&](int dim, int i, int j, auto &it, int tn) {
				vec2i pi(i,j);
				vec2d p(vec2i(i,j).face(dim));
				if( interpolate<Real>(solid,p) < 0.0 ) {
					Real derivative[DIM2];
					array_derivative2::derivative(solid,p,derivative);
					vec2d normal = vec2d(derivative)/m_dx;
					if( normal.norm2() ) {
						vec2d u = macarray_interpolator2::interpolate<Real>(velocity_save(),p);
						if( u * normal < 0.0 ) {
							it.set((u-normal*(u*normal))[dim]);
						}
					}
				}
			});
		}
		//
		velocity.parallel_actives([&](int dim, int i, int j, auto &it, int tn) {
				vec2i pi(i,j);
			if( dim == 1 ) {
				if( pi[dim]==0 && it() < 0.0 ) it.set(0.0);
				if( pi[dim]==m_shape[dim] && it() > 0.0 ) it.set(0.0);
			}
		});
	}
	virtual void extrapolate_and_constrain_velocity( const array2<Real> &solid, macarray2<Real> &velocity, int extrapolate_width ) const override {
		//
		macarray_extrapolator2::extrapolate(velocity,extrapolate_width);
		constrain_velocity(solid,velocity);
	}
	virtual void compute_area_fraction( const array2<Real> &solid, macarray2<Real> &areas, bool enclose_domain_boundary=true ) const override {
		//
		if( levelset_exist(solid) ) {
			//
			areas.clear(0.0);
			m_parallel.for_each( DIM2, [&]( size_t dim ) {
				areas[dim].activate_as(solid);
				areas[dim].activate_as(solid,-vec2i(dim!=0,dim!=1));
				areas[dim].set_as_fillable(1.0);
			});
			//
			areas.parallel_actives([&](int dim, int i, int j, auto &it, int tn) {
				double area;
				vec2i pi(i,j);
				if( enclose_domain_boundary && (pi[dim] == 0 || pi[dim] == m_shape[dim] )) area = 0.0;
				else area = 1.0-utility::fraction(solid(i,j),solid(i+(dim!=0),j+(dim!=1)));
				if( area && area < m_param.eps_solid ) area = m_param.eps_solid;
				it.set(area);
			});
			//
			m_parallel.for_each( DIM2, [&]( size_t dim ) {
				areas[dim].flood_fill();
			});
			//
		} else {
			//
			areas.clear(1.0);
			if( enclose_domain_boundary ) {
				for( int i=0; i<m_shape.w; ++i ) {
					areas[1].set(i,0,0.0);
					areas[1].set(i,m_shape.h,0.0);
				}
				for( int j=0; j<m_shape.h; ++j ) {
					areas[0].set(0,j,0.0);
					areas[0].set(m_shape.w,j,0.0);
				}
			}
		}
	}
	virtual void compute_fluid_fraction( const array2<Real> &fluid, macarray2<Real> &rhos ) const override {
		//
		if( levelset_exist(fluid)) {
			//
			rhos.clear(0.0);
			m_parallel.for_each( DIM2, [&]( size_t dim ) {
				rhos[dim].activate_as(fluid);
				rhos[dim].activate_as(fluid,vec2i(dim==0,dim==1));
				rhos[dim].set_as_fillable(1.0);
			});
			//
			rhos.parallel_actives([&](int dim, int i, int j, auto &it, int tn) {
				//
				double rho = utility::fraction(
					fluid(m_shape.clamp(i,j)),
					fluid(m_shape.clamp(i-(dim==0),j-(dim==1)))
				);
				if( rho && rho < m_param.eps_fluid ) rho = m_param.eps_fluid;
				it.set(rho);
			});
			//
			m_parallel.for_each( DIM2, [&]( size_t dim ) {
				rhos[dim].flood_fill();
			});
			//
		} else {
			rhos.clear(1.0);
		}
	}
	virtual void compute_face_density( const array2<Real> &solid, const array2<Real> &fluid, macarray2<Real> &density ) const override {
		//
		compute_fluid_fraction(fluid,density);
		if( levelset_exist(solid) ) {
			//
			shared_macarray2<Real> tmp_areas(density.type());
			compute_area_fraction(solid,tmp_areas());
			density.parallel_actives([&](int dim, int i, int j, auto &it, int tn) {
				it.multiply(tmp_areas()[dim](i,j));
			});
		}
	}
	double get_kinetic_energy( const macarray2<Real> &areas, const macarray2<Real> &rhos, const macarray2<Real> &velocity ) const {
		//
		std::vector<double> results (velocity.get_thread_num(),0.0);
		velocity.const_parallel_actives([&]( int dim, int i, int j, const auto &it, int tn ) {
			double area = areas[dim](i,j);
			if( area ) {
				double rho = rhos[dim](i,j);
				if( rho ) {
					double u = velocity[dim](i,j);
					double dA = (m_dx*m_dx) * (area*rho);
					results[tn] += 0.5*(u*u)*dA;
				}
			}
		});
		//
		double result (0.0);
		for( const auto &e : results ) result += e;
		return result;
	}
	virtual double get_kinetic_energy( const array2<Real> &solid, const array2<Real> &fluid, const macarray2<Real> &velocity ) const override {
		//
		shared_macarray2<Real> tmp_areas(velocity.type());
		shared_macarray2<Real> tmp_rhos(velocity.type());
		//
		compute_area_fraction(solid,tmp_areas());
		compute_fluid_fraction(fluid,tmp_rhos());
		//
		return get_kinetic_energy(tmp_areas(),tmp_rhos(),velocity);
	}
	double get_gravitational_potential_energy( const macarray2<Real> &areas, const macarray2<Real> &rhos, vec2d gravity ) const {
		//
		double sum (0.0);
		rhos.const_serial_actives([&]( int dim, int i, int j, const auto &it ) {
			const Real area = areas[dim](i,j);
			if( area ) {
				vec2d p = m_dx*vec2d(i,j).face(dim);
				sum += area * p[1] * it();
			}
		});
		for( int dim : DIMS2 ) {
			rhos[dim].const_serial_inside([&]( int i, int j, const auto &it ) {
				if( ! it.active()) {
					const Real area = areas[dim](i,j);
					if( area ) {
						vec2d p = m_dx*vec2d(i,j).face(dim);
						sum += area * p[1] * it();
					}
				}
			});
		}
		return -sum * gravity[1] * (m_dx*m_dx) / 2.0;
	}
	virtual double get_gravitational_potential_energy( const array2<Real> &solid, const array2<Real> &fluid, vec2d gravity ) const override {
		//
		shared_macarray2<Real> tmp_areas(fluid.shape());
		shared_macarray2<Real> tmp_rhos(fluid.shape());
		//
		compute_area_fraction(solid,tmp_areas());
		compute_fluid_fraction(fluid,tmp_rhos());
		//
		return get_gravitational_potential_energy(tmp_areas(),tmp_rhos(),gravity);
	}
	virtual double get_surfacetension_potential_energy( const array2<Real> &solid, const array2<Real> &fluid, double tension_coeff ) const override {
		//
		if( tension_coeff ) {
			std::vector<vec2d> vertices;
			std::vector<std::vector<size_t> > faces;
			m_mesher->generate_contour(fluid,vertices,faces);
			//
			return tension_coeff * utility::compute_length(vertices,faces,[&]( const vec2d &p ) {
				return interpolate<Real>(solid,m_dx*p) > 0.0;
			});
		} else {
			return 0.0;
		}
	}
	virtual double get_total_energy( const array2<Real> &solid, const array2<Real> &fluid, const macarray2<Real> &velocity, vec2d gravity, double tension_coeff ) const override {
		//
		const auto energy_list = get_all_kinds_of_energy(solid,fluid,velocity,gravity,tension_coeff);
		return std::get<0>(energy_list)+std::get<1>(energy_list)+std::get<2>(energy_list);
	}
	virtual std::tuple<double,double,double> get_all_kinds_of_energy( const array2<Real> &solid, const array2<Real> &fluid, const macarray2<Real> &velocity, vec2d gravity, double tension_coeff ) const override {
		//
		shared_macarray2<Real> tmp_areas(velocity.type());
		shared_macarray2<Real> tmp_rhos(velocity.type());
		//
		compute_area_fraction(solid,tmp_areas());
		compute_fluid_fraction(fluid,tmp_rhos());
		//
		return {
			get_gravitational_potential_energy(tmp_areas(),tmp_rhos(),gravity),
			get_kinetic_energy(tmp_areas(),tmp_rhos(),velocity),
			get_surfacetension_potential_energy(solid,fluid,tension_coeff)
		};
	}
	virtual void get_velocity_jacobian( const vec2d &p, const macarray2<Real> &velocity, vec2r jacobian[DIM2] ) const override {
		for( unsigned dim : DIMS2 ) {
			array_derivative2::derivative(velocity[dim],vec2d(p[0]/m_dx-0.5*(dim!=0),p[1]/m_dx-0.5*(dim!=1)),jacobian[dim].v);
			jacobian[dim] /= m_dx;
		}
	}
	virtual void assign_initial_variables( const dylibloader &dylib, 
					array2<Real> *solid, macarray2<Real> *solid_velocity,
					array2<Real> *fluid, macarray2<Real> *fluid_velocity,
					array2<Real> *density ) const override {
		//
		// Assign initial velocity
		if( fluid_velocity ) {
			auto velocity_func = reinterpret_cast<vec2d(*)(const vec2d &)>(dylib.load_symbol("velocity"));
			auto fluid_func = reinterpret_cast<double(*)(const vec2d &)>(dylib.load_symbol("fluid"));
			fluid_velocity->parallel_all([&](int dim, int i, int j, auto &it) {
				bool skip (false);
				if( fluid_func ) {
					skip = (*fluid_func)(m_dx*vec2i(i,j).face(dim)) > m_param.narrowband_dist*m_dx;
				}
				if( ! skip ) {
					vec2d value = velocity_func ? (*velocity_func)(m_dx*vec2i(i,j).face(dim)) : 0.0;
					it.set(value[dim]);
				}
			});
		}
		//
		// Assign fluid levelset
		if( fluid ) {
			//
			auto fluid_func = reinterpret_cast<double(*)(const vec2d &)>(dylib.load_symbol("fluid"));
			auto solid_func = reinterpret_cast<double(*)(const vec2d &)>(dylib.load_symbol("solid"));
			auto moving_solid_func = reinterpret_cast<std::pair<double,vec2d>(*)(double time, const vec2d &p)>(dylib.load_symbol("moving_solid"));
			if( fluid_func ) {
				if( solid_func || moving_solid_func ) {
					fluid->parallel_all([&](int i, int j, auto &it) {
						vec2d p = m_dx*vec2i(i,j).cell();
						double value = (*fluid_func)(p);
						double solid_value0 = solid_func ? (*solid_func)(p)+m_dx : 1.0;
						double solid_value1 = moving_solid_func ? (*moving_solid_func)(0.0,p).first+m_dx : 1.0;
						value = std::max(value,-solid_value0);
						value = std::max(value,-solid_value1);
						if( ! m_param.narrowband_dist || std::abs(value) < m_param.narrowband_dist*m_dx ) it.set(value);
					});
				} else {
					fluid->parallel_all([&](int i, int j, auto &it) {
						double value = (*fluid_func)(m_dx*vec2i(i,j).cell());
						if( ! m_param.narrowband_dist || std::abs(value) < m_param.narrowband_dist*m_dx ) it.set(value);
					});
				}
			} else {
				fluid->clear(-1.0);
			}
			//
			if( m_param.narrowband_dist ) {
				fluid->set_as_levelset(m_param.narrowband_dist*m_dx);
				fluid->flood_fill();
			}
		}
		//
		// Assign density
		if( density ) {
			auto density_func = reinterpret_cast<double(*)(const vec2d &)>(dylib.load_symbol("density"));
			if( density_func ) {
				density->parallel_all([&](int i, int j, auto &it) {
					it.set((*density_func)(m_dx*vec2i(i,j).cell()));
				});
			}
		}
		//
		// Assign solid variables
		update_solid_variables(dylib,0.0,solid,solid_velocity);
	}
	virtual void update_solid_variables( const dylibloader &dylib, double time, array2<Real> *solid, macarray2<Real> *solid_velocity ) const override {
		//
		auto moving_solid_func = reinterpret_cast<std::pair<double,vec2d>(*)(double time, const vec2d &p)>(dylib.load_symbol("moving_solid"));
		auto solid_func = reinterpret_cast<double(*)(const vec2d &)>(dylib.load_symbol("solid"));
		auto set_boundary_flux = reinterpret_cast<void(*)( double, Real [DIM2][2] )>(dylib.load_symbol("set_boundary_flux"));
		//
		if( solid ) {
			solid->clear(1.0);
			if( solid_func || moving_solid_func ) {
				//
				solid->parallel_all([&](int i, int j, auto &it) {
					const vec2d &p = m_dx*vec2i(i,j).nodal();
					double value (1.0);
					if( solid_func ) value = std::min(value,solid_func(p));
					if( moving_solid_func ) value = std::min(value,moving_solid_func(time,p).first);
					if( ! m_param.narrowband_dist || std::abs(value) < m_param.narrowband_dist*m_dx ) it.set(value);
				});
				if( m_param.narrowband_dist ) {
					solid->set_as_levelset(m_param.narrowband_dist*m_dx);
					solid->flood_fill();
				}
				//
				if( solid_velocity && moving_solid_func ) {
					//
					shared_array2<vec2d> tmp_nodal_velocity(m_shape.nodal());
					tmp_nodal_velocity->parallel_all([&](int i, int j, auto &it) {
						if( (*solid)(i,j) < 0.0 ) {
							it.set(moving_solid_func(time,m_dx*vec2i(i,j).nodal()).second);
						}
					});
					//
					const auto get_face_vertices = [&]( int dim, const vec2i &pi ) {
						std::vector<vec2i> result(2);
						result[0] = pi;
						result[1] = pi+vec2i(dim!=0,dim!=1);
						return result;
					};
					//
					solid_velocity->clear();
					solid_velocity->parallel_all([&]( int dim, int i, int j, auto &it ) {
						std::vector<vec2i> results = get_face_vertices(dim,vec2i(i,j));
						double wsum (0.0), sum (0.0);
						for( const auto &pi : results ) {
							if( (*solid)(pi) < 0.0 ) {
								wsum += 1.0;
								sum += tmp_nodal_velocity()(pi)[dim];
							}
						}
						if( wsum ) {
							it.set(sum/wsum);
						}
					});
				}
			}
		}
		//
		if( set_boundary_flux && solid_velocity ) {
			Real flux[DIM2][2] = {{0.0,0.0},{0.0,0.0}};
			set_boundary_flux(time,flux);
			for( int i=0; i<m_shape.w; ++i ) {
				(*solid_velocity)[1].set(i,0,flux[1][0]);
				(*solid_velocity)[1].set(i,m_shape.h,flux[1][1]);
			}
			for( int j=0; j<m_shape.h; ++j ) {
				(*solid_velocity)[0].set(0,j,flux[0][0]);
				(*solid_velocity)[0].set(m_shape.w,j,flux[0][1]);
			}
		}
		//
		auto velocity_func = reinterpret_cast<vec2d(*)(double, const vec2d &)>(dylib.load_symbol("boundary_velocity"));
		if( velocity_func && solid_velocity ) {
			solid_velocity->parallel_all([&]( int dim, int i, int j, auto &it ) {
				const double value = velocity_func(time,m_dx*vec2i(i,j).face(dim))[dim];
				if( value ) it.set(value);
				else it.set_off();
			});
		}
	}
	virtual void add_force( vec2d p, vec2d f, macarray2<Real> &external_force ) const override {
		for( unsigned dim : DIMS2 ) {
			vec2d index_coord = p/m_dx-vec2d(0.5,0.5);
			external_force[dim].set(m_shape.face(dim).clamp(index_coord),f[dim]);
		}
	}
	//
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
	virtual void configure( configuration &config ) override {
		config.get_double("EpsFluid",m_param.eps_fluid,"Minimal bound for fluid fraction");
		config.get_double("EpsSolid",m_param.eps_solid,"Minimal bound for solid fraction");
		config.get_bool("WENO",m_param.weno_interpolation,"Whether to use WENO interpolation");
		config.get_double("NarrowBandDist",m_param.narrowband_dist);
	}
	//
	struct Parameters {
		//
		double eps_fluid {1e-2};
		double eps_solid {1e-2};
		bool weno_interpolation {false};
		double narrowband_dist {sqrt(2.0)};
	};
	//
	Parameters m_param;
	double m_dx;
	shape2 m_shape;
	parallel_driver m_parallel{this};
	cellmesher2_driver m_mesher{this,"marchingsquare"};
};
//
extern "C" module * create_instance() {
	return new macutility2();
}
//
extern "C" const char *license() {
	return "MIT";
}
//