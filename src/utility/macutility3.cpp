/*
**	macutility3.cpp
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
#include <shiokaze/utility/macutility3_interface.h>
#include <shiokaze/array/array_interpolator3.h>
#include <shiokaze/array/array_extrapolator3.h>
#include <shiokaze/array/macarray_interpolator3.h>
#include <shiokaze/array/macarray_extrapolator3.h>
#include <shiokaze/array/shared_bitarray3.h>
#include <shiokaze/array/array_derivative3.h>
#include <shiokaze/array/array_utility3.h>
#include <shiokaze/cellmesher/cellmesher3_interface.h>
#include <shiokaze/core/dylibloader.h>
#include <shiokaze/core/console.h>
#include <shiokaze/core/scoped_timer.h>
#include <shiokaze/utility/utility.h>
#include <shiokaze/math/WENO3.h>
#include <algorithm>
//
SHKZ_USING_NAMESPACE
using namespace array_utility3;
using namespace array_interpolator3;
//
class macutility3 : public macutility3_interface {
protected:
	//
	virtual double compute_max_u ( const macarray3<Real> &velocity ) const override {
		//
		shared_array3<vec3r> cell_velocity(m_shape);
		velocity.convert_to_full(cell_velocity());
		//
		std::vector<double> max_u_t(cell_velocity->get_thread_num(),0.0);
		cell_velocity->parallel_actives([&]( int i, int j, int k, auto &it, int tn ) {
			max_u_t[tn] = std::max(max_u_t[tn],it().len());
		});
		double max_u (0.0);
		for( double u : max_u_t ) max_u = std::max(max_u,u);
		return max_u;
	}
	virtual void constrain_velocity( const array3<Real> &solid, macarray3<Real> &velocity ) const override {
		//
		shared_macarray3<Real> velocity_save = shared_macarray3<Real>(velocity.type());
		velocity_save->copy(velocity);
		//
		if( levelset_exist(solid) ) {
			//
			velocity.parallel_actives([&](int dim, int i, int j, int k, auto &it, int tn) {
				vec3i pi(i,j,k);
				vec3d p(vec3i(i,j,k).face(dim));
				if( interpolate<Real>(solid,p) < 0.0 ) {
					Real derivative[DIM3];
					array_derivative3::derivative(solid,p,derivative);
					vec3d normal = vec3d(derivative)/m_dx;
					if( normal.norm2() ) {
						vec3d u = macarray_interpolator3::interpolate<Real>(velocity_save(),p);
						if( u * normal < 0.0 ) {
							it.set((u-normal*(u*normal))[dim]);
						}
					}
				}
			});
		}
		//
		velocity.parallel_actives([&](int dim, int i, int j, int k, auto &it, int tn) {
				vec3i pi(i,j,k);
			if( dim == 1 ) {
				if( pi[dim]==0 && it() < 0.0 ) it.set(0.0);
				if( pi[dim]==m_shape[dim] && it() > 0.0 ) it.set(0.0);
			}
		});
	}
	virtual void extrapolate_and_constrain_velocity( const array3<Real> &solid, macarray3<Real> &velocity, int extrapolate_width ) const override {
		//
		macarray_extrapolator3::extrapolate(velocity,extrapolate_width);
		constrain_velocity(solid,velocity);
	}
	virtual void compute_area_fraction( const array3<Real> &solid, macarray3<Real> &areas, bool enclose_domain_boundary=true ) const override {
		//
		if( levelset_exist(solid) ) {
			//
			areas.clear(0.0);
			m_parallel.for_each( DIM3, [&]( size_t dim ) {
				areas[dim].activate_as(solid);
				if( dim == 0 ) {
					for( int jj=-1; jj<=0; ++jj ) for( int kk=-1; kk<=0; ++kk ) {
						areas[dim].activate_as(solid,vec3i(0,jj,kk));
					}
				} else if( dim == 1 ) {
					for( int ii=-1; ii<=0; ++ii ) for( int kk=-1; kk<=0; ++kk ) {
						areas[dim].activate_as(solid,vec3i(ii,0,kk));
					}
				} else if( dim == 2 ) {
					for( int ii=-1; ii<=0; ++ii ) for( int jj=-1; jj<=0; ++jj ) {
						areas[dim].activate_as(solid,vec3i(ii,jj,0));
					}
				}
				areas[dim].set_as_fillable(1.0);
			});
			//
			areas.parallel_actives([&](int dim, int i, int j, int k, auto &it, int tn) {
				double area;
				vec3i pi(i,j,k);
				if( enclose_domain_boundary && (pi[dim] == 0 || pi[dim] == m_shape[dim] )) area = 0.0;
				else {
					double quadsolid[2][2];
					if( dim == 0 ) {
						quadsolid[0][0] = solid(i,j,k);
						quadsolid[1][0] = solid(i,j+1,k);
						quadsolid[1][1] = solid(i,j+1,k+1);
						quadsolid[0][1] = solid(i,j,k+1);
					} else if( dim == 1 ) {
						quadsolid[0][0] = solid(i,j,k);
						quadsolid[1][0] = solid(i+1,j,k);
						quadsolid[1][1] = solid(i+1,j,k+1);
						quadsolid[0][1] = solid(i,j,k+1);
					} else if( dim == 2 ) {
						quadsolid[0][0] = solid(i,j,k);
						quadsolid[1][0] = solid(i+1,j,k);
						quadsolid[1][1] = solid(i+1,j+1,k);
						quadsolid[0][1] = solid(i,j+1,k);
					}
					area = 1.0-utility::get_area(quadsolid);
				}
				if( area && area < m_param.eps_solid ) area = m_param.eps_solid;
				it.set(area);
			});
			//
			m_parallel.for_each( DIM3, [&]( size_t dim ) {
				areas[dim].flood_fill();
			});
			//
		} else {
			//
			areas.clear(1.0);
			if( enclose_domain_boundary ) {
				for( int j=0; j<m_shape.h; ++j ) for( int k=0; k<m_shape.d; ++k ) {
					areas[0].set(0,j,k,0.0);
					areas[0].set(m_shape.w,j,k,0.0);
				}
				for( int i=0; i<m_shape.w; ++i ) for( int k=0; k<m_shape.d; ++k ) {
					areas[1].set(i,0,k,0.0);
					areas[1].set(i,m_shape.h,k,0.0);
				}
				for( int i=0; i<m_shape.w; ++i ) for( int j=0; j<m_shape.h; ++j ) {
					areas[2].set(i,j,0,0.0);
					areas[2].set(i,j,m_shape.d,0.0);
				}
			}
		}
	}
	virtual void compute_fluid_fraction( const array3<Real> &fluid, macarray3<Real> &rhos ) const override {
		//
		if( levelset_exist(fluid)) {
			//
			rhos.clear(0.0);
			m_parallel.for_each( DIM3, [&]( size_t dim ) {
				rhos[dim].activate_as(fluid);
				rhos[dim].activate_as(fluid,vec3i(dim==0,dim==1,dim==2));
				rhos[dim].set_as_fillable(1.0);
			});
			//
			rhos.parallel_actives([&](int dim, int i, int j, int k, auto &it, int tn) {
				//
				double rho = utility::fraction(
					fluid(m_shape.clamp(i,j,k)),
					fluid(m_shape.clamp(i-(dim==0),j-(dim==1),k-(dim==2)))
				);
				if( rho && rho < m_param.eps_fluid ) rho = m_param.eps_fluid;
				it.set(rho);
			});
			//
			m_parallel.for_each( DIM3, [&]( size_t dim ) {
				rhos[dim].flood_fill();
			});
			//
		} else {
			rhos.clear(1.0);
		}
	}
	virtual void compute_face_density( const array3<Real> &solid, const array3<Real> &fluid, macarray3<Real> &density ) const override {
		//
		compute_fluid_fraction(fluid,density);
		if( levelset_exist(solid) ) {
			//
			shared_macarray3<Real> tmp_areas(density.type());
			compute_area_fraction(solid,tmp_areas());
			density.parallel_actives([&](int dim, int i, int j, int k, auto &it, int tn) {
				it.multiply(tmp_areas()[dim](i,j,k));
			});
		}
	}
	double get_kinetic_energy( const macarray3<Real> &areas, const macarray3<Real> &rhos, const macarray3<Real> &velocity ) const {
		//
		std::vector<Real> results(velocity.get_thread_num(),0.0);
		velocity.const_parallel_actives([&]( int dim, int i, int j, int k, const auto &it, int tn ) {
			double area = areas[dim](i,j,k);
			if( area ) {
				double rho = rhos[dim](i,j,k);
				if( rho ) {
					double u = velocity[dim](i,j,k);
					double dV = (m_dx*m_dx*m_dx) * (area*rho);
					results[tn] += 0.5*(u*u)*dV;
				}
			}
		});
		//
		double result (0.0);
		for( const auto &e : results ) result += e;
		return result;
	}
	virtual double get_kinetic_energy( const array3<Real> &solid, const array3<Real> &fluid, const macarray3<Real> &velocity ) const override {
		//
		shared_macarray3<Real> tmp_areas(velocity.type());
		shared_macarray3<Real> tmp_rhos(velocity.type());
		//
		compute_area_fraction(solid,tmp_areas());
		compute_fluid_fraction(fluid,tmp_rhos());
		//
		return get_kinetic_energy(tmp_areas(),tmp_rhos(),velocity);
	}
	double get_gravitational_potential_energy( const macarray3<Real> &areas, const macarray3<Real> &rhos, vec3d gravity ) const {
		//
		double sum (0.0);
		rhos.const_serial_actives([&]( int dim, int i, int j, int k, const auto &it ) {
			const Real area = areas[dim](i,j,k);
			if( area ) {
				vec3d p = m_dx*vec3d(i,j,k).face(dim);
				sum += area * p[1] * it();
			}
		});
		for( int dim : DIMS3 ) {
			rhos[dim].const_serial_inside([&]( int i, int j, int k, const auto &it ) {
				if( ! it.active()) {
					const Real area = areas[dim](i,j,k);
					if( area ) {
						vec3d p = m_dx*vec3d(i,j,k).face(dim);
						sum += area * p[1] * it();
					}
				}
			});
		}
		return -sum * gravity[1] * (m_dx*m_dx*m_dx) / 3.0;
	}
	virtual double get_gravitational_potential_energy( const array3<Real> &solid, const array3<Real> &fluid, vec3d gravity ) const override {
		//
		shared_macarray3<Real> tmp_areas(fluid.shape());
		shared_macarray3<Real> tmp_rhos(fluid.shape());
		//
		compute_area_fraction(solid,tmp_areas());
		compute_fluid_fraction(fluid,tmp_rhos());
		//
		return get_gravitational_potential_energy(tmp_areas(),tmp_rhos(),gravity);
	}
	virtual double get_surfacetension_potential_energy( const array3<Real> &solid, const array3<Real> &fluid, double tension_coeff ) const override {
		//
		if( tension_coeff ) {
			std::vector<vec3d> vertices;
			std::vector<std::vector<size_t> > faces;
			m_mesher->generate_mesh(fluid,vertices,faces);
			//
			return tension_coeff * utility::compute_area(vertices,faces,[&]( const vec3d &p ) {
				return interpolate<Real>(solid,m_dx*p) > 0.0;
			});
		} else {
			return 0.0;
		}
	}
	virtual double get_total_energy( const array3<Real> &solid, const array3<Real> &fluid, const macarray3<Real> &velocity, vec3d gravity, double tension_coeff ) const override {
		//
		const auto energy_list = get_all_kinds_of_energy(solid,fluid,velocity,gravity,tension_coeff);
		return std::get<0>(energy_list)+std::get<1>(energy_list)+std::get<2>(energy_list);
	}
	virtual std::tuple<double,double,double> get_all_kinds_of_energy( const array3<Real> &solid, const array3<Real> &fluid, const macarray3<Real> &velocity, vec3d gravity, double tension_coeff ) const override {
		//
		shared_macarray3<Real> tmp_areas(velocity.type());
		shared_macarray3<Real> tmp_rhos(velocity.type());
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
	virtual void get_velocity_jacobian( const vec3d &p, const macarray3<Real> &velocity, vec3r jacobian[DIM3] ) const override {
		for( unsigned dim : DIMS3 ) {
			array_derivative3::derivative(velocity[dim],vec3d(p[0]/m_dx-0.5*(dim!=0),p[1]/m_dx-0.5*(dim!=1),p[2]/m_dx-0.5*(dim!=2)),jacobian[dim].v);
			jacobian[dim] /= m_dx;
		}
	}
	virtual void assign_initial_variables( const dylibloader &dylib, 
					array3<Real> *solid, macarray3<Real> *solid_velocity,
					array3<Real> *fluid, macarray3<Real> *fluid_velocity,
					array3<Real> *density ) const override {
		//
		// Scoped timer
		scoped_timer timer(this,"assign_initial_variables");
		//
		timer.tick(); console::dump( ">>> Assigining variables...\n" );
		//
		// Assign velocity
		const double sqrt3 = sqrt(3.0);
		if( fluid_velocity ) {
			fluid_velocity->set_touch_only_actives(true);
			auto velocity_func = reinterpret_cast<vec3d(*)(const vec3d &)>(dylib.load_symbol("velocity"));
			timer.tick(); console::dump( "Assigining velocity..." );
			auto fluid_func = reinterpret_cast<double(*)(const vec3d &)>(dylib.load_symbol("fluid"));
			fluid_velocity->parallel_all([&](int dim, int i, int j, int k, auto &it) {
				bool skip (false);
				if( fluid_func ) {
					skip = (*fluid_func)(m_dx*vec3i(i,j,k).face(dim)) > sqrt3*m_dx;
				}
				if( ! skip ) {
					vec3d value = velocity_func ? (*velocity_func)(m_dx*vec3i(i,j,k).face(dim)) : 0.0;
					it.set(value[dim]);
				}
			});
			console::dump( "Done. Took %s.\n", timer.stock("assign_velocity").c_str());
		}
		//
		// Assign fluid levelset
		if( fluid ) {
			//
			auto fluid_func = reinterpret_cast<double(*)(const vec3d &)>(dylib.load_symbol("fluid"));
			auto solid_func = reinterpret_cast<double(*)(const vec3d &)>(dylib.load_symbol("solid"));
			auto moving_solid_func = reinterpret_cast<std::pair<double,vec3d>(*)(double time, const vec3d &p)>(dylib.load_symbol("moving_solid"));
			if( fluid_func ) {
				timer.tick(); console::dump( "Assigining fluid levelset..." );
				if( solid_func || moving_solid_func ) {
					fluid->parallel_all([&](int i, int j, int k, auto &it) {
						vec3d p = m_dx*vec3i(i,j,k).cell();
						double value = (*fluid_func)(p);
						double solid_value0 = solid_func ? (*solid_func)(p)+m_dx : 1.0;
						double solid_value1 = moving_solid_func ? (*moving_solid_func)(0.0,p).first+m_dx : 1.0;
						value = std::max(value,-solid_value0);
						value = std::max(value,-solid_value1);
						if( ! m_param.narrowband_dist || std::abs(value) < m_param.narrowband_dist*m_dx ) it.set(value);
					});
				} else {
					fluid->parallel_all([&](int i, int j, int k, auto &it) {
						double value = (*fluid_func)(m_dx*vec3i(i,j,k).cell());
						if( ! m_param.narrowband_dist || std::abs(value) < m_param.narrowband_dist*m_dx ) it.set(value);
					});
				}
				console::dump( "Done. Took %s.\n", timer.stock("assign_fluid").c_str());
			} else {
				fluid->clear(-1.0);
			}
			if( m_param.narrowband_dist ) {
				fluid->set_as_levelset(m_param.narrowband_dist*m_dx);
				fluid->flood_fill();
			}
		}
		//
		// Assign density
		if( density ) {
			auto density_func = reinterpret_cast<double(*)(const vec3d &)>(dylib.load_symbol("density"));
			if( density_func ) {
				timer.tick(); console::dump( "Assigning initial density..." );
				density->parallel_all([&](int i, int j, int k, auto &it) {
					it.set((*density_func)(m_dx*vec3i(i,j,k).cell()));
				});
				console::dump( "Done. Took %s.\n", timer.stock("evaluate_density").c_str());
			}
		}
		//
		console::dump( "<<< Done. Took %s.\n", timer.stock("assign_variables").c_str());
		//
		// Assign solid variables
		update_solid_variables(dylib,0.0,solid,solid_velocity);
	}
	virtual void update_solid_variables( const dylibloader &dylib, double time, array3<Real> *solid, macarray3<Real> *solid_velocity ) const override {
		//
		auto moving_solid_func = reinterpret_cast<std::pair<double,vec3d>(*)(double time, const vec3d &p)>(dylib.load_symbol("moving_solid"));
		auto solid_func = reinterpret_cast<double(*)(const vec3d &)>(dylib.load_symbol("solid"));
		auto set_boundary_flux = reinterpret_cast<void(*)( double, Real [DIM3][2] )>(dylib.load_symbol("set_boundary_flux"));
		//
		if( solid ) {
			solid->clear(1.0);
			if( solid_func || moving_solid_func ) {
				//
				solid->parallel_all([&](int i, int j, int k, auto &it) {
					const vec3d &p = m_dx*vec3i(i,j,k).nodal();
					double value (1.0);
					if( solid_func ) value = std::min(value,solid_func(p));
					if( moving_solid_func ) value = std::min(value,moving_solid_func(time,p).first);
					if( ! m_param.narrowband_dist || std::abs(value) < m_param.narrowband_dist*m_dx ) it.set(value);
				});
				//
				if( m_param.narrowband_dist ) {
					solid->set_as_levelset(m_param.narrowband_dist*m_dx);
					solid->flood_fill();
				}
				//
				if( solid_velocity && moving_solid_func ) {
					//
					shared_array3<vec3d> tmp_nodal_velocity(m_shape.nodal());
					tmp_nodal_velocity->parallel_all([&](int i, int j, int k, auto &it) {
						if( (*solid)(i,j,k) < 0.0 ) {
							const auto info = moving_solid_func(time,m_dx*vec3i(i,j,k).nodal());
							it.set(info.second);
						}
					});
					//
					const auto get_face_vertices = [&]( int dim, const vec3i &pi ) {
						std::vector<vec3i> result(4);
						int i(pi[0]), j(pi[1]), k(pi[2]);
						if( dim == 0 ) {
							result[0] = vec3i(i,j,k);
							result[1] = vec3i(i,j+1,k);
							result[2] = vec3i(i,j+1,k+1);
							result[3] = vec3i(i,j,k+1);
						} else if( dim == 1 ) {
							result[0] = vec3i(i,j,k);
							result[1] = vec3i(i+1,j,k);
							result[2] = vec3i(i+1,j,k+1);
							result[3] = vec3i(i,j,k+1);
						} else if( dim == 2 ) {
							result[0] = vec3i(i,j,k);
							result[1] = vec3i(i+1,j,k);
							result[2] = vec3i(i+1,j+1,k);
							result[3] = vec3i(i,j+1,k);
						}
						return result;
					};
					//
					solid_velocity->clear();
					solid_velocity->parallel_all([&]( int dim, int i, int j, int k, auto &it ) {
						std::vector<vec3i> results = get_face_vertices(dim,vec3i(i,j,k));
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
			Real flux[DIM3][2] = {{0.0,0.0},{0.0,0.0},{0.0,0.0}};
			set_boundary_flux(time,flux);
			for( int j=0; j<m_shape.h; ++j ) for( int k=0; k<m_shape.d; ++k ) {
				(*solid_velocity)[0].set(0,j,k,flux[0][0]);
				(*solid_velocity)[0].set(m_shape.w,j,k,flux[0][1]);
			}
			for( int i=0; i<m_shape.w; ++i ) for( int k=0; k<m_shape.d; ++k ) {
				(*solid_velocity)[1].set(i,0,k,flux[1][0]);
				(*solid_velocity)[1].set(i,m_shape.h,k,flux[1][1]);
			}
			for( int i=0; i<m_shape.w; ++i ) for( int j=0; j<m_shape.h; ++j ) {
				(*solid_velocity)[2].set(i,j,0,flux[2][0]);
				(*solid_velocity)[2].set(i,j,m_shape.d,flux[2][1]);
			}
		}
		//
		auto velocity_func = reinterpret_cast<vec3d(*)(double, const vec3d &)>(dylib.load_symbol("boundary_velocity"));
		if( velocity_func && solid_velocity ) {
			solid_velocity->parallel_all([&]( int dim, int i, int j, int k, auto &it ) {
				const double value = velocity_func(time,m_dx*vec3i(i,j,k).face(dim))[dim];
				if( value ) it.set(value);
				else it.set_off();
			});
		}
	}
	virtual void add_force( vec3d p, vec3d f, macarray3<Real> &external_force ) const override {
		for( unsigned dim : DIMS3 ) {
			vec3d index_coord = p/m_dx-vec3d(0.5,0.5,0.5);
			external_force[dim].set(m_shape.face(dim).clamp(index_coord),f[dim]);
		}
	}
	//
	virtual void configure( configuration &config ) override {
		//
		config.get_double("EpsFluid",m_param.eps_fluid,"Minimal bound for fluid fraction");
		config.get_double("EpsSolid",m_param.eps_solid,"Minimal bound for solid fraction");
		config.get_bool("WENO",m_param.weno_interpolation,"Whether to use WENO interpolation");
		config.get_double("NarrowBandDist",m_param.narrowband_dist);
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
	struct Parameters {
		//
		double eps_fluid {1e-2};
		double eps_solid {1e-2};
		bool weno_interpolation {false};
		double narrowband_dist {sqrt(3.0)};
	};
	Parameters m_param;
	//
	double m_dx;
	shape3 m_shape;
	parallel_driver m_parallel{this};
	cellmesher3_driver m_mesher{this,"marchingcubes"};
};
//
extern "C" module * create_instance() {
	return new macutility3();
}
//
extern "C" const char *license() {
	return "MIT";
}
//