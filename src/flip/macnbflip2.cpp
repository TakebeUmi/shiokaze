/*
**	macnbflip2.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on March 29, 2017.
**	APIC extension by Takahiro Sato.
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
#include "macnbflip2.h"
#include <shiokaze/array/shared_array2.h>
#include <shiokaze/array/shared_bitarray2.h>
#include <shiokaze/array/array_utility2.h>
#include <shiokaze/array/array_interpolator2.h>
#include <shiokaze/array/macarray_interpolator2.h>
#include <shiokaze/array/array_derivative2.h>
#include <shiokaze/graphics/graphics_utility.h>
#include <shiokaze/utility/utility.h>
#include <algorithm>
//
SHKZ_USING_NAMESPACE
//
static const double default_mass = 1.0 / 4.0;
//
double macnbflip2::grid_kernel( const vec2d &r, double dx ) {
	//
	const double x = std::abs(r[0]) / dx;
	const double y = std::abs(r[1]) / dx;
	return std::max(0.0,1.0-x) * std::max(0.0,1.0-y);
}
//
vec2d macnbflip2::grid_gradient_kernel( const vec2d &r, double dx ) {
	//
	const double x = std::abs(r[0]) / dx;
	const double y = std::abs(r[1]) / dx;
	if( x <= 1.0 && y <= 1.0 ) {
		const double u = std::copysign(1.0-y,r[0]);
		const double v = std::copysign(1.0-x,r[1]);
		return vec2d(u,v) / dx;
	} else {
		return vec2d();
	}
}
//
void macnbflip2::configure( configuration &config ) {
	//
	config.get_bool("APIC",m_param.use_apic,"Whether to use APIC");
	config.get_unsigned("Narrowband",m_param.narrowband,"Narrowband bandwidth");
	config.get_double("FitParticleDist",m_param.fit_particle_dist,"FLIP particle fitting threshold");
	config.get_integer("RK_Order",m_param.RK_order,"Order of accuracy for Runge-kutta integration");
	config.get_double("Erosion",m_param.erosion,"Rate of erosion for internal levelset");
	config.get_double("MinMassPerCell",m_param.min_mass_per_cell,"Minimal target mass per cell");
	config.get_double("MaxMassPerCell",m_param.max_mass_per_cell,"Maximal target mass per cell");
	config.get_unsigned("MiminalLiveCount",m_param.minimal_live_count,"Minimal step of particles to stay alive");
	config.get_double("CorrectStiff",m_param.stiff,"Position correction strength");
	config.get_bool("VelocityCorrection",m_param.velocity_correction,"Should do velocity correction");
	config.get_double("BulletMaximalTime",m_param.bullet_maximal_time,"Maximal time for bullet particles to scale smaller");
	config.get_bool("SplatBulletParticle",m_param.splat_bullet_particle,"Splat bullet particles");
	config.get_bool("DrawFLIPParticles",m_param.draw_particles,"Whether to draw FLIP particles.");
	config.get_double("DecayRate",m_param.decay_rate,"Decay rate for tracer particles");
	config.get_bool("CollisionDomainBoundary",m_param.collision_domain_boundary,"Collision domain boundary");
	config.get_bool("RescaleGradient",m_param.rescale_gradient,"Rescale gradient");
}
//
void macnbflip2::initialize( const shape2 &shape, double dx ) {
	//
	m_shape = shape;
	m_dx = dx;
}
//
void macnbflip2::post_initialize( bool initialized_from_file ) {
	//
	if( ! initialized_from_file ) m_particles.resize(0);
	sort_particles();
}
//
void macnbflip2::splat( double time, macarray2<macflip2_interface::mass_momentum2> &mass_and_momentum ) const {
	//
	if( m_particles.size()) {
		//
		mass_and_momentum.clear();
		for( size_t n=0; n<m_particles.size(); ++n ) {
			if( m_param.splat_bullet_particle || ! m_particles[n].bullet ) {
				for( int dim : DIMS2 ) {
					const vec2d p = m_particles[n].p;
					const vec2i pi = vec2i(p[0]/m_dx-0.5*(dim!=0),p[1]/m_dx-0.5*(dim!=1));
					for( int ii=0; ii<=1; ++ii ) for( int jj=0; jj<=1; ++jj ) {
						mass_and_momentum[dim].set(m_shape.clamp(pi+vec2i(ii,jj)),{0.0,0.0});
					}
				}
			}
		}
		mass_and_momentum.parallel_actives([&]( int dim, int i, int j, auto &it, int tn ) {
			//
			Real mom (0.0), m (0.0);
			vec2d pos = m_dx*vec2i(i,j).face(dim);
			std::vector<size_t> neighbors = m_pointgridhash->get_face_neighbors(vec2i(i,j),dim);
			for( size_t k : neighbors ) {
				if( m_param.splat_bullet_particle || ! m_particles[k].bullet ) {
					const Particle &p = m_particles[k];
					const double w = grid_kernel(p.p-pos,m_dx);
					if( w ) {
						mom += w*p.mass*p.velocity[dim];
						m += w*p.mass;
						if( m_param.use_apic ) {
							const vec2d &c = (p.c)[dim];
							const vec2d r = pos-p.p;
							mom += w*p.mass*(c*r);
						}
					}
				}
			}
			if( m ) it.set({m,mom});
			else it.set_off();
		});
		//
	} else {
		mass_and_momentum.clear();
	}
}
//
void macnbflip2::advect( std::function<double(const vec2d &p)> solid,
						 std::function<vec2d(const vec2d &p)> velocity,
						 double time, double dt ) {
	//
	if( m_particles.size()) {
		//
		m_parallel.for_each(m_particles.size(),[&]( size_t n ) {
			Particle &particle = m_particles[n];
			bool bullet = particle.bullet;
			const vec2d &u = particle.velocity;
			const vec2d &p = particle.p;
			if( bullet ) {
				particle.p += dt*u;
			} else {
				vec2d u1 = velocity(p);
				if( u1.norm2()) {
					if( m_param.RK_order==4 ) {
						vec2d u2 = velocity(p+0.5*dt*u1);
						vec2d u3 = velocity(p+0.5*dt*u2);
						vec2d u4 = velocity(p+dt*u3);
						particle.p += dt*(u1+2.0*u2+2.0*u3+u4)/6.0;
					} else if( m_param.RK_order==2 ) {
						vec2d u2 = velocity(p+dt*u1);
						particle.p += dt*0.5*(u1+u2);
					} else if( m_param.RK_order==1 ) {
						particle.p += dt*(u1);
					} else {
						printf( "Unsupported RK order (%d)\n", m_param.RK_order );
						exit(0);
					}
				} else {
					particle.p += dt*u;
				}
			}
			// Decay particle sizing value
			m_particles[n].sizing_value = std::max(0.0,m_particles[n].sizing_value-m_param.decay_rate*dt);
		});
		//
		sort_particles();
	}
	//
	// Perform collision
	collision(solid);
}
//
void macnbflip2::mark_bullet( std::function<double(const vec2d &p)> fluid, std::function<vec2d(const vec2d &p)> velocity, double time ) {
	//
	if( m_particles.size()) {
		//
		// Mark bullet
		m_parallel.for_each(m_particles.size(),[&]( size_t n ) {
			Particle &particle = m_particles[n];
			char new_status (0);
			if( fluid(particle.p) > 0.0 ) {
				new_status = 1;
				for( int dim : DIMS2 ) particle.c[dim] = vec2d();
			}
			if( new_status != particle.bullet ) {
				particle.bullet = new_status;
				particle.bullet_time = new_status ? time : 0.0;
				if( new_status == 0 ) {
					particle.velocity = velocity(particle.p);
				}
			}
		});
		//
		// Remove long-term ballistic particles
		std::vector<char> remove_flag(m_particles.size(),0);
		m_parallel.for_each(m_particles.size(),[&]( size_t n ) {
			Particle &particle = m_particles[n];
			if( particle.bullet ) {
				if( time-particle.bullet_time > m_param.bullet_maximal_time ) {
					remove_flag[n] = 1;
				} else {
					const double scale = std::max(0.01,1.0 - std::max(0.0,time-particle.bullet_time) / m_param.bullet_maximal_time);
					particle.r = 0.25 * m_dx * scale;
					particle.mass = (scale*scale) * default_mass;
				}
			}
		});
		//
		// Reconstruct a new particle array
		size_t removed_total (0);
		for( size_t i=0; i<m_particles.size(); ++i) {
			if(  remove_flag[i] ) ++ removed_total;
		}
		//
		if( removed_total ) {
			//
			std::vector<Particle> old_particles (m_particles);
			m_particles.clear();
			for( size_t i=0; i<old_particles.size(); ++i) {
				if( ! remove_flag[i] ) m_particles.push_back(old_particles[i]);
			}
			m_particles.shrink_to_fit();
			sort_particles();
		}
	}
}
//
void macnbflip2::update( std::function<double(const vec2d &p)> solid, array2<Real> &fluid, double time, bool add_active ) {
	//
	if( m_particles.size()) {
		//
		shared_array2<Real> save_fluid (fluid);
		fluid.parallel_actives([&]( int i, int j, auto &it ) {
			if( solid(m_dx*vec2d(i,j)) > m_dx ) {
				it.increment(m_param.erosion*m_dx);
			}
		});
		//
		shared_bitarray2 mask(fluid.shape());
		std::vector<particlerasterizer2_interface::Particle2> points;
		for( int n=0; n<m_particles.size(); ++n ) {
			particlerasterizer2_interface::Particle2 point;
			vec2d p = m_particles[n].p;
			point.p = p;
			point.r = m_particles[n].r;
			points.push_back(point);
			const vec2i pi = mask->shape().clamp(p/m_dx);
			if( add_active || fluid.active(pi)) {
				mask->set(pi);
			}
		}
		mask->dilate(2);
		if( add_active ) fluid.activate_as_bit(mask());
		shared_array2<Real> particle_levelset(m_shape,0.125*m_dx);
		m_particlerasterizer->build_levelset(particle_levelset(),mask(),points);
		//
		// Safety check
		particle_levelset->parallel_actives([&]( int i, int j, auto &it, int tn ) {
			bool edge_flag (false);
			for( int dim : DIMS2 ) for( int dir=-1; dir<=1; dir+=2 ) {
				const vec2i pi = vec2i(i,j)+dir*vec2i(dim==0,dim==1);
				if( pi[dim] >= 0 && pi[dim] < m_shape[dim] && ! particle_levelset->active(pi)) {
					edge_flag = true;
				}
			}
			if( edge_flag ) it.set(std::max((Real)(0.125*m_dx),it()));
		});
		//
		if( m_param.rescale_gradient ) {
			particle_levelset->parallel_actives([&](int i, int j, auto &it, int tn) {
				it.multiply(m_gridutility->get_upwind_levelset_gradient(fluid,vec2i(i,j)).len());
			});
		}
		//
		fluid.parallel_actives([&](int i, int j, auto &it, int tn) {
			//
			const double fluid_value (it());
			const double save_value (save_fluid()(i,j));
			//
			double sizing_sum (0.0);
			double sizing_weight (0.0);
			//
			std::vector<size_t> list = m_pointgridhash->get_points_in_cell(vec2i(i,j));
			for( const auto &n : list ) {
				double w (m_particles[n].mass);
				sizing_sum += w*std::min((Real)1.0,m_particles[n].sizing_value);
				sizing_weight += w;
			}
			if( sizing_weight ) sizing_sum /= sizing_weight;
			//
			double value = sizing_sum * particle_levelset()(i,j) + (1.0-sizing_sum) * save_value;
			for( const auto &n : list ) {
				const double theta = std::min((Real)1.0,m_particles[n].sizing_value);
				const double particle_value = (m_particles[n].p-m_dx*vec2i(i,j).cell()).len()-m_particles[n].r;
				value = std::min(value,theta*particle_value+(1.0-theta)*save_value);
			}
			it.set(std::min(fluid_value,value));
		});
	}
}
//
void macnbflip2::collision( std::function<double(const vec2d &p)> solid ) {
	//
	if( m_particles.size()) {
		//
		m_parallel.for_each( m_particles.size(), [&]( size_t pindex ) {
			Particle &particle = m_particles[pindex];
			vec2r &p = particle.p;
			vec2r &u = particle.velocity;
			const double &r = particle.r;
			double phi = solid(p)-r;
			if( phi < 0.0 ) {
				vec2d gradient = interpolate_solid_gradient(solid,p);
				p = p - phi*gradient;
				double dot = gradient * u;
				if( dot < 0.0 ) u = u - gradient*dot;
			}
			if( m_param.collision_domain_boundary ) {
				for( unsigned dim : DIMS2 ) {
					if( p[dim] < r ) {
						p[dim] = r;
						if( u[dim] < 0.0 ) u[dim] = 0.0;
					}
					if( p[dim] > m_dx*m_shape[dim]-r ) {
						p[dim] = m_dx*m_shape[dim]-r;
						if( u[dim] > 0.0 ) u[dim] = 0.0;
					}
				}
			}
		});
	}
	//
	sort_particles();
}
//
void macnbflip2::correct( std::function<double(const vec2d &p)> fluid, const macarray2<Real> &velocity ) {
	//
	if( m_particles.size()) {
		//
		// Compute the displacements
		std::vector<vec2d> displacements(m_particles.size());
		m_parallel.for_each(m_particles.size(),[&]( size_t n ) {
			const Particle &pi = m_particles[n];
			vec2d displacement;
			std::vector<size_t> neighbors = m_pointgridhash->get_cell_neighbors(m_shape.find_cell(pi.p/m_dx),pointgridhash2_interface::USE_NODAL);
			for( const size_t &j : neighbors ) {
				if( n != j ) {
					const Particle &pj = m_particles[j];
					double dist2 = (pi.p-pj.p).norm2();
					double target = pi.r+pj.r;
					if( dist2 < target*target ) {
						double diff = target-sqrt(dist2);
						const double &mi = pi.mass;
						const double &mj = pj.mass;
						displacement += m_param.stiff * diff * (pi.p-pj.p).normal() * mj / (mi+mj);
					}
				}
			}
			displacements[n] = displacement;
		});
		//
		// Kill the normal components to prevent volume gain
		m_parallel.for_each(m_particles.size(),[&]( size_t n ) {
			if( ! displacements[n].empty()) {
				const Particle &p = m_particles[n];
				vec2d new_pos = p.p+displacements[n];
				vec2d normal = this->interpolate_fluid_gradient(fluid,new_pos);
				double dot = displacements[n] * normal;
				if( dot > 0.0 ) displacements[n] -= dot * normal;
			}
		});
		//
		// Apply velocity correction
		if( m_param.velocity_correction ) {
			m_parallel.for_each(m_particles.size(),[&]( size_t n ) {
				vec2d mid_p = 0.5 * (m_particles[n].p+displacements[n]);
				for( int dim : DIMS2 ) {
					vec2d gradient = array_derivative2::derivative(velocity[dim],m_dx*vec2i().face(dim),m_dx,mid_p);
					m_particles[n].velocity[dim] += gradient * displacements[n];
				}
			});
		}
		//
		// Apply displacement
		m_parallel.for_each(m_particles.size(),[&]( size_t n ) { m_particles[n].p += displacements[n]; });
		//
		// Update hash table
		sort_particles();
	}
}
//
void macnbflip2::fit_particle( std::function<double(const vec2d &p)> fluid, Particle &particle, const vec2d &gradient ) const {
	//
	if( std::abs(fluid(particle.p)) < m_param.fit_particle_dist * particle.r ) {
		for( unsigned n=0; n<3; ++n ) {
			particle.p -= 0.5 * (fluid(particle.p)+particle.r) * gradient;
		}
	}
}
//
size_t macnbflip2::resample(const array2<Real> &fluid,
						std::function<double(const vec2d &p)> solid,
						const macarray2<Real> &velocity,
						std::function<bool(const vec2d &p)> mask ) {
	//
	// Compute narrowband
	shared_bitarray2 narrowband_mask(m_shape);
	bool smoke_simulation_flag (true);
	if( ! mask ) mask = []( const vec2d &p ) { return true; };
	//
	if( m_param.narrowband ) {
		fluid.const_serial_actives([&](int i, int j, const auto &it) {
			const vec2d p = m_dx*vec2i(i,j).cell();
			if( it() > 0 && solid(p) > 0 ) narrowband_mask->set(i,j);
			if( it() < 0.0 ) smoke_simulation_flag = false;
		});
		//
		if( smoke_simulation_flag ) {
			narrowband_mask->activate_all();
		} else {
			narrowband_mask->dilate([&](int i, int j, auto& it, int tn ) {
				const vec2d p = m_dx*vec2i(i,j).cell();
				if(fluid.active(i,j) && fluid(i,j) < 0 && solid(p) > 0.125*m_dx && mask(p)) it.set();
			},m_param.narrowband);
			//
			fluid.const_serial_actives([&](int i, int j, auto &it) {
				if( it() > m_dx ) narrowband_mask->set_off(i,j);
			});
		}
	} else {
		fluid.const_serial_actives([&](int i, int j, const auto &it) {
			const vec2d p = m_dx*vec2i(i,j).cell();
			if( solid(p) > 0 ) narrowband_mask->set(i,j);
			if( it() < 0.0 ) smoke_simulation_flag = false;
		});
		fluid.const_serial_inside([&](int i, int j, const auto &it) {
			const vec2d p = m_dx*vec2i(i,j).cell();
			if( solid(p) > 0 ) narrowband_mask->set(i,j);
			if( it() < 0.0 ) smoke_simulation_flag = false;
		});
		if( smoke_simulation_flag ) narrowband_mask->activate_all();
	}
	//
	// Compute sizing function
	shared_array2<Real> sizing_array(m_shape);
	compute_sizing_func(fluid,narrowband_mask(),velocity,sizing_array());
	//
	// Update sizing value on particles
	if( m_particles.size()) {
		m_parallel.for_each(m_particles.size(),[&]( size_t n, int tn ) {
			Particle &p = m_particles[n];
			if( ! p.bullet ) {
				p.sizing_value = std::max((Real)p.sizing_value,sizing_array()(m_shape.find_cell(p.p/m_dx)));
			}
		});
	}
	//
	std::vector<std::vector<Particle> > new_particles_t(m_parallel.get_thread_num());
	std::vector<char> remove_particles(m_particles.size(),0);
	//
	// Bucket cell method to remove too dense particles
	shared_array2<float> cell_bucket(m_shape);
	if( m_particles.size()) {
		for( size_t n=0; n<m_particles.size(); ++n ) {
			const Particle &p = m_particles[n];
			if( mask(p.p)) {
				vec2i pi = m_shape.clamp(p.p/m_dx);
				int i (pi[0]), j (pi[1]);
				m_shape.clamp(i,j);
				if( ! p.bullet ) {
					if( ! sizing_array()(i,j) || cell_bucket()(i,j) > m_param.max_mass_per_cell || p.sizing_value <= 0.0 ) {
						if( p.live_count > m_param.minimal_live_count ) {
							remove_particles[n] = 1;
						}
					}
					if( ! narrowband_mask()(i,j) ) {
						remove_particles[n] = 1;
					}
				}
				if( ! remove_particles[n] && solid(p.p) < -p.r ) {
					remove_particles[n] = 1;
				}
				if( ! m_param.collision_domain_boundary ) {
					if( utility::box(p.p,vec2d(0.0,0.0),m_dx*vec2d(m_shape[0],m_shape[1])) > -p.r ) {
						remove_particles[n] = 1;
					}
				}
				if( ! remove_particles[n] ) {
					cell_bucket().increment(i,j,p.mass);
				}
			}
		}
		//
		// Increment live count
		for( size_t n=0; n<m_particles.size(); ++n ) {
			m_particles[n].live_count ++;
		}
	}
	//
	// Particle reseeding...
	narrowband_mask->const_parallel_actives([&]( int i, int j, int tn ) {
		//
		double mass_added (0.0);
		double sizing_value = sizing_array()(i,j);
		auto attempt_resample = [&]( const vec2d &p ) {
			//
			double r = 0.25 * m_dx;
			if( mask(p) && cell_bucket()(i,j)+mass_added < m_param.min_mass_per_cell ) {
				//
				if( smoke_simulation_flag || this->interpolate_fluid(fluid,p) < -r ) {
					//
					// Put a FLIP particle here...
					Particle new_particle;
					new_particle.p = p;
					bool sparse (true);
					const std::vector<size_t> &indices = m_pointgridhash->get_points_in_cell(vec2i(i,j));
					for( auto it=indices.begin(); it!=indices.end(); ++it ) {
						if( (m_particles[*it].p-p).len() <= 2.0*r ) {
							sparse = false;
							break;
						}
					}
					if( sparse && solid(new_particle.p) > r ) {
						new_particle.mass = default_mass;
						new_particle.velocity = macarray_interpolator2::interpolate(velocity,vec2d(),m_dx,p,true);
						new_particle.r = r;
						new_particle.bullet = 0;
						new_particle.bullet_time = 0.0;
						new_particle.live_count = 0;
						new_particle.sizing_value = sizing_value;
						if( m_param.use_apic ) this->update_velocity_derivative(new_particle,velocity);
						this->fit_particle([&](const vec2d &p) {
							return this->interpolate_fluid(fluid,p);
						},new_particle,this->interpolate_fluid_gradient(fluid,new_particle.p));
						new_particles_t[tn].push_back(new_particle);
						mass_added += new_particle.mass;
					}
				}
			}
		};
		//
		if( sizing_value ) {
			if( m_param.narrowband && ! smoke_simulation_flag && fluid(i,j) < -1.25*m_dx ) {
				attempt_resample(m_dx*vec2i(i,j).cell());
			} else {
				for( unsigned ii=0; ii<2; ii++ ) for( unsigned jj=0; jj<2; jj++ ) {
					vec2d p = m_dx*vec2d(i,j)+0.25*m_dx*vec2d(1,1)+0.5*m_dx*vec2d(ii,jj);
					attempt_resample(p);
				}
			}
		}
	});
	//
	// Reconstruct a new particle array
	std::vector<Particle> old_particles(m_particles);
	m_particles.clear();
	size_t reseeded (0);
	for( size_t t=0; t<new_particles_t.size(); ++t ) {
		m_particles.insert(m_particles.end(),new_particles_t[t].begin(),new_particles_t[t].end());
		reseeded += new_particles_t[t].size();
	}
	for( size_t i=0; i<remove_particles.size(); i++) {
		if( ! remove_particles[i] ) m_particles.push_back(old_particles[i]);
	}
	m_particles.shrink_to_fit();
	//
	// Update hash table
	sort_particles();
	return reseeded;
}
//
size_t macnbflip2::remove(std::function<bool(const vec2r &p, bool bullet)> test_function ) {
	//
	size_t removed_count (0);
	if( m_particles.size()) {
		//
		std::vector<Particle> old_particles (m_particles);
		std::vector<char> remove_flag (m_particles.size(),0);
		//
		m_parallel.for_each(m_particles.size(),[&]( size_t n ) {
			remove_flag[n] = test_function(old_particles[n].p,old_particles[n].bullet);
		});
		//
		m_particles.clear();
		for( size_t i=0; i<old_particles.size(); ++i) {
			if( remove_flag[i] ) ++ removed_count;
			else m_particles.push_back(old_particles[i]);
		}
		//
		if( removed_count ) {
			m_particles.shrink_to_fit();
			sort_particles();
		}
	}
	return removed_count;
}
//
void macnbflip2::update( const macarray2<Real> &prev_velocity, const macarray2<Real> &new_velocity,
						 double dt, vec2d gravity, double PICFLIP ) {
	//
	if(m_particles.size()) {
		m_parallel.for_each(m_particles.size(),[&]( size_t n ) {
			Particle &particle = m_particles[n];
			if( particle.bullet ) {
				//
				// If the particle is a ballistic particle, just add gravity force
				particle.velocity += dt * gravity;
				//
			} else {
				if( m_param.use_apic ) {
					//
					// Fetch grid velocity
					particle.velocity = macarray_interpolator2::interpolate(new_velocity,vec2d(),m_dx,particle.p,true);
					//
					// Update particle velocity derivative (APIC)
					update_velocity_derivative(particle,new_velocity);
					//
				} else {
					//
					// Fetch grid velocity
					vec2d new_grid_velocity = macarray_interpolator2::interpolate(new_velocity,vec2d(),m_dx,particle.p,true);
					vec2d old_grid_velocity = macarray_interpolator2::interpolate(prev_velocity,vec2d(),m_dx,particle.p,true);
					//
					// Compute pure FLIP velocity
					vec2d FLIP_velocity = particle.velocity + (new_grid_velocity-old_grid_velocity);
					//
					// Compute PICFLIP velocity
					vec2d PICFLIP_velocity = PICFLIP * FLIP_velocity + (1.0-PICFLIP) * new_grid_velocity;
					//
					// Compute the final velocity of FLIP
					particle.velocity = PICFLIP_velocity;
				}
			}
		});
	}
}
//
void macnbflip2::update( std::function<void(const vec2r &p, vec2r &velocity, Real &mass, bool bullet )> func ) {
	m_parallel.for_each(m_particles.size(),[&]( size_t n ) {
		Particle &particle = m_particles[n];
		func(particle.p,particle.velocity,particle.mass,particle.bullet);
	});
}
//
std::vector<macflip2_interface::particle2> macnbflip2::get_particles() const {
	std::vector<macflip2_interface::particle2> result;
	for( size_t n=0; n<m_particles.size(); ++n ) {
		macflip2_interface::particle2 particle;
		particle.p = m_particles[n].p;
		particle.r = m_particles[n].r;
		particle.sizing_value = m_particles[n].sizing_value;
		particle.bullet = m_particles[n].bullet;
		particle.bullet_time = m_particles[n].bullet_time;
		result.push_back(particle);
	}
	return result;
}
//
void macnbflip2::sort_particles() {
	//
	if( m_particles.size()) {
		std::vector<vec2r> points(m_particles.size());
		m_parallel.for_each(m_particles.size(),[&]( size_t n) {
			points[n] = m_particles[n].p;
		});
		m_pointgridhash->sort_points(points);
	} else {
		m_pointgridhash->clear();
	}
}
//
void macnbflip2::update_velocity_derivative( Particle& particle, const macarray2<Real> &velocity ) {
	//
	// Written by Takahiro Sato updated by Ryoichi Ando
	for( int dim : DIMS2 ) {
		vec2r &c = particle.c[dim];
		c = vec2d();
		//
		// particle position
		const vec2d p_pos = particle.p;
		vec2d offset = 0.5*vec2d(dim!=0,dim!=1);
		//
		// cell index
		int i = std::floor(p_pos[0]/m_dx-offset[0]);
		int j = std::floor(p_pos[1]/m_dx-offset[1]);
		//
		// cell position
		vec2d cell_pos[4]; // 0: lower left, 1: lower right, 2: upper left, 3: upper right
		cell_pos[0] = m_dx*(vec2d(i,j)+offset);
		cell_pos[1] = m_dx*(vec2d(i+1,j)+offset);
		cell_pos[2] = m_dx*(vec2d(i,j+1)+offset);
		cell_pos[3] = m_dx*(vec2d(i+1,j+1)+offset);
		//
		bool actives[4];
		const auto v_shape = velocity[dim].shape();
		actives[0] = velocity[dim].active(v_shape.clamp(i,j));
		actives[1] = velocity[dim].active(v_shape.clamp(i+1,j));
		actives[2] = velocity[dim].active(v_shape.clamp(i,j+1));
		actives[3] = velocity[dim].active(v_shape.clamp(i+1,j+1));
		//
		bool APIC_enabled (true);
		for( char n=0; n<4; ++n ) {
			if( ! actives[n] ) {
				APIC_enabled = false;
				break;
			}
		}
		//
		if( APIC_enabled ) {
			vec2d dw[4];
			dw[0] = grid_gradient_kernel(cell_pos[0]-p_pos, m_dx);
			dw[1] = grid_gradient_kernel(cell_pos[1]-p_pos, m_dx);
			dw[2] = grid_gradient_kernel(cell_pos[2]-p_pos, m_dx);
			dw[3] = grid_gradient_kernel(cell_pos[3]-p_pos, m_dx);
			//
			Real u[4];
			u[0] = velocity[dim](v_shape.clamp(i,j));
			u[1] = velocity[dim](v_shape.clamp(i+1,j));
			u[2] = velocity[dim](v_shape.clamp(i,j+1));
			u[3] = velocity[dim](v_shape.clamp(i+1,j+1));
			//
			for( char n=0; n<4; ++n ) {
				c += dw[n] * u[n];
			}
		}
	}
}
//
double macnbflip2::interpolate_fluid( const array2<Real> &fluid, const vec2d &p ) const {
	return array_interpolator2::interpolate(fluid,p/m_dx-vec2d(0.5,0.5));
}
//
vec2d macnbflip2::interpolate_fluid_gradient( std::function<double(const vec2d &p)> fluid, const vec2d &p ) const {
	//
	vec2d result;
	for( int dim : DIMS2 ) result[dim] = fluid(p+0.25*m_dx*vec2d(dim==0,dim==1))-fluid(p-0.25*m_dx*vec2d(dim==0,dim==1));
	return result.normal();
}
//
vec2d macnbflip2::interpolate_fluid_gradient( const array2<Real> &fluid, const vec2d &p ) const {
	//
	Real derivative[DIM2];
	array_derivative2::derivative(fluid,p/m_dx-vec2d(0.5,0.5),derivative);
	return vec2d(derivative).normal();
}
//
vec2d macnbflip2::interpolate_solid_gradient( std::function<double(const vec2d &p)> solid, const vec2d &p ) const {
	//
	vec2d result;
	for( int dim : DIMS2 ) result[dim] = solid(p+0.5*m_dx*vec2d(dim==0,dim==1))-solid(p-0.5*m_dx*vec2d(dim==0,dim==1));
	return result.normal();
}
//
void macnbflip2::draw_flip_circle ( graphics_engine &g, const vec2d &p, double r, bool bullet, double sizing_value ) const {
	//
	const unsigned num_v = 10;
	if( bullet ) {
		g.color4(1.0,0.5,0.5,sizing_value);
	} else {
		g.color4(0.5,0.5,1.0,sizing_value);
	}
	//
	g.begin(graphics_engine::MODE::TRIANGLE_FAN);
	for( unsigned t=0; t<num_v; t++ ) {
		double theta = 2.0 * M_PI * t / (double)num_v;
		g.vertex2v((p+r*vec2d(cos(theta),sin(theta))).v);
	}
	g.end();
	g.color4(1.0,1.0,1.0,0.5);
	g.begin(graphics_engine::MODE::LINE_LOOP);
	for( int t=0; t<num_v; t++ ) {
		double theta = 2.0 * M_PI * t / (double)num_v;
		g.vertex2v((p+r*vec2d(cos(theta),sin(theta))).v);
	}
	g.end();
}
//
void macnbflip2::draw( graphics_engine &g, double time ) const {
	//
	if( m_param.draw_particles ) {
		for( const Particle &particle : m_particles ) {
			const vec2d &p = particle.p;
			draw_flip_circle(g,p,particle.r,particle.bullet,particle.sizing_value);
			graphics_utility::draw_arrow(g,p.v,(p+m_dx*particle.velocity).v);
		}
	}
}
//
extern "C" module * create_instance() {
	return new macnbflip2();
}
//
extern "C" const char *license() {
	return "MIT";
}
//
