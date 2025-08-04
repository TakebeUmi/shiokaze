/*
**	macoctreeliquid2.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on October 9, 2018.All rights reserved.
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
#include "macoctreeliquid2.h"
#include <shiokaze/graphics/graphics_utility.h>
#include <shiokaze/array/array_interpolator2.h>
#include <shiokaze/array/macarray_interpolator2.h>
#include <shiokaze/array/array_derivative2.h>
#include <shiokaze/array/shared_array2.h>
#include <shiokaze/utility/utility.h>
#include <shiokaze/array/array_extrapolator2.h>
#include <shiokaze/array/macarray_extrapolator2.h>
#include <shiokaze/core/filesystem.h>
#include "fastapprox/fastexp.h"
#include <cmath>
#include <limits>
#include <mutex>
#include <thread>
#include <numeric>
//
SHKZ_USING_NAMESPACE
using namespace macotreeliquid2_namespace;
//
void macoctreeliquid2::load( configuration &config ) {
	//
	std::string name("waterdrop2"); config.get_string("Name",name,"Scene file name");
	m_dylib.open_library(filesystem::resolve_libname(name));
	m_dylib.load(config);
	m_dylib.overwrite(config);
}
//
void macoctreeliquid2::configure( configuration &config ) {
	//
	m_dylib.configure(config);
	//
	m_shape = shape2(64,32);
	//
	config.set_default_string("LinSolver","amg");
	//
	config.get_vec2d("Gravity",m_param.gravity.v,"Gravity vector");
	config.get_bool("UseFLIP",m_param.use_FLIP,"Whether to use FLIP");
	config.get_unsigned("MinResolution",m_param.min_resolution,"Minimal resolution");
	config.get_unsigned("ResolutionX",m_shape[0],"Resolution towards X axis");
	config.get_unsigned("ResolutionY",m_shape[1],"Resolution towards Y axis");
	config.get_double("PICFLIP",m_param.PICFLIP,"PICFLIP blending factor");
	config.get_bool("VolumeCorrection",m_param.volume_correction,"Whether to perform volume correction");
	config.get_bool("RegionalVolumeCorrection",m_param.regional_volume_correction,"Regional volume correction");
	config.get_bool("MacCormack",m_param.maccormack,"Use MacCormack method");
	config.get_unsigned("ErodeWidth",m_param.erode_width,"Erosion width");
	config.get_double("SurfaceTension",m_param.surftens_k,"Surface tension coefficient");
	config.get_bool("UseSizingFunc",m_param.use_sizing_func,"Use sizing function");
	config.get_unsigned("InitialRefinement",m_param.initial_refinement,"Initial refinement count");
	config.get_double("MaxCFLAccumulation",m_param.maximal_CFL_accumulation,"CFL sum trigger for remeshing");
	//
	if( m_param.use_sizing_func ) {
		config.set_default_bool("SteepAdaptivity",false);
		config.set_default_unsigned("DilateCount",2);
	}
	//
	double scale (1.0);
	config.get_double("ResolutionScale",scale,"Resolution doubling scale");
	//
	double view_scale (1.0);
	config.get_double("ViewScale",view_scale,"View scale");
	//
	double resolution_scale (1.0);
	config.get_double("ResolutionScale",resolution_scale,"Resolution doubling scale");
	//
	m_shape *= resolution_scale;
	m_dx = view_scale * m_shape.dx();
}
//
void macoctreeliquid2::post_initialize( bool initialized_from_file ) {
	//
	auto initialize_func = reinterpret_cast<void(*)(const shape2 &shape, double dx)>(m_dylib.load_symbol("initialize"));
	if( initialize_func ) {
		initialize_func(m_shape,m_dx);
	}
	//
	// Get function pointers
	auto fluid_func = reinterpret_cast<double(*)(const vec2d &)>(m_dylib.load_symbol("fluid"));
	m_solid_func = reinterpret_cast<double(*)(const vec2d &)>(m_dylib.load_symbol("solid"));
	m_draw_func = reinterpret_cast<void(*)(graphics_engine &,double)>(m_dylib.load_symbol("draw"));
	m_moving_solid_func = reinterpret_cast<std::pair<double,vec2d>(*)(double time, const vec2d &p)>(m_dylib.load_symbol("moving_solid"));
	m_check_inject_func = reinterpret_cast<bool(*)(double, double, double, unsigned)>(m_dylib.load_symbol("check_inject"));
	m_inject_func = reinterpret_cast<bool(*)(const vec2d &, double, double, double, unsigned, double &, vec2d &)>(m_dylib.load_symbol("inject"));
	m_post_inject_func = reinterpret_cast<void(*)(double, double, double, unsigned, double&)>(m_dylib.load_symbol("post_inject"));
	auto velocity_func = reinterpret_cast<vec2d(*)(const vec2d &)>(m_dylib.load_symbol("velocity"));
	m_gravity_func = reinterpret_cast<vec2d(*)(double)>(m_dylib.load_symbol("gravity"));
	m_set_boundary_flux = reinterpret_cast<void(*)( double, Real [DIM2][2] )>(m_dylib.load_symbol("set_boundary_flux"));
	m_combined_solid_func = [&]( const vec2d &p ) {
		double value (1.0);
		if( m_solid_func ) value = std::min(value,m_solid_func(p));
		if( m_moving_solid_func ) value = std::min(value,m_moving_solid_func(m_timestepper->get_current_time(),p).first);
		return value;
	};
	//
	if( m_moving_solid_func ) {
		m_macoctreeproject.set_moving_solid([&]( const vec2d &p) {
			return m_moving_solid_func(m_timestepper->get_current_time(),p).first;
		});
	}
	//
	m_accumulated_CFL = 0.0;
	m_grid_0.clear();
	m_grid_1.clear();
	m_macoctreehelper.initialize();
	//
	// Get the narrowband size
	m_narrowband_depth = 3.0;
	if( m_param.use_FLIP ) {
		unsigned cells;
		assert(m_flip->const_send_message("narrowband",&cells));
		m_narrowband_depth = cells * m_dx;
	}
	//
	unsigned n = 1;
	while( (m_shape/n).min() >= m_param.min_resolution ) {
		shape2 shape = m_shape/n;
		m_grid_0.add_layer(shape,n*m_dx);
		m_grid_1.add_layer(shape,n*m_dx);
		n *= 2;
	}
	//
	if( m_set_boundary_flux ) {
		flux_boundary_condition2 boundary_cond;
		m_set_boundary_flux(0.0,boundary_cond.velocity);
		m_grid_0.set_flux_boundary_condition(boundary_cond);
		m_grid_1.set_flux_boundary_condition(boundary_cond);
	}
	//
	int refinement_count = m_param.use_sizing_func ? m_param.initial_refinement : 1;
	int count (0);
	while( refinement_count-- ) {
		//
		std::swap(m_grid,m_grid_prev);
		if( m_param.use_sizing_func ) {
			if( count ) {
				m_macoctreesizingfunc.activate_cells(*m_grid_prev,*m_grid,m_combined_solid_func,nullptr);
			} else {
				m_grid->activate_cells([&](char depth, const vec2d &p) {
					return depth > 3;
				});
			}
		} else {
			m_grid->activate_cells(fluid_func,m_combined_solid_func);
		}
		m_grid->balance_layers();
		m_grid->assign_indices();
		m_grid->assign_levelset(fluid_func,m_combined_solid_func);
		//
		if( m_param.use_sizing_func ) {
			m_macoctreesizingfunc.compute_sizing_function(*m_grid_prev,*m_grid,0.0,m_combined_solid_func,[&]( const vec2d &p ) {
				return m_moving_solid_func ? m_moving_solid_func(m_timestepper->get_current_time(),p).second : vec2d();
			});
		}
		count ++;
	}
	//
	if( velocity_func ) {
		m_grid->set_velocity([&]( const vec2d &p, char dim ) {
			return velocity_func(p)[dim];
		});
	}
	//
	if( m_param.use_FLIP ) {
		shared_array2<Real> highres_fluid(m_shape);
		shared_macarray2<Real> highres_velocity(m_shape);
		//
		if(m_grid->reconstruct_fluid(highres_fluid(),[&]( int i, int j, const Real &levelset_value ) {
			const double value = fluid_func(m_dx*vec2i(i,j).cell());
			return value < m_dx && value > -m_narrowband_depth;
		})) {
			//
			m_grid->reconstruct_velocity(highres_velocity(),[&]( int dim, int i, int j) {
				const double value = fluid_func(m_dx*vec2i(i,j).face(dim));
				return value < m_dx && value > -m_narrowband_depth;
			});
			//
			m_flip->resample(highres_fluid(),m_combined_solid_func,highres_velocity());
		}
	}
	//
	const double sqrt2 = sqrt(2.0);
	//
	m_solid.initialize(m_shape.nodal());
	m_solid.set_as_levelset(sqrt2*m_dx);
	if( m_solid_func ) {
		std::vector<vec2i> valid_positions[m_parallel.get_thread_num()];
		m_parallel.for_each(m_shape.nodal(),[&](int i, int j, int tid) {
			double value = m_solid_func(m_dx*vec2i(i,j).nodal());
			if( std::abs(value) < sqrt2*m_dx ) {
				valid_positions[tid].push_back(vec2i(i,j));
			}
		});
		for( const auto &e : valid_positions ) for( const auto &pi : e ) {
			double value = m_solid_func(m_dx*pi.nodal());
			m_solid.set(pi,value);
		}
		m_solid.flood_fill();
	}
	//
	// Segregate region
	if( m_param.volume_correction && m_param.regional_volume_correction ) {
		//
		m_region_count = m_macoctreesegregator.segregate(*m_grid,m_regions);
		m_y_list.resize(m_region_count);
		m_macoctreesegregator.compute(*m_grid,m_regions,m_volumes);
		console::dump("m_initial_volume = %.2e\n",
			std::accumulate(m_volumes.begin(),m_volumes.end(),0.0));
	}
	//
	// Assemble matrix for visualization
	m_macoctreeproject.assemble_matrix(*m_grid);
	m_macoctreehelper.reset_focus();
	//
	m_camera->set_bounding_box(vec2d().v,m_shape.box(m_dx).v);
}
//
void macoctreeliquid2::drag( double x, double y, double z, double u, double v, double w ) {
	//
	double scale (1e3);
	for( char depth=0; depth < m_grid->layers.size(); ++depth ) {
		auto &layer = *m_grid->layers[depth];
		for( char dim : DIMS2 ) {
			vec2i pi = layer.shape.find_face(vec2d(x,y)/layer.dx,dim);
			if( layer.active_faces[dim].active(pi)) {
				m_grid->velocity[layer.active_faces[dim](pi)] += scale * vec2d(u,v)[dim];
			}
		}
	}
	//
	m_flip->update([&](const vec2r &p, vec2r &velocity, Real &mass, bool bullet) {
		if( (vec2d(x,y)-p).len() < m_dx ) {
			velocity += scale * vec2r(u,v);
		}
	});
}
//
void macoctreeliquid2::idle() {
	//
	// Get the current step
	const unsigned step = m_timestepper->get_step_count()+1;
	//
	// Compute time step
	double max_u_per_unit (0.0);
	for( size_t n=0; n<m_grid->velocity.size(); ++n ) max_u_per_unit = std::max(max_u_per_unit,(double)std::abs(m_grid->velocity[n]));
	const double dt = m_timestepper->advance(max_u_per_unit,m_grid->get_finest_dx());
	const double time = m_timestepper->get_current_time();
	const double CFL = m_timestepper->get_current_CFL();
	//
	// Begin (check) injecting fluid
	begin_inject_external_fluid(dt,time,step);
	//
	if( m_param.use_sizing_func ) {
		m_macoctreesizingfunc.compute_sizing_function(*m_grid_prev,*m_grid,dt,m_combined_solid_func,[&]( const vec2d &p ) {
				return m_moving_solid_func ? m_moving_solid_func(m_timestepper->get_current_time(),p).second : vec2d();
		});
	}
	//
	// Swap grid and remesh
	m_accumulated_CFL += CFL;
	std::swap(m_grid,m_grid_prev);
	if( m_do_inject || m_accumulated_CFL >= m_param.maximal_CFL_accumulation ) {
		//
		m_accumulated_CFL = 0.0;
		if( m_param.use_sizing_func ) {
			if( m_do_inject ) {
				m_macoctreesizingfunc.activate_cells(*m_grid_prev,*m_grid,m_combined_solid_func,
					[&]( const vec2d &p ) {
						vec2d u; double value;
						m_inject_func(p,m_dx,dt,time,step,value,u);
						return value;
					});
			} else {
				m_macoctreesizingfunc.activate_cells(*m_grid_prev,*m_grid,m_combined_solid_func,nullptr);
			}
		} else {
			m_grid->activate_cells([&]( const vec2d &p ){
				//
				double inject_levelset (std::numeric_limits<double>::max());
				if( m_do_inject ) {
					vec2d u; double value;
					if( m_inject_func(p,m_dx,dt,time,step,value,u)) {
						inject_levelset = value;
					}
				}
				vec2d u (m_grid_prev->sample_velocity(p));
				return std::min(inject_levelset,m_grid_prev->sample_levelset(p-dt*u));
			},m_combined_solid_func);
		}
		//
		m_grid->balance_layers();
		m_grid->assign_indices();
		//
	} else {
		m_grid->copy(*m_grid_prev);
	}
	//
	// High-res fluid level set and velocity
	shared_array2<Real> highres_fluid(m_shape);
	shared_macarray2<Real> highres_velocity(m_shape);
	shared_macarray2<Real> save_velocity(m_shape);
	//
	// Define reconstruct test functions
	auto reconstruct_test_func_cell = [&]( int i, int j, const Real &levelset_value ) {
		return levelset_value < (m_param.erode_width+2)*m_dx && levelset_value > -m_narrowband_depth-(m_param.erode_width+2)*m_dx;
	};
	//
	auto reconstruct_test_func_face = [&]( int dim, int i, int j) {
		return highres_fluid->active(m_shape.clamp(i,j)) || highres_fluid->active(m_shape.clamp(i-(dim==0),j-(dim==1)));
	};
	//
	// Update fluid levelset
	if( m_param.use_FLIP && m_flip->get_particle_count()) {
		//
		if(m_grid_prev->reconstruct_fluid(highres_fluid(),reconstruct_test_func_cell)) {
			m_flip->update(m_combined_solid_func,highres_fluid(),time,false);
			highres_fluid->parallel_actives([&]( int i, int j, auto &it ) {
				const auto &highres_layer = m_grid_prev->layers[0];
				m_grid_prev->levelset[highres_layer->active_cells(i,j)] = it();
			});
		}
	}
	//
	// Set boundary flux
	if( m_set_boundary_flux ) {
		flux_boundary_condition2 boundary_cond;
		m_set_boundary_flux(time,boundary_cond.velocity);
		m_grid->set_flux_boundary_condition(boundary_cond);
	}
	//
	// Advect level set
	m_grid->assign_levelset([&]( const vec2d &p ) {
		vec2d u (m_grid_prev->sample_velocity(p));
		return m_grid_prev->sample_levelset(p-dt*u);
	},m_combined_solid_func);
	//
	// Advect region
	std::vector<uint_type> regions01; // Region after the advection
	if( m_param.volume_correction && m_param.regional_volume_correction ) {
		m_macoctreesegregator.backtrace(*m_grid_prev,*m_grid,dt,
			[&]( const vec2d &p ) {
				vec2d u (m_grid_prev->sample_velocity(p));
				return m_grid_prev->sample_levelset(p-dt*u);
			},m_regions,regions01);
		m_macoctreesegregator.extrapolate_jacobi(*m_grid,regions01);
		m_macoctreesegregator.prune(*m_grid,regions01);
	}
	//
	// Advect FLIP particles
	if( m_param.use_FLIP ) {
		m_flip->advect(
			m_combined_solid_func,
			[&](const vec2d &p){
				return vec2d(m_grid_prev->sample_velocity(p));
			},
			m_timestepper->get_current_time(),dt);
	}
	//
	// Advect velocity
	if( m_param.maccormack ) {
		//
		using Real2 = struct { Real v[2] = {0.0, 0.0}; };
		std::vector<Real2> min_max_values(m_grid->face_count);
		std::vector<Real> u0(m_grid->face_count), u1(m_grid->face_count), _u0(m_grid->face_count);
		std::vector<Real> u_x(m_grid->face_count), u_y(m_grid->face_count);
		std::vector<char> near_surface_flag (m_grid->face_count);
		//
		m_grid->iterate_active_faces([&]( const face_id2 &face_id, int tid ) {
			const vec2d &p = m_grid->get_face_position(face_id);
			uint_type index = face_id.index;
			u_x[index] = m_grid_prev->sample_velocity(p,0);
			u_y[index] = m_grid_prev->sample_velocity(p,1);
			u0[index] = m_grid_prev->sample_velocity(p,face_id.dim);
		});
		m_grid->iterate_active_faces([&]( const face_id2 &face_id, int tid ) {
			const vec2d &p = m_grid->get_face_position(face_id);
			uint_type index = face_id.index;
			vec2d u (u_x[index],u_y[index]);
			u1[index] = m_grid_prev->sample_face(p-dt*u,face_id.dim,m_grid_prev->velocity,min_max_values[face_id.index].v);
			const double dx = m_grid->get_face_dx(face_id);
			near_surface_flag[index] = m_grid_prev->sample_levelset(p-dt*u) > -dx || m_grid_prev->sample_levelset(p+dt*u) > -dx;
		});
		m_grid->iterate_active_faces([&]( const face_id2 &face_id, int tid ) {
			uint_type index = face_id.index;
			if( ! near_surface_flag[index] ) {
				const vec2d &p = m_grid->get_face_position(face_id);
				vec2d u (u_x[index],u_y[index]);
				_u0[index] = m_grid->sample_face(p+dt*u,face_id.dim,u1);
			}
		});
		//
		m_grid->iterate_active_faces([&]( const face_id2 &face_id, int tid ) {
			const uint_type index = face_id.index;
			const Real *v = min_max_values[index].v;
			if( near_surface_flag[index] ) {
				m_grid->velocity[index] = u1[index];
			} else {
				m_grid->velocity[index] = std::max(v[0],std::min(v[1],u1[index]+0.5f*(u0[index]-_u0[index])));
			}
		});
	} else {
		m_grid->set_velocity([&]( const vec2d &p, char dim ) {
			vec2d u (m_grid_prev->sample_velocity(p));
			return m_grid_prev->sample_velocity(p-dt*u,dim);
		});
	}
	//
	if( m_param.use_FLIP && m_flip->get_particle_count()) {
		//
		// Define interpolation functions
		auto interpolate_fluid = [&]( const vec2d &p) {
			return m_grid->sample_levelset(p);
		};
		auto interpolate_velocity = [&]( const vec2d &p) {
			return vec2d(m_grid->sample_velocity(p));
		};
		//
		// Mark bullet particles
		m_flip->mark_bullet(
			[&](const vec2d &p){ return interpolate_fluid(p); },
			[&](const vec2d &p){ return interpolate_velocity(p); },
			m_timestepper->get_current_time()
		);
		//
		// Correct positions
		m_flip->correct([&](const vec2d &p){ return interpolate_fluid(p); },highres_velocity());
		//
		// Splat momentum and mass of FLIP particles onto grids
		shared_macarray2<macflip2_interface::mass_momentum2> mass_and_momentum(m_shape);
		m_flip->splat(time,mass_and_momentum());
		//
		// Velocity overwrite
		mass_and_momentum->parallel_actives([&](char dim, int i, int j, auto &it, int tn ) {
			auto &highres_layer = m_grid->layers[0];
			if( highres_layer->active_faces[dim].active(i,j)) {
				const auto value = it();
				double grid_mass = std::max(0.0,1.0-value.mass);
				Real &grid_velocity = m_grid->velocity[highres_layer->active_faces[dim](i,j)];
				grid_velocity = (grid_mass*grid_velocity+value.momentum) / (grid_mass+value.mass);
			}
		});
		//
		// Save the current velocity
		m_grid->reconstruct_velocity(save_velocity(),reconstruct_test_func_face);
	}
	//
	// Add gravity force
	const vec2d gravity = m_gravity_func ? m_gravity_func(time) : m_param.gravity;
	m_grid->iterate_active_faces([&]( const face_id2 &face_id, int tid ) {
		m_grid->velocity[face_id.index] += dt*gravity[face_id.dim];
	});
	//
	if( m_param.volume_correction && m_param.regional_volume_correction ) {
		//
		std::vector<Real> target_volumes1;
		std::vector<Real> volumes1; // Volume after the advection
		m_macoctreesegregator.compute(*m_grid,regions01,volumes1);
		//
		std::vector<uint_type> regions11; // New region list
		m_region_count = m_macoctreesegregator.segregate(*m_grid,regions11);
		std::vector<Real> y_lists0(m_y_list);
		//
		// regions01 ... Region label array of state 0 after the advection
		// regions11 ... New region label array of state 1 after the advection
		// m_volumes ... List of volumes of state 0 before the advection
		// volumes1 .... List of volumes of state 0 after the advection
		m_y_list.resize(m_region_count);
		m_macoctreesegregator.compute_target_volume(*m_grid,regions01,m_volumes,regions11,volumes1,y_lists0,m_current_volumes,target_volumes1,m_y_list);
		m_regions = regions11;
		m_volumes = target_volumes1;
		//
		m_macoctreesegregator.compute(*m_grid,m_regions,volumes1);
	}
	//
	// Inject external fluid
	do_inject_external_fluid(dt,time,step);
	//
	// Ending injecting fluid
	end_inject_external_fluid(dt,time,step);
	//
	// Add surface tension force
	if( m_param.surftens_k ) m_grid->add_surfacetense_force(m_param.surftens_k,dt);
	//
	// Assemble matrix
	m_macoctreeproject.assemble_matrix(*m_grid);
	m_macoctreehelper.reset_focus();
	//
	// Define solid velocity
	auto solid_velocity_func = [&]( const vec2d &p ) {
		return m_moving_solid_func ? m_moving_solid_func(time,p).second : vec2d();
	};
	//
	// Project
	if( m_param.volume_correction && m_param.regional_volume_correction ) {
		m_macoctreeproject.project(*m_grid,dt,m_region_count,m_regions,m_current_volumes,m_volumes,m_y_list,solid_velocity_func);
	} else {
		m_macoctreeproject.project(*m_grid,dt,solid_velocity_func);
	}
	m_grid->extrapolate_toward_solid(m_combined_solid_func);
	m_grid->extrapolate(m_combined_solid_func);
	//
	if( m_param.use_FLIP ) {
		//
		// Copy new high-res velocity field and the level set
		if( m_grid->reconstruct_fluid(highres_fluid(),reconstruct_test_func_cell)) {
			//
			m_grid->reconstruct_velocity(highres_velocity(),reconstruct_test_func_face);
			if( m_param.erode_width ) highres_fluid->erode(m_param.erode_width);
			//
			// Reseed particles
			m_flip->resample(highres_fluid(),m_combined_solid_func,highres_velocity());
			//
			// Remove particles
			m_flip->remove([&](const vec2r &p, bool bullet) {
				if( bullet && m_grid->sample_levelset(p) < 0.0 ) {
					if( highres_fluid->active(m_shape.find_cell(p/m_dx))) return false;
					else return true;
				} else {
					if( ! bullet && ! highres_fluid->active(m_shape.find_cell(p/m_dx))) return true;
					else return false;
				}
				return false;
			});
			//
			// Update FLIP velocity
			m_flip->update(save_velocity(),highres_velocity(),dt,gravity,m_param.PICFLIP);
			//
		} else {
			m_flip->remove([](const vec2r &p, bool bullet){ return true; });
		}
	}
}
//
void macoctreeliquid2::begin_inject_external_fluid( double dt, double time, unsigned step ) {
	m_do_inject = m_check_inject_func && m_inject_func && m_check_inject_func(m_dx,dt,time,step);
}
//
void macoctreeliquid2::do_inject_external_fluid( double dt, double time, unsigned step ) {
	//
	if( m_do_inject ) {
		//
		std::vector<double> total_injected (m_grid->parallel.get_thread_num());
		m_grid->iterate_active_cells([&]( const cell_id2 &cell_id, int tid ) {
			//
			const vec2d p = m_grid->get_cell_position(cell_id);
			double value (0.0); vec2d u;
			//
			if( m_inject_func(p,m_dx,dt,time,step,value,u)) {
				const double v = m_grid->levelset[cell_id.index];
				m_grid->levelset[cell_id.index] = std::min(value,v);
				if( value < 0.0 && v > 0.0 ) {
					const double dx = m_grid->get_cell_dx(cell_id);
					total_injected[tid] += dx*dx;
				}
			}
		});
		m_grid->iterate_active_faces([&]( const face_id2 &face_id, int tid ) {
			const vec2d p = m_grid->get_face_position(face_id);
			const double dx = m_grid->get_face_dx(face_id);
			double value (0.0); vec2d u;
			if( m_inject_func(p,m_dx,dt,time,step,value,u)) {
				if( value < dx ) {
					m_grid->velocity[face_id.index] = u[face_id.dim];
				}
			}
		});
		m_injected_volume = std::accumulate(total_injected.begin(),total_injected.end(),0.0);
	}
}
//
void macoctreeliquid2::end_inject_external_fluid( double dt, double time, unsigned step ) {
	//
	auto extrapolate_with_volume_change = [&]() {
		//
		std::vector<std::vector<double> > volume_changes(m_region_count);
		for( auto &e : volume_changes ) e.resize(m_grid->parallel.get_thread_num());
		//
		while(true) {
			//
			std::vector<unsigned char> found_bucket(m_grid->parallel.get_thread_num());
			std::vector<uint_type> regions_save(m_regions);
			m_grid->iterate_active_cells([&]( const cell_id2 &cell_id, int tid ) {
				bool found (false);
				if( m_grid->levelset[cell_id.index] < 0.0 && ! regions_save[cell_id.index]) {
					uint_type new_index (0);
					m_grid->iterate_cell_neighbors(cell_id,[&]( char dim, const cell_id2& neighbor_cell_id ){
						if( m_grid->levelset[neighbor_cell_id.index] < 0.0 && regions_save[neighbor_cell_id.index] ) {
							new_index = regions_save[neighbor_cell_id.index];
						}
					});
					if( new_index ) {
						const double v = m_grid->get_cell_volume(cell_id);
						m_regions[cell_id.index] = new_index;
						volume_changes[new_index-1][tid] += v;
						found = true;
					}
				}
				found_bucket[tid] = found_bucket[tid] || found;
			});
			if( ! std::accumulate(found_bucket.begin(),found_bucket.end(),0) ) break;
		}
		//
		for( unsigned n=0; n<m_region_count; ++n ) {
			const double v = std::accumulate(volume_changes[n].begin(),volume_changes[n].end(),0.0);
			m_volumes[n] += v;
			m_current_volumes[n] += v;
		}
		//
		bool region_added (false);
		m_grid->serial_iterate_active_cells([&]( const cell_id2 &cell_id ) {
			if( m_grid->levelset[cell_id.index] < 0.0 && ! m_regions[cell_id.index]) {
				const double v = m_grid->get_cell_volume(cell_id);
				if( ! region_added ) {
					m_y_list.push_back(0.0);
					m_current_volumes.push_back(v);
					m_volumes.push_back(v);
					m_region_count ++;
					region_added = true;
				} else {
					m_current_volumes[m_region_count-1] += v;
					m_volumes[m_region_count-1] += v;
				}
				m_regions[cell_id.index] = m_region_count;
			}
		});
	};
	//
	if( m_do_inject ) {
		if( m_post_inject_func ) {
			m_post_inject_func(m_dx,dt,time,step,m_injected_volume);
		}
		if( m_param.volume_correction && m_injected_volume ) {
			if( m_param.regional_volume_correction ) {
				extrapolate_with_volume_change();
			} else {
				if( ! m_macoctreeproject.m_initial_volume ) {
					m_macoctreeproject.m_initial_volume = m_grid->get_volume();
				} else {
					m_macoctreeproject.m_initial_volume += m_injected_volume;
				}
			}
		}
	}
}
//
void macoctreeliquid2::setup_window( std::string &name, int &width, int &height ) const {
	height = width * (m_shape[1]/ (double)m_shape[0]);
}
//
void macoctreeliquid2::cursor( double x, double y, double z ) {
	m_macoctreehelper.cursor(*m_grid,x,y,z);
}
//
bool macoctreeliquid2::keyboard( int key, int action, int mods ) {
	return m_macoctreehelper.keyboard(key,action,mods);
}
//
void macoctreeliquid2::draw( graphics_engine &g ) const {
	//
	const double time = m_timestepper->get_current_time();
	//
	g.color4(0.9,0.6,0.3,0.5);
	m_solid_gridvisualizer->draw_levelset(g,m_solid);
	//
	// Draw from scene library
	if( m_draw_func ) m_draw_func(g,time);
	//
	// Draw octree grid
	m_grid->draw_fluid(g);
	m_grid->draw_grid(g);
	//
	// Draw FLIP
	m_flip->draw(g,time);
	//
	// Draw debug info
	m_macoctreehelper.draw_debug(*m_grid,m_macoctreeproject.m_matrix,g);
	//
	if( m_param.regional_volume_correction ) {
		//
		// Draw region
		m_macoctreesegregator.draw_region(g,*m_grid,m_regions);
	}
}
//
extern "C" module * create_instance() {
	return new macoctreeliquid2;
}
//

// --- 適応的グリッド（AMR: Adaptive Mesh Refinement）の実装ポイント ---
//
// 1. グリッド構造
//   - m_grid, m_grid_prev, m_grid_0, m_grid_1 などが octree ベースの多層グリッドを保持。
//   - m_grid->layers で各解像度レベルのグリッドを管理（octreeの各深さが異なる解像度）。
//
// 2. 初期化・リファインメント
//   - post_initialize() 内で
//     - m_shape/n の最小値が min_resolution 以上になるまで add_layer() で多層グリッドを構築。
//     - m_param.initial_refinement 回数だけ activate_cells() で細分化（リファインメント）を実施。
//     - use_sizing_func=true の場合は m_macoctreesizingfunc.activate_cells() でセルごとに細分化基準を適用。
//     - それ以外は fluid_func, m_combined_solid_func で流体/固体領域に応じて細分化。
//     - balance_layers() で隣接セルのレベル差を1以内に調整（octreeのバランス化）。
//
// 3. 時間発展中のリメッシュ
//   - idle() 内で
//     - m_accumulated_CFL が閾値を超える or 流体注入時に std::swap(m_grid, m_grid_prev) でグリッドを入れ替え、
//       activate_cells() で再細分化（リファインメント）を実施。
//     - use_sizing_func=true の場合は sizing function に基づきセルごとに細分化。
//     - それ以外は流体/固体領域や流体注入に応じて細分化。
//     - balance_layers(), assign_indices() でグリッドを再構成。
//
// 4. 細分化基準
//   - use_sizing_func=true の場合: m_macoctreesizingfunc.compute_sizing_function() で物理量やユーザー関数に基づき細分化レベルを決定。
//   - それ以外: fluid_func, m_combined_solid_func で流体/固体の存在や注入条件に応じてセルを細分化。
//   - activate_cells() のラムダでセルごとに細分化/縮退の判定を行う。
//
// 5. その他
//   - balance_layers() でoctreeのバランスを保ち、隣接セルのレベル差が大きくならないようにしている。
//   - assign_levelset(), assign_indices() で各セルの物理量やインデックスを再割り当て。
//   - m_grid->reconstruct_fluid(), m_grid->reconstruct_velocity() で高解像度フィールドを再構築。

// --- まとめ ---
// ・octreeベースの多層グリッドを使い、セルごとに細分化/縮退をactivate_cells()で判定。
// ・細分化基準はsizing functionまたは流体/固体領域・注入条件など。
// ・balance_layers()でoctreeのバランスを維持。
// ・時間発展ごとに必要に応じてグリッドをリファイン・リメッシュしている。