/*
**	macflipliquid2.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on April 17, 2017.
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
#include "macflipliquid2.h"
#include <shiokaze/array/shared_array2.h>
#include <shiokaze/array/shared_bitarray2.h>
#include <shiokaze/array/macarray_extrapolator2.h>
#include <shiokaze/array/macarray_interpolator2.h>
#include <shiokaze/array/array_upsampler2.h>
#include <shiokaze/array/array_interpolator2.h>
#include <shiokaze/utility/utility.h>
#include <cmath>
//
SHKZ_USING_NAMESPACE
//
macflipliquid2::macflipliquid2 () {
	m_param.PICFLIP = 0.95;
}
//
void macflipliquid2::configure( configuration &config ) {
	//
	config.get_double("PICFLIP",m_param.PICFLIP,"PICFLIP blending factor");
	config.get_bool("DisableResample",m_param.disable_resample,"Disable resampling");
	config.get_unsigned("LevelsetHalfWidth",m_param.levelset_half_bandwidth_count,"Level set half bandwidth");
	config.get_bool("PreupdateFLIP",m_param.preupdate_FLIP,"Pre-update FLIP");
	assert( m_param.PICFLIP >= 0.0 && m_param.PICFLIP <= 1.0 );
	//
	macliquid2::configure(config);
}
//
void macflipliquid2::post_initialize ( bool initialized_from_file ) {
	//
	macliquid2::post_initialize(initialized_from_file);
	if( ! initialized_from_file ) {
		extend_both();
		m_flip->resample(m_fluid,[&](const vec2d &p){ return interpolate_solid(p); },m_velocity);
	}
}
//
size_t macflipliquid2::do_inject_external_fluid( array2<Real> &fluid, macarray2<Real> &velocity, double dt, double time, unsigned step ) {
	//
	size_t count = macliquid2::do_inject_external_fluid(fluid,velocity,dt,time,step);
	m_flip->resample(m_fluid,[&](const vec2d &p){ return interpolate_solid(p); },m_velocity,[&]( const vec2d &p ) {
		double value (0.0); vec2d u;
		if( this->m_inject_func(p,m_dx,dt,time,step,value,u)) {
			return value < 0.0;
		} else {
			return false;
		}
	});
	return count;
}
//
void macflipliquid2::idle() {
	//
	// Add to graph
	add_to_graph();
	//
	// Compute the timestep size
	const double dt = m_timestepper->advance(m_macutility->compute_max_u(m_velocity),m_dx);
	const double time = m_timestepper->get_current_time();
	const unsigned step = m_timestepper->get_step_count();
	//
	shared_macarray2<Real> save_velocity(m_shape);
	shared_macarray2<macflip2_interface::mass_momentum2> mass_and_momentum(m_shape);
	//
	// Update solid
	m_macutility->update_solid_variables(m_dylib,time,&m_solid,&m_solid_velocity);
	//
	// Update fluid levelset (pre)
	if( m_param.preupdate_FLIP ) {
		m_flip->update([&](const vec2d &p){ return interpolate_solid(p); },m_fluid,time);
	}
	//
	// Advect fluid levelset
	shared_array2<Real> fluid_save(m_fluid);
	m_macadvection->advect_scalar(m_fluid,m_velocity,fluid_save(),dt);
	m_fluid.flood_fill();
	//
	// Advect FLIP particles
	m_flip->advect(
		[&](const vec2d &p){ return interpolate_solid(p); },
		[&](const vec2d &p){ return interpolate_velocity(p); },
		time,dt);
	//
	// Update fluid levelset (post)
	if( ! m_param.preupdate_FLIP ) {
		m_flip->update([&](const vec2d &p){ return interpolate_solid(p); },m_fluid,time);
	}
	//
	// Redistance
	m_redistancer->redistance(m_fluid,m_param.levelset_half_bandwidth_count);
	m_gridutility->extrapolate_levelset(m_solid,m_fluid);
	//
	// Grid velocity advection
	m_macadvection->advect_vector(m_velocity,m_velocity,m_fluid,dt);
	//
	// Splat momentum and mass of FLIP particles onto grids
	m_flip->splat(time,mass_and_momentum());
	//
	// Compute the combined grid velocity
	shared_macarray2<Real> overwritten_velocity(m_shape);
	overwritten_velocity->activate_as(mass_and_momentum());
	overwritten_velocity->parallel_actives([&](int dim, int i, int j, auto &it, int tn ) {
		const auto value = mass_and_momentum()[dim](i,j);
		const Real grid_mass = m_velocity[dim].active(i,j) ? std::max(0.0,1.0-value.mass) : 0.0;
		it.set((grid_mass*m_velocity[dim](i,j)+value.momentum)/(grid_mass+value.mass));
	});
	//
	// Velocity overwrite
	overwritten_velocity->const_serial_actives([&](int dim, int i, int j, auto &it) {
		m_velocity[dim].set(i,j,it());
	});
	//
	// Mark bullet particles
	m_flip->mark_bullet(
		[&](const vec2d &p){ return interpolate_fluid(p); },
		[&](const vec2d &p){ return interpolate_velocity(p); },
		time
	);
	//
	// Save the current velocity
	save_velocity->copy(m_velocity);
	//
	// Add external force
	inject_external_force(m_velocity,dt);
	//
	// Inject external fluid
	inject_external_fluid(m_fluid,m_velocity,dt,time);
	//
	// Set volume correction
	set_volume_correction(m_macproject.get());
	//
	// Project
	m_macproject->project(dt,m_velocity,m_solid,m_solid_velocity,m_fluid,(macliquid2::m_param).surftens_k);
	//
	// Extend both the level set and velocity
	extend_both();
	//
	// Correct positions
	m_flip->correct([&](const vec2d &p){ return interpolate_fluid(p); },m_velocity);
	//
	// Reseed particles
	if( ! m_param.disable_resample ) {
		m_flip->resample(m_fluid,[&](const vec2d &p){ return interpolate_solid(p); },m_velocity);
	}
	//
	// Update FLIP velocity
	m_flip->update(save_velocity(),m_velocity,dt,(macliquid2::m_param).gravity,m_param.PICFLIP);
	//
	// Report stats
	m_macstats->dump_stats(m_solid,m_fluid,m_velocity,m_timestepper.get());
}
//
double macflipliquid2::interpolate_fluid( const vec2d &p ) const {
	return array_interpolator2::interpolate(m_fluid,p/m_dx-vec2d(0.5,0.5));
}
//
double macflipliquid2::interpolate_solid( const vec2d &p ) const {
	return array_interpolator2::interpolate(m_solid,p/m_dx);
}
//
vec2d macflipliquid2::interpolate_velocity( const vec2d &p ) const {
	return macarray_interpolator2::interpolate(m_velocity,vec2d(),m_dx,p);
}
//
void macflipliquid2::draw( graphics_engine &g ) const {
	//
	macliquid2::draw(g);
	//
	// Draw FLIP
	m_flip->draw(g,m_timestepper->get_current_time());
}
//
extern "C" module * create_instance() {
	return new macflipliquid2;
}
//
extern "C" const char *license() {
	return "MIT";
}
//