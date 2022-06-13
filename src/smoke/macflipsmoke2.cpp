/*
**	macflipsmoke2.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on Aug 22, 2017.
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
#include "macflipsmoke2.h"
#include <shiokaze/array/shared_array2.h>
#include <shiokaze/array/macarray_extrapolator2.h>
#include <shiokaze/array/macarray_interpolator2.h>
#include <shiokaze/array/array_interpolator2.h>
#include <shiokaze/graphics/graphics_utility.h>
#include <cmath>
//
SHKZ_USING_NAMESPACE
//
macflipsmoke2::macflipsmoke2 () {
	//
	m_param.PICFLIP = 0.95;
	m_param.gridmass = 1.0;
}
//
void macflipsmoke2::configure( configuration &config ) {
	//
	config.get_double("GridMass",m_param.gridmass,"Mass of grid cell");
	config.get_double("PICFLIP",m_param.PICFLIP,"PICFLIP blending factor");
	assert( m_param.PICFLIP >= 0.0 && m_param.PICFLIP <= 1.0 );
	//
	macsmoke2::configure(config);
}
//
void macflipsmoke2::post_initialize ( bool initialized_from_file ) {
	//
	macsmoke2::post_initialize(initialized_from_file);
	//
	shared_array2<Real> fluid(m_shape);
	fluid->set_as_levelset(m_dx);
	m_flip->resample(fluid(),[&](const vec2d &p){ return interpolate_solid(p); },m_velocity);
}
//
void macflipsmoke2::idle() {
	//
	// Add to graph
	add_to_graph();
	//
	// Compute the timestep size
	const double dt = m_timestepper->advance(m_macutility->compute_max_u(m_velocity),m_dx);
	const double time = m_timestepper->get_current_time();
	//
	// Update solid
	m_macutility->update_solid_variables(m_dylib,time,&m_solid,&m_solid_velocity);
	//
	// Advect FLIP particles and get the levelset after the advection
	m_flip->advect(
		[&](const vec2d &p){ return interpolate_solid(p); },
		[&](const vec2d &p){ return interpolate_velocity(p); },
		time,dt);
	//
	// Correct positions
	m_flip->correct([&](const vec2d &p){ return -1.0; },m_velocity);
	//
	// Reseed particles
	shared_array2<Real> fluid(m_shape);
	fluid->set_as_levelset(m_dx);
	m_flip->resample(m_fluid,
		[&](const vec2d &p){ return interpolate_solid(p); },
		m_velocity
	);
	//
	// Advect density and velocity
	if( (macsmoke2::m_param).use_dust ) advect_dust_particles(m_velocity,dt);
	else {
		m_density.dilate(std::ceil(m_timestepper->get_current_CFL()));
		m_macadvection->advect_scalar(m_density,m_velocity,m_fluid,dt);
		const double minimal_density = (macsmoke2::m_param).minimal_density;
		m_density.parallel_actives([&](auto &it) {
			if( std::abs(it()) <= minimal_density ) it.set_off();
		});
	}
	//
	// Splat momentum and mass of FLIP particles onto grids
	shared_macarray2<macflip2_interface::mass_momentum2> mass_and_momentum(m_shape);
	m_flip->splat(time,mass_and_momentum());
	//
	// Overwrite grid velocity
	m_velocity.parallel_actives([&]( int dim, int i, int j, auto &it, int tn) {
		const auto value = mass_and_momentum()[dim](i,j);
		if( value.mass ) it.set(value.momentum / value.mass);
	});
	//
	// Save the current velocity
	shared_macarray2<Real> save_velocity(m_velocity);
	//
	// Add external force
	inject_external_force(m_velocity);
	//
	// Add buoyancy force
	add_buoyancy_force(m_velocity,m_density,dt);
	//
	// Add source
	add_source(m_velocity,m_density,m_timestepper->get_current_time(),dt);
	//
	// Project
	m_macproject->project(dt,m_velocity,m_solid,m_solid_velocity,m_fluid,0.0);
	//
	// Extrapolate
	m_macutility->extrapolate_and_constrain_velocity(m_solid,m_velocity,(macsmoke2::m_param).extrapolated_width);
	//
	// Update FLIP momentum
	m_flip->update(save_velocity(),m_velocity,dt,vec2d(),m_param.PICFLIP);
	//
	// Report stats
	m_macstats->dump_stats(m_solid,m_fluid,m_velocity,m_timestepper.get());
}
//
double macflipsmoke2::interpolate_solid( const vec2d &p ) const {
	return array_interpolator2::interpolate(m_solid,p/m_dx);
}
//
vec2d macflipsmoke2::interpolate_velocity( const vec2d &p ) const {
	return macarray_interpolator2::interpolate(m_velocity,vec2d(),m_dx,p);
}
//
void macflipsmoke2::draw( graphics_engine &g ) const {
	//
	macsmoke2::draw(g);
	//
	// Draw FLIP
	m_flip->draw(g,m_timestepper->get_current_time());
}
//
extern "C" module * create_instance() {
	return new macflipsmoke2;
}
//
extern "C" const char *license() {
	return "MIT";
}
//