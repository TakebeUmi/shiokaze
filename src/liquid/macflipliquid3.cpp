/*
**	macflipliquid3.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on June 2, 2017.
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
#include "macflipliquid3.h"
#include <shiokaze/array/shared_array3.h>
#include <shiokaze/array/shared_bitarray3.h>
#include <shiokaze/array/macarray_extrapolator3.h>
#include <shiokaze/array/macarray_interpolator3.h>
#include <shiokaze/array/array_upsampler3.h>
#include <shiokaze/array/array_interpolator3.h>
#include <shiokaze/array/array_utility3.h>
#include <shiokaze/graphics/graphics_utility.h>
#include <shiokaze/core/console.h>
#include <shiokaze/core/timer.h>
#include <shiokaze/core/filesystem.h>
#include <cmath>
//
SHKZ_USING_NAMESPACE
//
macflipliquid3::macflipliquid3 () {
	m_param.PICFLIP = 0.95;
}
//
void macflipliquid3::configure( configuration &config ) {
	//
	config.get_double("PICFLIP",m_param.PICFLIP,"PICFLIP blending factor");
	config.get_bool("DisableResample",m_param.disable_resample,"Disable resampling");
	config.get_bool("PreupdateFLIP",m_param.preupdate_FLIP,"Pre-update FLIP");
	assert( m_param.PICFLIP >= 0.0 && m_param.PICFLIP <= 1.0 );
	//
	macliquid3::configure(config);
}
//
void macflipliquid3::post_initialize ( bool initialized_from_file ) {
	//
	macliquid3::post_initialize(initialized_from_file);
	//
	if( ! initialized_from_file ) {
		scoped_timer timer(this);
		timer.tick(); console::dump( ">>> Started FLIP initialization\n" );
		//
		extend_both();
		m_flip->resample(m_fluid,[&](const vec3d &p){ return interpolate_solid(p); },m_velocity);
		//
		console::dump( "<<< Initialization finished. Took %s\n", timer.stock("initialization").c_str());
	}
}
//
size_t macflipliquid3::do_inject_external_fluid( array3<Real> &fluid, macarray3<Real> &velocity, double dt, double time, unsigned step ) {
	//
	const size_t count = macliquid3::do_inject_external_fluid(fluid,velocity,dt,time,step);
	m_flip->resample(m_fluid,[&](const vec3d &p){ return interpolate_solid(p); },m_velocity,[&]( const vec3d &p ) {
		double value (0.0); vec3d u;
		if( this->m_inject_func(p,m_dx,dt,time,step,value,u)) {
			return value < 0.0;
		} else {
			return false;
		}
	});
	return count;
}
//
void macflipliquid3::idle() {
	//
	scoped_timer timer(this);
	//
	// Add to graph
	add_to_graph();
	//
	// Compute the timestep size
	const double dt = m_timestepper->advance(m_macutility->compute_max_u(m_velocity),m_dx);
	const double time = m_timestepper->get_current_time();
	const double CFL = m_timestepper->get_current_CFL();
	const unsigned step = m_timestepper->get_step_count();
	timer.tick(); console::dump( ">>> %s step started (dt=%.2e,CFL=%.2f)...\n", dt, CFL, console::nth(step).c_str());
	//
	shared_macarray3<Real> save_velocity(m_shape);
	shared_macarray3<macflip3_interface::mass_momentum3> mass_and_momentum(m_shape);
	//
	// Update solid
	m_macutility->update_solid_variables(m_dylib,time,&m_solid,&m_solid_velocity);
	//
	// Update fluid levelset
	if( m_param.preupdate_FLIP ) {
		m_flip->update([&](const vec3d &p){ return interpolate_solid(p); },m_fluid,time);
	}
	//
	// Advect fluid levelset
	shared_array3<Real> fluid_save(m_fluid);
	m_macadvection->advect_scalar(m_fluid,m_velocity,fluid_save(),dt,"levelset");
	m_fluid.flood_fill();
	//
	// Advect FLIP particles
	m_flip->advect(
		[&](const vec3d &p){ return interpolate_solid(p); },
		[&](const vec3d &p){ return interpolate_velocity(p); },
		time,dt);
	//
	// Update fluid levelset
	if( ! m_param.preupdate_FLIP ) {
		m_flip->update([&](const vec3d &p){ return interpolate_solid(p); },m_fluid,time);
	}
	//
	// Redistance
	timer.tick(); console::dump( "Re-distancing fluid levelsets..." );
	m_redistancer->redistance(m_fluid,m_param.levelset_half_bandwidth_count);
	m_gridutility->extrapolate_levelset(m_solid,m_fluid);
	console::dump( "Done. Took %s\n", timer.stock("redistance_levelset").c_str());
	//
	// Grid velocity advection
	m_macadvection->advect_vector(m_velocity,m_velocity,m_fluid,dt,"velocity");
	//
	// Splat momentum and mass of FLIP particles onto grids
	m_flip->splat(time,mass_and_momentum());
	//
	// Compute the combined grid velocity
	timer.tick(); console::dump( "Computing combined grid velocity..." );
	//
	shared_macarray3<Real> overwritten_velocity(m_shape);
	overwritten_velocity->activate_as(mass_and_momentum());
	overwritten_velocity->parallel_actives([&](int dim, int i, int j, int k, auto &it, int tn ) {
		const auto value = mass_and_momentum()[dim](i,j,k);
		const Real grid_mass = m_velocity[dim].active(i,j,k) ? std::max(0.0,1.0-value.mass) : 0.0;
		it.set((grid_mass*m_velocity[dim](i,j,k)+value.momentum)/(grid_mass+value.mass));
	});
	//
	// Velocity overwrite
	overwritten_velocity->const_serial_actives([&](int dim, int i, int j, int k, auto &it) {
		m_velocity[dim].set(i,j,k,it());
	});
	//
	// Mark bullet particles
	m_flip->mark_bullet(
		[&](const vec3d &p){ return interpolate_fluid(p); },
		[&](const vec3d &p){ return interpolate_velocity(p); },
		time
	);
	//
	console::dump( "Done. Took %s\n", timer.stock("compute_combined_velocity").c_str());
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
	m_macproject->project(dt,m_velocity,m_solid,m_solid_velocity,m_fluid,(macliquid3::m_param).surftens_k);
	//
	// Extend both the level set and velocity
	extend_both();
	//
	// Correct positions
	m_flip->correct([&](const vec3d &p){ return interpolate_fluid(p); },m_velocity);
	//
	// Reseed particles
	if( ! m_param.disable_resample ) {
		m_flip->resample(m_fluid,
			[&](const vec3d &p){ return interpolate_solid(p); },
			m_velocity
		);
	}
	//
	// Update FLIP momentum
	m_flip->update(save_velocity(),m_velocity,dt,(macliquid3::m_param).gravity,m_param.PICFLIP);
	//
	console::dump( "<<< %s step done. Took %s\n", console::nth(step).c_str(), timer.stock("simstep").c_str());
	//
	// Export mesh
	export_mesh();
	//
	// Save status
	save_state();
	//
	// Report stats
	m_macstats->dump_stats(m_solid,m_fluid,m_velocity,m_timestepper.get());
}
//
void macflipliquid3::do_export_mesh( unsigned frame ) const {
	//
	scoped_timer timer(this);
	std::string particle_path = console::format_str("%s/%d_particles.dat",m_export_path.c_str(),frame);
	timer.tick(); console::dump( "Writing ballistic particles..." );
	//
	std::vector<particlerasterizer3_interface::Particle3> ballistic_points;
	std::vector<macflip3_interface::particle3> particles = m_flip->get_particles();
	for( int n=0; n<particles.size(); ++n ) {
		//
		particlerasterizer3_interface::Particle3 point;
		point.p = particles[n].p;
		point.r = particles[n].r;
		//
		if( interpolate_fluid(particles[n].p) > point.r ) {
			ballistic_points.push_back(point);
		}
	}
	//
	FILE *fp = fopen(particle_path.c_str(),"wb");
	const size_t size = ballistic_points.size();
	fwrite(&size,1,sizeof(unsigned),fp);
	for( size_t n=0; n<size; ++n ) {
		float position[3] = { (float)ballistic_points[n].p.v[0],
							  (float)ballistic_points[n].p.v[1],
							  (float)ballistic_points[n].p.v[2] };
		float radius = ballistic_points[n].r;
		fwrite(position,3,sizeof(float),fp);
		fwrite(&radius,1,sizeof(float),fp);
	}
	fclose(fp);
	console::dump( "Done. Size=%d. Took %s\n", size, timer.stock("write_ballistic").c_str());
	//
	macliquid3::do_export_mesh(frame);
}
//
void macflipliquid3::render_mesh( unsigned frame ) const {
	//
	scoped_timer timer(this);
	global_timer::pause();
	//
	assert(console::get_root_path().size());
	//
	std::string mitsuba_path = console::get_root_path() + "/flipliquid_mitsuba";
	std::string copy_from_path = filesystem::find_resource_path("flipliquid","mitsuba");
	if( ! filesystem::is_exist(mitsuba_path)) {
		if( filesystem::is_exist(copy_from_path)) {
			console::run( "cp -r %s %s", copy_from_path.c_str(), mitsuba_path.c_str());
		} else {
			console::dump( "Could not lcoate mitsuba files (%s).\n", copy_from_path.c_str());
			exit(0);
		}
	}
	//
	// Write variables
	FILE *fp = fopen((m_export_path+"/common.ini").c_str(),"w");
	assert(fp);
	fprintf(fp,"[Common]\n");
	fprintf(fp,"OriginPos = [%f,%f,%f]\n",(macliquid3::m_param).origin[0],(macliquid3::m_param).origin[1],(macliquid3::m_param).origin[2]);
	fprintf(fp,"TargetPos = [%f,%f,%f]\n",(macliquid3::m_param).target[0],(macliquid3::m_param).target[1],(macliquid3::m_param).target[2]);
	fprintf(fp,"SampleCount = %d\n", (macliquid3::m_param).render_sample_count );
	fprintf(fp,"LiquidColor = [%f,%f,%f]\n", 0.5, 0.5, 1.0);
	fclose(fp);
	//
	std::string render_command = console::format_str("cd %s; /usr/bin/python render.py %d mesh",
				mitsuba_path.c_str(),
				frame);
	//
	console::dump("Running command: %s\n", render_command.c_str());
	console::system(render_command.c_str());
	//
	if( (macliquid3::m_param).render_transparent ) {
		//
		std::string render_command = console::format_str("cd %s; /usr/bin/python render.py %d transparent",
				mitsuba_path.c_str(),
				frame);
		//
		console::dump("Running command: %s\n", render_command.c_str());
		console::system(render_command.c_str());
	}
	//
	global_timer::resume();
}
//
double macflipliquid3::interpolate_fluid( const vec3d &p ) const {
	return array_interpolator3::interpolate(m_fluid,p/m_dx-vec3d(0.5,0.5,0.5));
}
//
double macflipliquid3::interpolate_solid( const vec3d &p ) const {
	return array_interpolator3::interpolate(m_solid,p/m_dx);
}
//
vec3d macflipliquid3::interpolate_velocity( const vec3d &p ) const {
	return macarray_interpolator3::interpolate(m_velocity,vec3d(),m_dx,p);
}
//
void macflipliquid3::draw( graphics_engine &g ) const {
	//
	macliquid3::draw(g);
	//
	// Draw FLIP particles
	m_flip->draw(g,m_timestepper->get_current_time());
}
//
extern "C" module * create_instance() {
	return new macflipliquid3;
}
//
extern "C" const char *license() {
	return "MIT";
}
//