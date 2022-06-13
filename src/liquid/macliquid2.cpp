/*
**	macliquid2.cpp
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
#include "macliquid2.h"
#include <shiokaze/core/filesystem.h>
#include <shiokaze/array/macarray_extrapolator2.h>
#include <shiokaze/array/macarray_interpolator2.h>
#include <shiokaze/array/shared_bitarray2.h>
#include <shiokaze/utility/utility.h>
#include <shiokaze/core/dylibloader.h>
#include <cmath>
#include <numeric>
//
SHKZ_USING_NAMESPACE
//
macliquid2::macliquid2 () {
	//
	m_param.gravity = vec2d(0.0,-9.8);
	m_param.surftens_k = 0.0;
	m_param.volume_correction = true;
	m_param.volume_change_tol_ratio = 0.03;
	m_param.show_graph = false;
	//
	m_shape = shape2{64,32};
}
//
void macliquid2::load( configuration &config ) {
	//
	std::string name("waterdrop2"); config.get_string("Name",name,"Scene file name");
	m_dylib.open_library(filesystem::resolve_libname(name));
	m_dylib.load(config);
	m_dylib.overwrite(config);
}
//
void macliquid2::configure( configuration &config ) {
	//
	m_dylib.configure(config);
	//
	config.get_vec2d("Gravity",m_param.gravity.v,"Gravity vector");
	config.get_bool("VolumeCorrection", m_param.volume_correction,"Should perform volume correction");
	config.get_double("VolumeChangeTolRatio",m_param.volume_change_tol_ratio,"Volume change tolerance ratio");
	config.get_double("SurfaceTension",m_param.surftens_k,"Surface tenstion coefficient");
	config.get_bool("ShowGraph",m_param.show_graph,"Show graph");
	//
	config.get_unsigned("ResolutionX",m_shape[0],"Resolution towards X axis");
	config.get_unsigned("ResolutionY",m_shape[1],"Resolution towards Y axis");
	//
	double view_scale (1.0);
	config.get_double("ViewScale",view_scale,"View scale");
	//
	double resolution_scale (1.0);
	config.get_double("ResolutionScale",resolution_scale,"Resolution doubling scale");
	//
	m_shape *= resolution_scale;
	m_dx = view_scale / m_shape[0];
}
//
void macliquid2::setup_window( std::string &name, int &width, int &height ) const {
	double ratio = m_shape[1] / (double)m_shape[0];
	height = width * ratio;
}
//
void macliquid2::post_initialize ( bool initialized_from_file ) {
	//
	auto initialize_func = reinterpret_cast<void(*)(const shape2 &m_shape, double m_dx)>(m_dylib.load_symbol("initialize"));
	if( initialize_func ) {
		initialize_func(m_shape,m_dx);
	}
	//
	// Get functions
	m_check_inject_func = reinterpret_cast<bool(*)(double, double, double, unsigned)>(m_dylib.load_symbol("check_inject"));
	m_inject_func = reinterpret_cast<bool(*)(const vec2d &, double, double, double, unsigned, double &, vec2d &)>(m_dylib.load_symbol("inject"));
	m_post_inject_func = reinterpret_cast<void(*)(double, double, double, unsigned, double&)>(m_dylib.load_symbol("post_inject"));
	m_gravity_func = reinterpret_cast<vec2d(*)(double)>(m_dylib.load_symbol("gravity"));
	//
	if( initialized_from_file ) {
		m_fluid.flood_fill();
		m_solid.flood_fill();
	} else {
		//
		m_velocity.initialize(m_shape);
		m_solid_velocity.initialize(m_shape);
		m_external_force.initialize(m_shape);
		m_solid.initialize(m_shape.nodal());
		m_fluid.initialize(m_shape.cell());
		//
		// Initialize arrays
		m_force_exist = false;
		//
		// Assign initial variables from script
		m_macutility->assign_initial_variables(m_dylib,&m_solid,&m_solid_velocity,&m_fluid,&m_velocity);
		//
		// Compute the initial volume
		m_target_volume = m_gridutility->get_area(m_solid,m_fluid);
		//
		double max_u = m_macutility->compute_max_u(m_velocity);
		if( max_u ) {
			//
			// Project to make sure that the velocity field is divergence free at the beggining
			const double CFL = m_timestepper->get_target_CFL();
			m_macproject->project(CFL*m_dx/max_u,m_velocity,m_solid,m_solid_velocity,m_fluid);
		}
	}
	//
	m_camera->set_bounding_box(vec2d().v,m_shape.box(m_dx).v);
	//
	if( m_param.show_graph ) {
		m_graphplotter->clear();
		if( m_param.gravity.norm2() ) m_graph_lists[0] = m_graphplotter->create_entry("Gravitational Energy");
		m_graph_lists[1] = m_graphplotter->create_entry("Kinetic Energy");
		if( m_param.surftens_k ) m_graph_lists[2] = m_graphplotter->create_entry("Surface Area Energy");
		m_graph_lists[3] = m_graphplotter->create_entry("Total Energy");
	}
}
//
void macliquid2::drag( double x, double y, double z, double u, double v, double w ) {
	//
	double scale (1e3);
	m_macutility->add_force(vec2d(x,y),scale*vec2d(u,v),m_external_force);
	m_force_exist = true;
}
//
void macliquid2::inject_external_force( macarray2<Real> &velocity, double dt, bool clear ) {
	//
	if( m_force_exist ) {
		velocity += m_external_force;
		if( clear ) {
			m_force_exist = false;
			m_external_force.clear();
		}
	}
	// Add gravity force
	velocity += dt*m_param.gravity;
}
//
void macliquid2::inject_external_fluid( array2<Real> &fluid, macarray2<Real> &velocity, double dt, double time ) {
	//
	const unsigned step = m_timestepper->get_step_count();
	//
	if( m_check_inject_func && m_check_inject_func(m_dx,dt,time,step)) {
		size_t total_injected (0);
		if( m_inject_func ) {
			total_injected = do_inject_external_fluid(fluid,velocity,dt,time,step);
			double volume_change = (m_dx*m_dx) * total_injected;
			if( m_post_inject_func ) {
				m_post_inject_func(m_dx,dt,time,step,volume_change);
			}
			if( volume_change ) {
				m_target_volume += volume_change;
			}
		}
	}
}
//
size_t macliquid2::do_inject_external_fluid( array2<Real> &fluid, macarray2<Real> &velocity, double dt, double time, unsigned step ) {
	//
	auto interp_vel = [&]( const vec2d &p ) {
		return macarray_interpolator2::interpolate(velocity,vec2d(),m_dx,p);
	};
	//
	std::vector<size_t> inject_count(fluid.get_thread_num(),0);
	std::vector<std::vector<vec2i> > injected_positions(fluid.get_thread_num());
	double current_CFL = m_timestepper->get_current_CFL();
	//
	size_t total_injected (0);
	fluid.parallel_all([&]( int i, int j, auto &it, int tid ) {
		//
		vec2d p = m_dx*vec2i(i,j).cell();
		double value (0.0); vec2d u;
		if( m_inject_func(p,m_dx,dt,time,step,value,u)) {
			if( value < 0.0 ) {
				injected_positions[tid].push_back(vec2i(i,j));
				if( it() >= 0.0 ) inject_count[tid] ++;
			}
			if( std::abs(value) < fluid.get_background_value() ||
				(value < fluid.get_background_value() && it.active())) {
				it.set(std::min((Real)value,it()));
			}
		}
	});
	fluid.flood_fill();
	total_injected = std::accumulate(inject_count.begin(),inject_count.end(),0);
	//
	shared_bitarray2 eval_cells(fluid.shape());
	for( const auto &e : injected_positions ) for( const auto &pi : e ) {
		eval_cells->set(pi);
	}
	eval_cells->dilate(1);
	eval_cells->const_serial_actives([&]( int i, int j ) {
		double value (0.0); vec2d u;
		if( m_inject_func(m_dx*vec2i(i,j).cell(),m_dx,dt,time,step,value,u)) {
			for( int dim : DIMS2 ) {
				vec2d p0 = m_dx*vec2i(i,j).face(dim);
				vec2d p1 = m_dx*vec2i(i+dim==0,j+dim==1).face(dim);
				m_inject_func(p0,m_dx,dt,time,step,value,u);
				velocity[dim].set(i,j,u[dim]); u = vec2d();
				m_inject_func(p1,m_dx,dt,time,step,value,u);
				velocity[dim].set(i+dim==0,j+dim==1,u[dim]); u = vec2d();
			}
		}
	});
	return total_injected;
}
//
void macliquid2::set_volume_correction( macproject2_interface *macproject ) {
	//
	// Set volume correction if requested
	if( m_param.volume_correction ) {
		double volume = m_gridutility->get_area(m_solid,m_fluid);
		if( std::abs(1.0-volume/m_target_volume) > m_param.volume_change_tol_ratio ) {
			double target_volume;
			if( volume > m_target_volume ) {
				target_volume = (1.0+m_param.volume_change_tol_ratio) * m_target_volume;
			} else {
				target_volume = (1.0-m_param.volume_change_tol_ratio) * m_target_volume;
			}
			macproject->set_target_volume(volume,target_volume);
		}
	}
}
//
void macliquid2::extend_both( int w ) {
	//
	unsigned width = w+m_timestepper->get_current_CFL();
	macarray_extrapolator2::extrapolate<Real>(m_velocity,width);
	m_macutility->constrain_velocity(m_solid,m_velocity);
	m_fluid.dilate(width);
}
//
void macliquid2::idle() {
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
	// Extend both the velocity field and the level set
	extend_both();
	//
	// Advect surface
	m_macsurfacetracker->advect(m_fluid,m_solid,m_velocity,dt);
	//
	// Advect velocity
	shared_macarray2<Real> velocity_save(m_velocity);
	m_macadvection->advect_vector(m_velocity,velocity_save(),m_fluid,dt);
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
	m_macproject->project(dt,m_velocity,m_solid,m_solid_velocity,m_fluid,m_param.surftens_k);
	//
	// Report stats
	m_macstats->dump_stats(m_solid,m_fluid,m_velocity,m_timestepper.get());
}
//
void macliquid2::add_to_graph() {
	//
	if( m_param.show_graph ) {
		//
		// Compute total energy
		const double time = m_timestepper->get_current_time();
		const auto energy_list = m_macutility->get_all_kinds_of_energy(m_solid,m_fluid,m_velocity,m_param.gravity,m_param.surftens_k);
		const double total_energy = std::get<0>(energy_list)+std::get<1>(energy_list)+std::get<2>(energy_list);
		//
		// Add to graph
		if( m_param.gravity.norm2() ) m_graphplotter->add_point(m_graph_lists[0],time,std::get<0>(energy_list));
		m_graphplotter->add_point(m_graph_lists[1],time,std::get<1>(energy_list));
		if( m_param.surftens_k ) m_graphplotter->add_point(m_graph_lists[2],time,std::get<2>(energy_list));
		m_graphplotter->add_point(m_graph_lists[3],time,total_energy);
	}
}
//
void macliquid2::draw( graphics_engine &g ) const {
	//
	// Draw solid levelset
	shared_array2<Real> solid_to_visualize(m_solid.shape());
	m_gridutility->assign_visualizable_solid(m_dylib,m_dx,solid_to_visualize());
	m_gridvisualizer->draw_solid(g,solid_to_visualize());
	//
	// Visualize moving solid
	auto draw_func = reinterpret_cast<void(*)(graphics_engine &,double)>(m_dylib.load_symbol("draw"));
	if( draw_func ) {
		g.color4(1.0,0.8,0.5,0.3);
		draw_func(g,m_timestepper->get_current_time());
	}
	//
	// Draw grid lines
	m_gridvisualizer->draw_grid(g);
	//
	// Draw fluid
	m_gridvisualizer->draw_fluid(g,m_solid,m_fluid);
	//
	// Draw projection component
	m_macproject->draw(g);
	//
	// Draw velocity
	m_macvisualizer->draw_velocity(g,m_velocity);
	//
	// Draw graph
	m_graphplotter->draw(g);
}
//
extern "C" module * create_instance() {
	return new macliquid2;
}
//
extern "C" const char *license() {
	return "MIT";
}