/*
**	macsmoke2_oc.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on April 10, 2017.
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
#include "macsmoke2_oc.h"
#include <shiokaze/core/filesystem.h>
#include <shiokaze/graphics/graphics_utility.h>
#include <shiokaze/array/shared_array2.h>
#include <shiokaze/array/array_derivative2.h>
#include <shiokaze/array/array_interpolator2.h>
#include <shiokaze/array/macarray_interpolator2.h>
#include <cmath>
#include <random>
//
SHKZ_USING_NAMESPACE
//
macsmoke2_oc::macsmoke2_oc () {
	//
	m_shape = shape2(64,32);
	m_dx = m_shape.dx();
}
//
void macsmoke2_oc::setup_window( std::string &name, int &width, int &height ) const {
	//
	double ratio = m_shape[1] / (double)m_shape[0];
	height = width * ratio;
}
//
void macsmoke2_oc::load( configuration &config ) {
	//
	std::string name("plume2"); config.get_string("Name",name,"Scene file name");
	m_dylib.open_library(filesystem::resolve_libname(name));
	m_dylib.load(config);
	m_dylib.overwrite(config);
}
//
void macsmoke2_oc::configure( configuration &config ) {
	//
	m_dylib.configure(config);
	//
	config.get_bool("UseDustParticles",m_param.use_dust,"Whether to use dust particles instead of density field");
	if( m_param.use_dust ) {
		config.get_unsigned("DustSampleNum",m_param.r_sample,"Subsampling number for dust particles per dimension divided by 2");
	} else {
		config.get_double("MinimalActiveDensity",m_param.minimal_density,"Minimal density to trim active cells");
	}
	config.get_bool("ShowGraph",m_param.show_graph,"Show graph");
	config.get_double("BuoyancyFactor",m_param.buoyancy_factor,"Buoyancy force rate");
	config.get_unsigned("SolidExtrapolationDepth",m_param.extrapolated_width,"Solid extrapolation depth");
	config.get_unsigned("ResolutionX",m_shape[0],"Resolution towards X axis");
	config.get_unsigned("ResolutionY",m_shape[1],"Resolution towards Y axis");
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
void macsmoke2_oc::post_initialize ( bool initialized_from_file ) {
	//
	// Initialize scene
	auto initialize_func = reinterpret_cast<void(*)(const shape2 &shape, double dx)>(m_dylib.load_symbol("initialize"));
	if( initialize_func ) initialize_func(m_shape,m_dx);
	//
	// Get functions
	m_solid_func = reinterpret_cast<double(*)(const vec2d &)>(m_dylib.load_symbol("solid"));
	m_draw_func = reinterpret_cast<void(*)(graphics_engine &,double)>(m_dylib.load_symbol("draw"));
	m_moving_solid_func = reinterpret_cast<std::pair<double,vec2d>(*)(double time, const vec2d &p)>(m_dylib.load_symbol("moving_solid"));
	
	// m_check_inject_func = reinterpret_cast<bool(*)(double, double, double, unsigned)>(m_dylib.load_symbol("check_inject"));
	// m_inject_func = reinterpret_cast<bool(*)(const vec2d &, double, double, double, unsigned, double &, vec2d &)>(m_dylib.load_symbol("inject"));
	// m_post_inject_func = reinterpret_cast<void(*)(double, double, double, unsigned, double&)>(m_dylib.load_symbol("post_inject"));
	// auto velocity_func = reinterpret_cast<vec2d(*)(const vec2d &)>(m_dylib.load_symbol("velocity"));
	// m_gravity_func = reinterpret_cast<vec2d(*)(double)>(m_dylib.load_symbol("gravity"));
	m_set_boundary_flux = reinterpret_cast<void(*)( double, Real [DIM2][2] )>(m_dylib.load_symbol("set_boundary_flux"));
	m_combined_solid_func = [&]( const vec2d &p ) {
		double value (1.0);
		if( m_solid_func ) value = std::min(value,m_solid_func(p));
		if( m_moving_solid_func ) value = std::min(value,m_moving_solid_func(m_timestepper->get_current_time(),p).first);
		return value;
	};
	if( m_moving_solid_func ) {
		m_macoctreeproject.set_moving_solid([&]( const vec2d &p) {
			return m_moving_solid_func(m_timestepper->get_current_time(),p).first;
		});
	}
	
	//ここはmacoctreeliquidにあってmacsmokeにない部分。比較のために書き起こしてみた
	//macoctreeproject.set_moving_solidに引数を入れて何かをしているみたい。
	
	//added
	m_accumulated_CFL = 0.0;
	m_grid_0.clear();
	m_grid_1.clear();
	// m_macoctreehelper.initialize();
	//octreehelperというもので何かを初期化している（octreeの構造？）

	//added



	
	// Initialize arrays
	m_force_exist = false;
	m_velocity.initialize(m_shape);
	m_solid_velocity.initialize(m_shape);
	m_external_force.initialize(m_shape);
	//
	m_solid.initialize(m_shape.nodal());
	m_fluid.initialize(m_shape.cell(),-1.0);
	m_density.initialize(m_shape.cell(),0.0);
	//
	if( m_param.use_dust ) {
		m_accumulation.initialize(m_shape.cell(),0.0);
	}
	m_dust_particles.clear();
	//
	// Assign initial variables from scene
	m_macutility->assign_initial_variables(m_dylib,&m_solid,&m_solid_velocity,nullptr,&m_velocity,&m_density);
	//
    unsigned n = 1;
    while( (m_shape/n).min() >= m_param.min_resolution ) {
        shape2 shape = m_shape/n;
        m_grid_0.add_layer(shape,n*m_dx);
        m_grid_1.add_layer(shape,n*m_dx);
        n *= 2;
    }

	int refinement_count = m_param.use_sizing_func ? m_param.initial_refinement : 1;
    while( refinement_count-- ) {
        std::swap(m_grid,m_grid_prev);

        if( m_param.use_sizing_func ) {
            if( refinement_count ) {
                m_grid->activate_cells([&](char depth, const vec2d &p) {
                    return depth > 3;
                });
            }
        } else {
            // 密度に基づいてセルを活性化
            m_grid->activate_cells([&](const vec2d &p) {
                return array_interpolator2::interpolate<Real>(m_density,p/m_dx) > m_param.minimal_density;
            }, m_combined_solid_func);
        }

        m_grid->balance_layers();
        m_grid->assign_indices();
        m_grid->assign_levelset([&](const vec2d &p) {
            return array_interpolator2::interpolate<Real>(m_density,p/m_dx);
        }, m_combined_solid_func);
    }
	// Project to make sure that the velocity field is divergence free at the beggining
	double max_u = m_macutility->compute_max_u(m_velocity);
	if( max_u ) {
		double CFL = m_timestepper->get_target_CFL();
		m_macproject->project(CFL*m_dx/max_u,m_velocity,m_solid,m_solid_velocity,m_fluid);
	}
	//
	// Seed dust particles if requested
	if( m_param.use_dust ) {
		//
		shared_array2<Real> density_copy(m_density);
		density_copy->dilate();
		//
		double space = 1.0 / m_param.r_sample;
		density_copy->const_serial_actives([&]( int i, int j, const auto &it ) {
			for( int ii=0; ii<m_param.r_sample; ++ii ) for( int pjj=0; pjj<m_param.r_sample; ++pjj ) {
				int jj = ii % 2 == 0 ? pjj : m_param.r_sample-pjj-1;
				vec2d unit_pos = 0.5*vec2d(space,space)+vec2d(ii*space,jj*space);
				vec2d pos = m_dx*(unit_pos+vec2d(i,j));
				if( array_interpolator2::interpolate<Real>(m_solid,pos/m_dx) > 0.0 && 
					array_interpolator2::interpolate<Real>(m_density,pos/m_dx-vec2d(0.5,0.5))) {
					m_dust_particles.push_back(pos);
				}
			}
		});
		rasterize_dust_particles(m_density);
	}
	//
	m_camera->set_bounding_box(vec2d().v,m_shape.box(m_dx).v);
	//
	if( m_param.show_graph ) {
		m_graphplotter->clear();
		m_graph_id = m_graphplotter->create_entry("Kinetic Energy");
	}
}
//
void macsmoke2_oc::drag( double x, double y, double z, double u, double v, double w ) {
	//
	double scale (1e3);
	m_macutility->add_force(vec2d(x,y),scale*vec2d(u,v),m_external_force);
	m_force_exist = true;
}
//
void macsmoke2_oc::inject_external_force( macarray2<Real> &velocity ) {
	//
	if( m_force_exist ) {
		velocity += m_external_force;
		m_external_force.clear();
		m_force_exist = false;
	}
}
//
void macsmoke2_oc::add_source ( macarray2<Real> &velocity, array2<Real> &density, double time, double dt ) {
	auto add_func = reinterpret_cast<void(*)(const vec2d &, vec2d &, double &, double, double)>(m_dylib.load_symbol("add"));
	if( add_func ) {
		//
		// Velocity
		velocity.parallel_all([&](int dim, int i, int j, auto &it) {
			vec2d p = m_dx*vec2i(i,j).face(dim);
			double dummy; vec2d u;
			add_func (p,u,dummy,time,dt);
			if( u[dim] ) it.increment(u[dim]);
		});
		//
		// Density
		auto add_density = [&]( array2<Real> &density ) {
			density.parallel_all([&](int i, int j, auto &it) {
				vec2d p = m_dx*vec2i(i,j).cell();
				double d(0.0); vec2d dummy;
				add_func (p,dummy,d,time,dt);
				if( d ) it.increment(d);
			});
		};
		//
		if( m_param.use_dust ) {
			//
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_real_distribution<> dis(-1.0,1.0);
			//
			add_density(m_accumulation);
			//
			double scale = 1.0 / pow(m_param.r_sample,DIM2);
			bool should_re_rasterize (false);
			m_accumulation.serial_op([&]( int i, int j, auto &it) {
				double d = it();
				while( d > scale ) {
					vec2d p = m_dx*vec2i(i,j).cell()+0.5*m_dx*vec2d(dis(gen),dis(gen));
					m_dust_particles.push_back(p);
					should_re_rasterize = true;
					d -= scale;
				}
				it.set(d);
			});
			//
			if( should_re_rasterize ) {
				rasterize_dust_particles(density);
			}
			//
		} else {
			add_density(density);
		}
	}
}
//
void macsmoke2_oc::rasterize_dust_particles( array2<Real> &rasterized_density ) {
	//
	rasterized_density.clear();
	double scale = 1.0 / pow(m_param.r_sample,DIM2);
	for( const vec2d &p : m_dust_particles ) {
		vec2i pi = p/m_dx;
		if( ! rasterized_density.shape().out_of_bounds(pi)) {
			rasterized_density.increment(pi,scale);
		}
	}
}
//
void macsmoke2_oc::add_buoyancy_force( macarray2<Real> &velocity, const array2<Real> &density, double dt ) {
	//
	velocity[1].parallel_all([&]( int i, int j, auto &it, int tn ) {
		vec2i pi = vec2i(i,j).face(1);
		Real d = array_interpolator2::interpolate<Real>(density,(pi-vec2d(0.5,0.5)));
		it.increment(m_param.buoyancy_factor*dt*d);
	});
}
//
void macsmoke2_oc::idle() {
	//
	// Add to graph
	add_to_graph();
	//
	// Compute the timestep size
	const double dt = m_timestepper->advance(m_macutility->compute_max_u(m_velocity),m_dx);
	const double time = m_timestepper->get_current_time();
	const double CFL = m_timestepper->get_current_CFL();
	//added

    // octreeグリッドのリメッシュ判定
    m_accumulated_CFL += CFL;
    std::swap(m_grid,m_grid_prev);
    if( m_accumulated_CFL >= m_param.maximal_CFL_accumulation ) {
        m_accumulated_CFL = 0.0;

        if( m_param.use_sizing_func ) {
			//ここでactivatecellsの引数がて不適切らしい。->よく考えたらactivate cellsを行う基準を考えてなかったのでそれはそうかも→あくまで練習だし適当に決めてもいいんじゃないかな～～～
			//解決！
           m_grid->activate_cells([&](char depth, const vec2d &p) {
            double density = array_interpolator2::interpolate<Real>(m_density,p/m_dx);
            return density > m_param.minimal_density;
        	});
        } else {
            m_grid->activate_cells([&](const vec2d &p) {
                return array_interpolator2::interpolate<Real>(m_density,p/m_dx) > m_param.minimal_density;
            }, m_combined_solid_func);
        }

        m_grid->balance_layers();
        m_grid->assign_indices();
    } else {
        m_grid->copy(*m_grid_prev);
    }


    m_grid->assign_levelset([&](const vec2d &p) {
        vec2d u (m_grid_prev->sample_velocity(p));
        double d = array_interpolator2::interpolate<Real>(m_density,p/m_dx);
        return d > m_param.minimal_density ? -1.0 : 1.0;
    }, m_combined_solid_func);
	//added

	// Update solid
	m_macutility->update_solid_variables(m_dylib,time,&m_solid,&m_solid_velocity);
	//
	// Advect density and velocity
	if( m_param.use_dust ) advect_dust_particles(m_velocity,dt);
	else {
		m_density.dilate(std::ceil(m_timestepper->get_current_CFL()));
		m_macadvection->advect_scalar(m_density,m_velocity,m_fluid,dt);
		m_density.parallel_actives([&](auto &it) {
			if( std::abs(it()) <= m_param.minimal_density ) it.set_off();
		});
	}
	//
	shared_macarray2<Real> velocity_save(m_velocity);
	m_macadvection->advect_vector(m_velocity,velocity_save(),m_fluid,dt);
	//
	// Add buoyancy force
	add_buoyancy_force(m_velocity,m_density,dt);
	//
	// Add source
	add_source(m_velocity,m_density,m_timestepper->get_current_time(),dt);
	//
	// Add external force
	inject_external_force(m_velocity);
	//
	// Project
	m_macproject->project(dt,m_velocity,m_solid,m_solid_velocity,m_fluid,0.0);
	m_macutility->extrapolate_and_constrain_velocity(m_solid,m_velocity,m_param.extrapolated_width);
	//
	// Report stats
	m_macstats->dump_stats(m_solid,m_fluid,m_velocity,m_timestepper.get());
}
//
void macsmoke2_oc::advect_dust_particles( const macarray2<Real> &velocity, double dt ) {
	//
	m_parallel.for_each( m_dust_particles.size(), [&]( size_t n, int tn ) {
		vec2d &p = m_dust_particles[n];
		vec2d u0 = macarray_interpolator2::interpolate<Real>(velocity,p/m_dx);
		vec2d u1 =  macarray_interpolator2::interpolate<Real>(velocity,(p+dt*u0)/m_dx);
		p += 0.5 * dt * (u0+u1);
	});
	//
	std::vector<char> remove_flags (m_dust_particles.size());
	m_parallel.for_each( m_dust_particles.size(), [&]( size_t n, int tn ) {
		vec2d &p = m_dust_particles[n];
		Real phi = array_interpolator2::interpolate<Real>(m_solid,p/m_dx);
		if( phi < 0.0 ) {
			Real derivative[DIM2];
			array_derivative2::derivative(m_solid,p/m_dx,derivative);
			p = p - phi*vec2d(derivative).normal();
		}
		for( unsigned dim : DIMS2 ) {
			if( p[dim] < 0.0 || p[dim] > m_dx*m_shape[dim] ) remove_flags[n] = true;
		}
	});
	//
	std::vector<vec2d> save_dust_particles (m_dust_particles);
	m_dust_particles.clear();
	for( size_t n=0; n<save_dust_particles.size(); ++n ) {
		if( ! remove_flags[n] ) m_dust_particles.push_back(save_dust_particles[n]);
	}
	rasterize_dust_particles(m_density);
}
//
void macsmoke2_oc::draw_dust_particles( graphics_engine &g ) const {
	//
	using ge = graphics_engine;
	double r = m_dx * 0.5 / m_param.r_sample;
	for( const vec2d &p : m_dust_particles ) {
		g.color4(1.0,1.0,1.0,1.0);
		graphics_utility::draw_circle(g,p.v,r,ge::MODE::LINE_LOOP);
		g.color4(1.0,1.0,1.0,0.3);
		graphics_utility::draw_circle(g,p.v,r,ge::MODE::TRIANGLE_FAN);
	}
}
//
void macsmoke2_oc::add_to_graph() {
	//
	if( m_param.show_graph ) {
		//
		// Compute total energy
		const double time = m_timestepper->get_current_time();
		const double total_energy = m_macutility->get_kinetic_energy(m_solid,m_fluid,m_velocity);
		//
		// Add to graph
		m_graphplotter->add_point(m_graph_id,time,total_energy);
	}
}
//
void macsmoke2_oc::draw( graphics_engine &g ) const {
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
	// Draw grid edges
	//m_gridvisualizer->draw_grid(g);
	//

    m_grid->draw_grid(g);
    //m_grid->draw_fluid(g);
	//
	// Draw projection component
	// m_macproject->draw(g);
	//
	// Draw velocity
	// m_macvisualizer->draw_velocity(g,m_velocity);
	//
	// Draw density
	if( m_param.use_dust ) draw_dust_particles(g);
	else m_gridvisualizer->draw_density(g,m_density);
	//
	// Draw graph
	//m_graphplotter->draw(g);
}
//
extern "C" module * create_instance() {
	return new macsmoke2_oc;
}
//
extern "C" const char *license() {
	return "MIT";
}
//