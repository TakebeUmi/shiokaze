/*
**	cloud_oc.cpp
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
#include "cloud_oc.h"
#include <openvdb/openvdb.h>
#include <shiokaze/core/console.h>
#include <shiokaze/core/timer.h>
#include <shiokaze/core/filesystem.h>
#include <shiokaze/array/shared_array3.h>
#include <shiokaze/array/array_derivative3.h>
#include <shiokaze/array/array_interpolator3.h>
#include <shiokaze/array/macarray_interpolator3.h>
#include <shiokaze/utility/utility.h>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <iomanip>

//
SHKZ_USING_NAMESPACE
//

cloud_oc::cloud_oc () {
	//
	m_shape = shape3(64,64,64);
	//m_dx = m_shape.dx();
}
//
void cloud_oc::setup_window( std::string &name, int &width, int &height ) const {
	height = width;
}
//
void cloud_oc::load( configuration &config ) {
	//
	std::string name("plume3"); config.get_string("Name",name,"Scene file name");
	m_dylib.open_library(filesystem::resolve_libname(name));
	m_dylib.load(config);
	m_dylib.overwrite(config);
	//
	m_param.render_density = console::system("mitsuba > /dev/null 2>&1") == 0;
	config.get_bool("RenderDensity",m_param.render_density,"Whether to render density");
}

//added

void cloud_oc::write_to_txt(std::string type) const {
	std::time_t t  = std::time(nullptr);
	std::tm* now = std::localtime(&t);

	std::ostringstream oss;
    oss << "density_output_"
        << type << "_"
        << std::put_time(now, "%Y-%m-%d_%H-%M-%S") // ← 日時フォーマット
        << ".txt";

	std::string filename = oss.str();
	std::ofstream outfile(filename);
	if (outfile.is_open()) {
		for (size_t i = 0; i < m_grid->density.size(); ++i) {
			
		char buffer[64];
		std::snprintf(buffer, sizeof(buffer), "%04zu:%.6f ", i, m_grid->density[i]);
		outfile << buffer << " ";
			if ((i+1) % m_shape[0] == 0) outfile << "\n"; // New line for each z-slice (assuming 64x64x64 grid)
		}
		outfile.close();
		console::dump("Density data written to %s\n", filename.c_str());
	} else {
		console::dump("Failed to open file %s for writing.\n", filename.c_str());
	}
}

//
bool cloud_oc::should_quit() const {
	return m_should_quit_on_save || m_timestepper->should_quit();
}
//added
//
void cloud_oc::configure( configuration &config ) {
	//
	// Configure the set of tools
	m_dylib.configure(config);
	//
	config.get_bool("UseDustParticles",m_param.use_dust,"Whether to use dust particles instead of density field");
	if( m_param.use_dust ) {
		config.get_unsigned("DustSampleNum",m_param.r_sample,"Subsampling number for dust particles per dimension divided by 2");
	} else {
		config.get_double("MinimalActiveDensity",m_param.minimal_density,"Minimal density to trim active cells");
	}
	config.get_bool("MouseInteration",m_param.mouse_interaction, "Enable mouse interaction");
	config.get_bool("ShowGraph",m_param.show_graph,"Show graph");
	config.get_double("BuoyancyFactor",m_param.buoyancy_factor,"Buoyancy force rate");
	config.get_unsigned("SolidExtrapolationDepth",m_param.extrapolated_width,"Solid extrapolation depth");
	config.get_unsigned("ResolutionX",m_shape[0],"Resolution towards X axis");
	config.get_unsigned("ResolutionY",m_shape[1],"Resolution towards Y axis");
	config.get_unsigned("ResolutionZ",m_shape[2],"Resolution towards Z axis");
	config.get_unsigned("RenderSampleCount",m_param.render_sample_count,"Sample count for rendering");
	config.get_double("VolumeScale",m_param.volume_scale,"Volume scaling for rendering");
	//
	//added
	config.get_vec3d("Gravity",m_param.gravity.v,"Gravity vector");
	config.get_bool("UseFLIP",m_param.use_FLIP,"Whether to use FLIP");
	config.get_bool("RenderMesh",m_param.render_mesh,"Whether to render mesh files");
	config.get_bool("RenderWireframe",m_param.render_wireframe,"Whether to render wireframe view");
	config.get_bool("RenderGrid",m_param.render_grid,"Render grid");
	config.get_bool("RemoveQuarter",m_param.remove_quater,"Remove front-left meshes");
	config.get_double("ZPosition",m_param.z,"Z coordinate position");
	config.get_bool("ExportSVG",m_param.export_svg,"Export cutaway SVG");
	config.get_unsigned("MinResolution",m_param.min_resolution,"Minimal resolution");
	config.get_unsigned("ErodeWidth",m_param.erode_width,"Erosion width");
	config.get_double("SurfaceTension",m_param.surftens_k,"Surface tension coefficient");
	config.get_unsigned("ResolutionX",m_shape[0],"Resolution towards X axis");
	config.get_unsigned("ResolutionY",m_shape[1],"Resolution towards Y axis");
	config.get_unsigned("ResolutionZ",m_shape[2],"Resolution towards Z axis");
	config.get_vec3d("TargetPos",m_param.target.v,"Camera target position");
	config.get_vec3d("OriginPos",m_param.origin.v,"Camera origin position");
	config.get_double("PICFLIP",m_param.PICFLIP,"PICFLIP blending factor");
	config.get_bool("RenderTransparent",m_param.render_transparent,"Whether to render transparent view");
	config.get_unsigned("RenderSampleCount",m_param.render_sample_count,"Sample count for rendering");
	config.get_unsigned("RenderTransparentSampleCount",m_param.render_transparent_sample_count,"Sample count for transparent rendering");
	config.get_bool("TransferFile",m_param.transfer_file,"Tranfer file via rsync");
	config.get_unsigned("SaveInterval",m_param.save_interval,"Saving state interval time steps");
	config.get_bool("VolumeCorrection",m_param.volume_correction,"Whether to perform volume correction");
	config.get_bool("RegionalVolumeCorrection",m_param.regional_volume_correction,"Regional volume correction");
	config.get_bool("MacCormack",m_param.maccormack,"Use MacCormack method");
	config.get_bool("UseSizingFunc",m_param.use_sizing_func,"Use sizing function");
	config.get_unsigned("InitialRefinement",m_param.initial_refinement,"Initial refinement count");
	config.get_double("MaxCFLAccumulation",m_param.maximal_CFL_accumulation,"CFL sum trigger for remeshing");
	config.get_integer("DebugMode",m_param.debug_mode,"Debug mode");
	//
	if( m_param.use_sizing_func ) {
		config.set_default_bool("SteepAdaptivity",false);
		config.set_default_unsigned("DilateCount",2);
	}
	//added
	double view_scale (1.0);
	config.get_double("ViewScale",view_scale,"View scale");
	//
	double resolution_scale (1.0);
	config.get_double("ResolutionScale",resolution_scale,"Resolution doubling scale");
	//
	m_shape *= resolution_scale;
	m_dx = view_scale * m_shape.dx();
	//console::dump("Computed dx: %e\n",m_dx);
	//added
	unsigned solid_max_resolution (512);
	config.get_unsigned("SolidMaxResolution",solid_max_resolution,"Max resolution for solid level set for visualization");
	m_solid_shape = m_shape;
	m_solid_dx = m_dx;
	while( m_solid_shape.max() > solid_max_resolution ) {
		m_solid_shape = m_solid_shape/2;
		m_solid_dx = 2.0 * m_solid_dx;
	}
	// m_solid_gridvisualizer->set_environment("shape",&m_solid_shape);
	// m_solid_gridvisualizer->set_environment("dx",&m_solid_dx);
	//
	m_solid_mesher->set_environment("shape",&m_solid_shape);
	m_solid_mesher->set_environment("dx",&m_solid_dx);
	//added
}
//
void cloud_oc::post_initialize ( bool initialized_from_file ) {
	//
	scoped_timer timer(this);
	//
	timer.tick(); console::dump( ">>> Started initialization (%dx%dx%d)\n", m_shape[0], m_shape[1], m_shape[2] );
	//
	auto initialize_func = reinterpret_cast<void(*)(const shape3 &shape, double dx)>(m_dylib.load_symbol("initialize"));
	if( initialize_func ) initialize_func(m_shape,m_dx);
	//
	// Get functions
	m_set_boundary_flux = reinterpret_cast<void(*)( double, Real [DIM3][2] )>(m_dylib.load_symbol("set_boundary_flux"));
	m_draw_func = reinterpret_cast<void(*)(graphics_engine &,double)>(m_dylib.load_symbol("draw"));
	m_moving_solid_func = reinterpret_cast<std::pair<double,vec3d>(*)(double time, const vec3d &p)>(m_dylib.load_symbol("moving_solid"));
	//
	//added
	//m_inject_func = reinterpret_cast<bool(*)(const vec3d &, double, double, double, unsigned, double &, vec3d &)>(m_dylib.load_symbol("inject"));
	auto fluid_func = reinterpret_cast<double(*)(const vec3d &)>(m_dylib.load_symbol("fluid"));
	m_solid_func = reinterpret_cast<double(*)(const vec3d &)>(m_dylib.load_symbol("solid"));
	//added
	auto velocity_func = reinterpret_cast<vec3d(*)(const vec3d &)>(m_dylib.load_symbol("velocity"));
	m_combined_solid_func = [&]( const vec3d &p ) {
		double value (1.0);
		if( m_solid_func ) value = std::min(value,m_solid_func(p));
		if( m_moving_solid_func ) value = std::min(value,m_moving_solid_func(m_timestepper->get_current_time(),p).first);
		return value;
	};
	if( m_moving_solid_func ) {
		m_macoctreeproject.set_moving_solid([&]( const vec3d &p) {
			return m_moving_solid_func(m_timestepper->get_current_time(),p).first;
		});
	}
	m_accumulated_CFL = 0.0;
	m_grid_0.clear();
	m_grid_1.clear();
	//added


	// Initialize arrays
	// m_force_exist = false;
	// m_velocity.initialize(m_shape);
	// m_solid_velocity.initialize(m_shape);
	// m_external_force.initialize(m_shape);
	//
	// m_solid.initialize(m_shape.nodal());
	// m_fluid.initialize(m_shape.cell(),-1.0);
	// m_density.initialize(m_shape.cell(),0.0);
	//
	// if( m_param.use_dust ) {
	// 	m_accumulation.initialize(m_shape.cell(),0.0);
	// }
	m_dust_particles.clear();
	//
	// Assign initial variables from script
	//m_macutility->assign_initial_variables(m_dylib,&m_solid,&m_solid_velocity,nullptr,&m_velocity,&m_density);
	//
    unsigned n = 1;
    while( (m_shape/n).min() >= m_param.min_resolution ) {
        shape3 shape = m_shape/n;
        m_grid_0.add_layer(shape,n*m_dx);
        m_grid_1.add_layer(shape,n*m_dx);
        n *= 2;
    }

	int refinement_count = m_param.use_sizing_func ? m_param.initial_refinement : 1;
	int count (0);
		while( refinement_count-- ) {
			timer.tick(); console::dump( ">>> Refinement #%d...\n", count+1 );
			//
			std::swap(m_grid,m_grid_prev);
			if( m_param.use_sizing_func ) {
				if( count ) {
					m_macoctreesizingfunc.activate_cells(*m_grid_prev,*m_grid,m_combined_solid_func,nullptr);
				} else {
					m_grid->activate_cells([&](char depth, const vec3d &p) {
						return depth > 3;
					});
				}
			} else {
				m_grid->activate_cells(fluid_func,m_combined_solid_func);
			}

			m_grid->balance_layers();
			m_grid->assign_indices();
			//m_grid->assign_levelset(fluid_func,m_combined_solid_func);
			//m_grid->assign_density(fluid_func); 
			//
			//write_to_txt("assigned_density"); //added
			if( m_param.use_sizing_func ) {
				m_macoctreesizingfunc.compute_sizing_function(*m_grid_prev,*m_grid,0.0,m_combined_solid_func,[&]( const vec3d &p ) {
					return m_moving_solid_func ? m_moving_solid_func(m_timestepper->get_current_time(),p).second : vec3d();
				});
			}
			console::dump( "<<< Done. Took %s.\n", timer.stock("refinement"+std::to_string(count)).c_str());
			count++;
		}
		

	// if( m_set_boundary_flux ) {
	// 	flux_boundary_condition3 boundary_cond;
	// 	m_set_boundary_flux(0.0,boundary_cond.velocity);
	// 	m_grid_0.set_flux_boundary_condition(boundary_cond);
	// 	m_grid_1.set_flux_boundary_condition(boundary_cond);
	// }
		// int refinement_count = m_param.use_sizing_func ? m_param.initial_refinement : 1;
		// int count (0);
		// while( refinement_count-- ) {
		// 	timer.tick(); console::dump( ">>> Refinement #%d...\n", count+1 );
		// 	//
		// 	std::swap(m_grid,m_grid_prev);
		// 	if( m_param.use_sizing_func ) {
		// 		if( count ) {
		// 			m_macoctreesizingfunc.activate_cells(*m_grid_prev,*m_grid,m_combined_solid_func,nullptr);
		// 		} else {
		// 			m_grid->activate_cells([&](char depth, const vec3d &p) {
		// 				return depth > 3;
		// 			});
		// 		}
		// 	} else {
		// 		m_grid->activate_cells(fluid_func,m_combined_solid_func);
		// 	}
		// 	m_grid->balance_layers();
		// 	m_grid->assign_indices();
		// 	m_grid->assign_levelset(fluid_func,m_combined_solid_func);
		// 	m_grid->assign_density(fluid_func,m_combined_solid_func);
		// 	//
		// 	if( m_param.use_sizing_func ) {
		// 		m_macoctreesizingfunc.compute_sizing_function(*m_grid_prev,*m_grid,0.0,m_combined_solid_func,[&]( const vec3d &p ) {
		// 			return m_moving_solid_func ? m_moving_solid_func(m_timestepper->get_current_time(),p).second : vec3d();
		// 		});
		// 	}
		// 	console::dump( "<<< Done. Took %s.\n", timer.stock("refinement"+std::to_string(count)).c_str());
		// 	count ++;
		// }

	// velocity関数があれば速度場を初期化
	if( velocity_func ) {
    m_grid->set_velocity([&]( const vec3d &p, char dim ) {
        return velocity_func(p)[dim];
    	});
	}

	//added

	// Ensure divergence free
	//
	// Seed dust particles if requested
	// if( m_param.use_dust ) {
	// 	timer.tick(); console::dump( "Seeding dust particles..." );
	// 	//
	// 	shared_array3<Real> density_copy(m_density);
	// 	density_copy->dilate();
	// 	//
	// 	double space = 1.0 / m_param.r_sample;
	// 	density_copy->const_serial_actives([&]( int i, int j, int k, const auto &it ) {
	// 		for( int ii=0; ii<m_param.r_sample; ++ii ) for( int jj=0; jj<m_param.r_sample; ++jj ) for( int kk=0; kk<m_param.r_sample; ++kk ) {
	// 			vec3d unit_pos = 0.5*vec3d(space,space,space)+vec3d(ii*space,jj*space,kk*space);
	// 			vec3d pos = m_dx*(unit_pos+vec3d(i,j,k));
	// 			if( array_interpolator3::interpolate<Real>(m_solid,pos/m_dx) > 0.0 &&
	// 				array_interpolator3::interpolate<Real>(m_density,pos/m_dx-vec3d(0.5,0.5,0.5))) {
	// 				m_dust_particles.push_back(pos);
	// 			}
	// 		}
	// 	});
	//rasterize_dust_particles(m_density);
	// 	console::dump( "Done. Seeded=%d. Took %s.\n", m_dust_particles.size(), timer.stock("seed_m_dust_particles").c_str());
	// }
	//
	m_camera->set_bounding_box(vec3d().v,m_shape.box(m_dx).v);
	console::dump( "<<< Initialization finished. Took %s\n", timer.stock("initialization").c_str());
	//
	if( m_param.show_graph ) {
		m_graphplotter->clear();
		m_graph_id = m_graphplotter->create_entry("Kinetic Energy");
	}
}
//
// void cloud_oc::drag( double x, double y, double z, double u, double v, double w ) {
	
// 	if( m_param.mouse_interaction ) {
// 		double scale (1e3);
// 		m_macutility->add_force(vec3d(x,y,z),scale*vec3d(u,v,w),m_external_force);
// 		m_force_exist = true;
// 	}
// }
//
// void cloud_oc::inject_external_force( macarray3<Real> &velocity ) {
// 	//
// 	if( m_force_exist ) {
// 		velocity += m_external_force;
// 		m_external_force.clear();
// 		m_force_exist = false;
// 	}
// }
//
// void cloud_oc::add_source ( macarray3<Real> &velocity, array3<Real> &density, double time, double dt ) {
// 	//
// 	scoped_timer timer(this);
// 	//
// 	auto add_func = reinterpret_cast<void(*)(const vec3d &, vec3d &, double &, double, double)>(m_dylib.load_symbol("add"));
// 	if( add_func ) {
// 		timer.tick(); console::dump( "Adding sources..." );
// 		//
// 		// Velocity
// 		velocity.parallel_all([&](int dim, int i, int j, int k, auto &it) {
// 			vec3d p = m_dx*vec3i(i,j,k).face(dim);
// 			double dummy; vec3d u;
// 			add_func (p,u,dummy,time,dt);
// 			if( u[dim] ) it.increment(u[dim]);
// 		});
// 		//
// 		// Density
// 		auto add_density = [&]( array3<Real> &density ) {
// 			density.parallel_all([&](int i, int j, int k, auto &it) {
// 				vec3d p = m_dx*vec3i(i,j,k).cell();
// 				double d(0.0); vec3d dummy;
// 				add_func (p,dummy,d,time,dt);
// 				density.increment(i,j,k,d);
// 			});
// 		};
// 		//
// 		// Density
// 		unsigned seeded (0);
// 		if( m_param.use_dust ) {
// 			//
// 			std::random_device rd;
// 			std::mt19937 gen(rd());
// 			std::uniform_real_distribution<> dis(-1.0,1.0);
// 			//
// 			add_density(m_accumulation);
// 			//
// 			double scale = 1.0 / pow(m_param.r_sample,DIM3);
// 			bool should_re_rasterize (false);
// 			m_accumulation.serial_op([&]( int i, int j, int k, auto &it) {
// 				double d = it();
// 				while( d > scale ) {
// 					vec3d p = m_dx*vec3i(i,j,k).cell()+0.5*m_dx*vec3d(dis(gen),dis(gen),dis(gen));
// 					m_dust_particles.push_back(p); ++ seeded;
// 					should_re_rasterize = true;
// 					d -= scale;
// 				}
// 				it.set(d);
// 			});
// 			//
// 			if( should_re_rasterize ) {
// 				rasterize_dust_particles(density);
// 			}
// 			//
// 		} else {
// 			add_density(density);
// 		}
// 		//
// 		if( m_param.use_dust ) console::dump( "Done. Seeded=%d. Took %s.\n", seeded, timer.stock("add_func").c_str());
// 		else console::dump( "Done. Took %s.\n", timer.stock("add_func").c_str());
// 	}
// }
//これを

void cloud_oc::microphysics_cloud(grid3 &grid, double dt) {
	grid.iterate_active_cells([&](const cell_id3 &cell_id, int tid) {
		double qv = grid.qv[cell_id.index];
		double qc = grid.qc[cell_id.index];
		double qr = grid.qr[cell_id.index];
		double theta = grid.theta[cell_id.index];
		double newtheta, newqv, newqc, newqr;
		vec3d z = grid.get_cell_position(cell_id);
		grid.kessler_model(dt, grid.param.T0, grid.param.p0, grid.param.gamma, grid.param.g, grid.param.alphaCE, grid.param.alphaA, grid.param.alphaK, z[1], qv, qc, qr, theta, newqv, newqc, newqr, newtheta);
		grid.qv[cell_id.index] = newqv;
		grid.qc[cell_id.index] = newqc;
		grid.qr[cell_id.index] = newqr;
		grid.theta[cell_id.index] = newtheta;
	});
}

void cloud_oc::add_source_oc (double time, double dt , grid3 &grid) {
	//
	scoped_timer timer(this);
	//
	auto add_func = reinterpret_cast<void(*)(const vec3d &, vec3d &, double &, double, double)>(m_dylib.load_symbol("add"));
	if( add_func ) {
		timer.tick(); console::dump( "Adding sources..." );
		//
		// Velocity
		grid.iterate_active_faces([&](const face_id3 &face_id, int tid) {
			const vec3d p = grid.get_face_position(face_id);
			double d(0.0); vec3d u;
			add_func (p,u,d,time,dt);
			// density.increment(i,j,k,d);
			grid.velocity[face_id.index] = u[face_id.dim];
		});
		//
		// Density
		double min_x, min_y, min_z = 100.0;
		double max_x, max_y, max_z = -100.0;
		auto add_density = [&]() {
			grid.iterate_active_cells([&](const cell_id3 &cell_id, int tid) {
				const vec3d p = grid.get_cell_position(cell_id);
				min_x = std::min(min_x,p[0]); min_y = std::min(min_y,p[1]); min_z = std::min(min_z,p[2]);
				max_x = std::max(max_x,p[0]); max_y = std::max(max_y,p[1]); max_z = std::max(max_z,p[2]);
				//console::dump( "Adding source at (%f,%f,%f)\n", p[0], p[1], p[2] );
				//こっから
				vec3d center (0.15, 0.15, 0.5);
				if ((p-center).len() < 0.1) {
				double d(0.0); vec3d dummy;       
				add_func (p,dummy,d,time,dt);
				// density.increment(i,j,k,d);
				grid.density[cell_id.index] += d;
				grid.qc[cell_id.index] += d; // Example: 90% of added density goes to water vapor
				//console::dump("add_source");
				//if (d != 0.0) console::dump( "add_source density[%d] = %f\n", cell_id.index, m_grid->density[cell_id.index] );
				}

			});
		};
		//
		// Density
		unsigned seeded (0);
		// if( m_param.use_dust ) {
		// 	//
		// 	std::random_device rd;
		// 	std::mt19937 gen(rd());
		// 	std::uniform_real_distribution<> dis(-1.0,1.0);
		// 	//
		// 	add_density(m_accumulation);
		// 	//
		// 	double scale = 1.0 / pow(m_param.r_sample,DIM3);
		// 	bool should_re_rasterize (false);
		// 	m_accumulation.serial_op([&]( int i, int j, int k, auto &it) {
		// 		double d = it();
		// 		while( d > scale ) {
		// 			vec3d p = m_dx*vec3i(i,j,k).cell()+0.5*m_dx*vec3d(dis(gen),dis(gen),dis(gen));
		// 			m_dust_particles.push_back(p); ++ seeded;
		// 			should_re_rasterize = true;
		// 			d -= scale;
		// 		}
		// 		it.set(d);
		// 	});
		// 	//
		// 	if( should_re_rasterize ) {
		// 		rasterize_dust_particles(density);
		// 	}
			//
		//} else {
			add_density();
			//console::dump( "Source added in region x:(%f,%f) y:(%f,%f) z:(%f,%f)\n", min_x, max_x, min_y, max_y, min_z, max_z );
		//}
		//
		if( m_param.use_dust ) console::dump( "Done. Seeded=%d. Took %s.\n", seeded, timer.stock("add_func").c_str());
		else console::dump( "Done. Took %s.\n", timer.stock("add_func").c_str());
	}
}
//
// void cloud_oc::rasterize_dust_particles( array3<Real> &rasterized_density ) {
// 	//
// 	rasterized_density.clear();
// 	double scale = 1.0 / pow(m_param.r_sample,DIM3);
// 	for( const vec3d &p : m_dust_particles ) {
// 		vec3i pi = p/m_dx;
// 		if( ! m_shape.out_of_bounds(pi)) {
// 			rasterized_density.increment(pi[0],pi[1],pi[2],scale);
// 		}
// 	}
// }
//
// void cloud_oc::add_buoyancy_force( macarray3<Real> &velocity, const array3<Real> &density, double dt ) {
// 	//
// 	velocity[1].parallel_all([&]( int i, int j, int k, auto &it, int tn ) {
// 		vec3d pi = vec3i(i,j,k).face(1);
// 		Real d = array_interpolator3::interpolate<Real>(density,(pi-vec3d(0.5,0.5,0.5)));
// 		it.increment(m_param.buoyancy_factor*dt*d);
// 	});
// }
//
void cloud_oc::idle() {
	//
	//write_to_txt("before_timestep"); //added
	scoped_timer timer(this);
	//
	// Add to graph
	//add_to_graph();
	//
	// Compute the timestep size
	timer.tick(); console::dump( "Computing time step...");
	double max_u_per_unit (0.0);
	//get_finest_dxが0しか返していないのでdtがほぼ0になってしまう→たぶんm_gridまわりでlayerを追加する処理をしてないから
	for( size_t n=0; n<m_grid->velocity.size(); ++n ) max_u_per_unit = std::max(max_u_per_unit,(double)std::abs(m_grid->velocity[n]));
	const double dt = m_timestepper->advance(max_u_per_unit,m_grid->get_finest_dx());
	const double time = m_timestepper->get_current_time();
	const double CFL = m_timestepper->get_current_CFL();
	//console::dump( "Done. time=%.2e. max_u_perunit=%.2e, dx=%2.e, dt=%.2e,CFL=%.2f. Took %s\n", time, max_u_per_unit, m_grid->get_finest_dx(), dt, CFL, timer.stock("compute_timestep").c_str());
	//timer.tick(); console::dump( ">>> %s step started (dt=%.2e,CFL=%.2f)...\n", CFL, console::nth(step).c_str());
	//
    m_accumulated_CFL += CFL;
    
	std::swap(m_grid,m_grid_prev);
    if( m_accumulated_CFL >= m_param.maximal_CFL_accumulation ) {
        m_accumulated_CFL = 0.0;

		timer.tick(); console::dump( ">>> Remeshing...\n");
		if( m_param.use_sizing_func ) {
			// if( m_do_inject ) {
			// 	m_macoctreesizingfunc.activate_cells(*m_grid_prev,*m_grid,m_combined_solid_func,
			// 		[&]( const vec3d &p ) {
			// 			vec3d u; double value;
			// 			m_inject_func(p,m_dx,dt,time,step,value,u);
			// 			return value;
			// 		});
			// } else {
				m_macoctreesizingfunc.activate_cells(*m_grid_prev,*m_grid,m_combined_solid_func,nullptr);
			//}
		} else {
			m_grid->activate_cells([&]( const vec3d &p ){
				//
				double inject_levelset (std::numeric_limits<double>::max());
				// if( m_do_inject ) {
				// 	vec3d u; double value;
				// 	if( m_inject_func(p,m_dx,dt,time,step,value,u)) {
				// 		inject_levelset = value;
				// 	}
				// }
				vec3d u (m_grid_prev->sample_velocity(p));
				return std::min(inject_levelset,m_grid_prev->sample_levelset(p-dt*u));
			},m_combined_solid_func);
		}

        m_grid->balance_layers();
        m_grid->assign_indices();
		console::dump( "<<< Done. Took %s.\n", timer.stock("remeshing").c_str());

    } else {
		console::dump( "Accumulated CFL=%.2f. Copying previous grid...", m_accumulated_CFL );
		timer.tick();
		m_grid->copy(*m_grid_prev);
		console::dump( "Done. Took %s.\n", timer.stock("copy_grid").c_str());
    }
	
	//added
	// Advect density
	m_grid->assign_density([&]( const vec3d &p ) {
		vec3d u (m_grid_prev->sample_velocity(p));
		vec3d center (0.15,0.15,0.5);
		vec3d one (1.0, 1.0, 1.0);
		//if (m_grid_prev->sample_density(-p) > 0.000001 && (p-center).len() >= 0.1) console::dump( "density at (%f,%f,%f): (%f) advected by (%f,%f,%f), dt = %f, len = %f\n", p[0], p[1], p[2], m_grid->sample_density(p-dt*u), u[0], u[1], u[2], dt, (p-center).len() );
		return m_grid_prev->sample_density(p-dt*u);
	});
	//write_to_txt("after_density_advection"); //added
	// Update solid
	// m_macutility->update_solid_variables(m_dylib,time,&m_solid,&m_solid_velocity);

	//added
	if( m_param.maccormack ) {
		//
		timer.tick(); console::dump( ">>> MacCormack velocity advection...\n");
		//
		using Real2 = struct { Real v[2] = {0.0, 0.0}; };
		std::vector<Real2> min_max_values(m_grid->face_count);
		std::vector<Real> u0(m_grid->face_count), u1(m_grid->face_count), _u0(m_grid->face_count);
		std::vector<Real> u_x(m_grid->face_count), u_y(m_grid->face_count), u_z(m_grid->face_count);
		std::vector<char> near_surface_flag (m_grid->face_count);
		//
		timer.tick(); console::dump( "Mapping initial velocity...");
		m_grid->iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
			const vec3d &p = m_grid->get_face_position(face_id);
			uint_type index = face_id.index;
			u_x[index] = m_grid_prev->sample_velocity(p,0);
			u_y[index] = m_grid_prev->sample_velocity(p,1);
			u_z[index] = m_grid_prev->sample_velocity(p,2);
			u0[index] = m_grid_prev->sample_velocity(p,face_id.dim);
		});
		console::dump( "Done. Took %s.\n", timer.stock("maccormack_initial_mapping").c_str());
		//
		timer.tick(); console::dump( "Forward advection...");
		m_grid->iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
			const vec3d &p = m_grid->get_face_position(face_id);
			uint_type index = face_id.index;
			vec3d u (u_x[index],u_y[index],u_z[index]);
			u1[index] = m_grid_prev->sample_face(p-dt*u,face_id.dim,m_grid_prev->velocity,min_max_values[index].v);
			const double dx = m_grid->get_face_dx(face_id);
			near_surface_flag[index] = m_grid_prev->sample_levelset(p-dt*u) > -dx || m_grid_prev->sample_levelset(p+dt*u) > -dx;
		});
		console::dump( "Done. Took %s.\n", timer.stock("maccormack_forward_advection").c_str());
		//
		timer.tick(); console::dump( "Backward advection...");
		m_grid->iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
			uint_type index = face_id.index;
			if( ! near_surface_flag[index] ) {
				const vec3d &p = m_grid->get_face_position(face_id);
				vec3d u (u_x[index],u_y[index],u_z[index]);
				_u0[index] = m_grid->sample_face(p+dt*u,face_id.dim,u1);
			}
		});
		console::dump( "Done. Took %s.\n", timer.stock("maccormack_backward_advection").c_str());
		//
		timer.tick(); console::dump( "Combining the values...");
		std::vector<size_t> reverted_count (m_grid->parallel.get_thread_num());
		std::vector<size_t> filled_count (m_grid->parallel.get_thread_num());
		m_grid->iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
			const uint_type index = face_id.index;
			const Real *v = min_max_values[index].v;
			bool in_the_air (true);
			m_grid->get_gradient(face_id,[&]( const cell_id3 &cell_id, double value, const grid3::gradient_info3 &info ) {
				if( m_grid->levelset[cell_id.index] < 0.0 ) {
					in_the_air = false;
				}
			});
			if( ! in_the_air ) filled_count[tid] ++;
			if( near_surface_flag[index] ) {
				if( ! in_the_air ) reverted_count[tid] ++;
				m_grid->velocity[index] = u1[index];
			} else {
				m_grid->velocity[index] = std::max(v[0],std::min(v[1],u1[index]+0.5f*(u0[index]-_u0[index])));
			}
		});
		size_t num_inside = std::accumulate(filled_count.begin(),filled_count.end(),0);
		size_t num_reverted = std::accumulate(reverted_count.begin(),reverted_count.end(),0);
		double reverted_ratio = num_reverted / (double) num_inside;
		console::dump( "Done. Reverted=%u/%u (%.2f%%). Took %s.\n", num_reverted, num_inside, 100.0*reverted_ratio, timer.stock("maccormack_combine").c_str());
		console::dump( "<<< Done. Took %s.\n", timer.stock("maccormack_velocity_advection").c_str());
	} else {
		m_grid->set_velocity([&]( const vec3d &p, char dim ) {
			vec3d u (m_grid_prev->sample_velocity(p));
			return m_grid_prev->sample_velocity(p-dt*u,dim);
		});
	}

	m_macoctreeproject.assemble_matrix(*m_grid);
	auto solid_velocity_func = [&]( const vec3d &p ) {
		return m_moving_solid_func ? m_moving_solid_func(time,p).second : vec3d();
	};
	// Project(added)
	m_macoctreeproject.project_density(*m_grid,dt,solid_velocity_func);
	// //added
	// m_grid->extrapolate_toward_solid(m_combined_solid_func);
	// m_grid->extrapolate(m_combined_solid_func);
	//added
	//
	// shared_macarray3<Real> velocity_save(m_velocity);
	// m_macadvection->advect_vector(m_velocity,velocity_save(),m_fluid,dt,"velocity");
	//
	// Add buoyancy force
	//add_buoyancy_force(m_velocity,m_density,dt);
	//
	// Add source
	//add_source(m_velocity,m_density,m_timestepper->get_current_time(),dt);
	
	add_source_oc(m_timestepper->get_current_time(),dt, *m_grid);
	//
	m_grid->vorticity_confinement(dt);
	m_grid->add_buoyancy(dt);
	//microphysics_cloud(*m_grid, dt);
	// Add external force
	//inject_external_force(m_velocity);
	//
	// Projection
	//m_macproject->project(dt,m_velocity,m_solid,m_solid_velocity,m_fluid,0.0);
	//m_macutility->extrapolate_and_constrain_velocity(m_solid,m_velocity,m_param.extrapolated_width);
	//
	//console::dump( "<<< %s step done. Took %s\n", console::nth(step).c_str(), timer.stock("simstep").c_str());
	//
	// Export density
	//export_density();
	//
	// Report stats
	//m_macstats->dump_stats(m_solid,m_fluid,m_velocity,m_timestepper.get());
}
//
// void cloud_oc::advect_dust_particles( const macarray3<Real> &velocity, double dt ) {
// 	//
// 	m_parallel.for_each( m_dust_particles.size(), [&]( size_t n, int tn ) {
// 		vec3d &p = m_dust_particles[n];
// 		vec3d u0 = macarray_interpolator3::interpolate<Real>(velocity,p/m_dx);
// 		vec3d u1 =  macarray_interpolator3::interpolate<Real>(velocity,(p+dt*u0)/m_dx);
// 		p += 0.5 * dt * (u0+u1);
// 	});
// 	//
// 	m_parallel.for_each( m_dust_particles.size(), [&]( size_t n, int tn ) {
// 		vec3d &p = m_dust_particles[n];
// 		double phi = array_interpolator3::interpolate<Real>(m_solid,p/m_dx);
// 		if( phi < 0.0 ) {
// 			Real derivative[DIM3];
// 			array_derivative3::derivative(m_solid,p/m_dx,derivative);
// 			p = p - phi*vec3d(derivative).normal();
// 		}
// 		for( unsigned dim : DIMS3 ) {
// 			if( p[dim] < 0.0 ) p[dim] = 0.0;
// 			if( p[dim] > m_dx*m_shape[dim] ) p[dim] = m_dx*m_shape[dim];
// 		}
// 	});
// 	//
// 	rasterize_dust_particles(m_density);
// }
//
// void cloud_oc::add_to_graph() {
// 	//
// 	if( m_param.show_graph ) {
// 		//
// 		// Compute total energy
// 		const double time = m_timestepper->get_current_time();
// 		const double total_energy = m_macutility->get_kinetic_energy(m_solid,m_fluid,m_velocity);
// 		//
// 		// Add to graph
// 		m_graphplotter->add_point(m_graph_id,time,total_energy);
// 	}
// }
//
// void cloud_oc::draw_dust_particles( graphics_engine &g ) const {
// 	using ge = graphics_engine;
// 	g.color4(1.0,1.0,1.0,1.0);
// 	g.begin(ge::MODE::POINTS);
// 	for( const vec3d &p : m_dust_particles ) {
// 		g.vertex3v(p.v);
// 	}
// 	g.end();
// }
//
void cloud_oc::draw( graphics_engine &g ) const {
	//
	const double time = m_timestepper->get_current_time();
	//
	// Visualize solid
	//shared_array3<Real> solid_to_visualize(m_solid.shape());
	//m_gridutility->assign_visualizable_solid(m_dylib,m_dx,solid_to_visualize());
	//m_gridvisualizer->draw_solid(g,solid_to_visualize());
	//
	// Visualize moving solid
	auto draw_func = reinterpret_cast<void(*)(graphics_engine &,double)>(m_dylib.load_symbol("draw"));
	if( draw_func ) {
		g.color4(1.0,0.8,0.5,0.3);
		draw_func(g,m_timestepper->get_current_time());
	}
	//
	// Draw velocity
	//m_macvisualizer->draw_velocity(g,m_velocity);
	//
	// Draw projection component
	//m_macproject->draw(g);
	//
	// Draw concentration
	// if( m_param.use_dust ) draw_dust_particles(g);
	// else m_gridvisualizer->draw_density(g,m_density);
	m_gridvisualizer_oc->draw_qc(g,*m_grid);
	//write_to_txt("before_draw"); //added
	//m_gridvisualizer->draw_density_oc(g, *m_grid);
	// //densityを点で描画
	// m_grid->draw_density_oc(g);
    //         g.point_size(3.0);
    //         g.begin(graphics_engine::MODE::POINTS);
            
    //         //grid.densityを使用して密度を描画
    //         m_grid->iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
    //             double density_value = m_grid->density[cell_id.index];
    //             vec3d pos = m_grid->get_cell_position(cell_id);
    //             g.color4(1.0,1.0,1.0,density_value);
    //             g.vertex3v(pos.v);
    //         });
	// Draw graph
	m_graphplotter->draw(g);
}
//
void cloud_oc::export_density () const {
	//
	scoped_timer timer(this);
	if( console::get_root_path().size()) {
		int frame = m_timestepper->should_export_frame();
		if( frame ) {
			timer.tick(); console::dump( "Exporting %s density...", console::nth(frame).c_str());
			do_export_density(frame);
			console::dump( "Done. Took %s\n", timer.stock("export_mesh").c_str());
			if( m_param.render_density ) {
				render_density(frame);
			}
		}
	}
}
//
void cloud_oc::do_export_density( int frame ) const {
	//
	std::string dir_path = console::get_root_path()+"/density";
	if( ! filesystem::is_exist(dir_path)) filesystem::create_directory(dir_path);
	//
	std::string path = console::format_str("%s/%d_density.vol",dir_path.c_str(),frame);
	FILE *fp = fopen(path.c_str(),"wb");
	assert(fp);
	//
	const char *vol_str = "VOL";
	const char version (3);
	const int value (1), xn (m_shape[0]), yn (m_shape[1]), zn (m_shape[2]);
	//
	const Real minX (-0.5*xn*m_dx);
	const Real minY (-0.5*yn*m_dx);
	const Real minZ (-0.5*zn*m_dx);
	const Real maxX (0.5*xn*m_dx);
	const Real maxY (0.5*yn*m_dx);
	const Real maxZ (0.5*zn*m_dx);
	//
	fwrite(vol_str,3,1,fp);
	fwrite(&version,sizeof(version),1,fp);
	fwrite(&value,sizeof(value),1,fp);
	fwrite(&xn,sizeof(xn),1,fp);
	fwrite(&yn,sizeof(yn),1,fp);
	fwrite(&zn,sizeof(zn),1,fp);
	fwrite(&value,sizeof(value),1,fp);
	fwrite(&minX,sizeof(minX),1,fp);
	fwrite(&minY,sizeof(minY),1,fp);
	fwrite(&minZ,sizeof(minZ),1,fp);
	fwrite(&maxX,sizeof(maxX),1,fp);
	fwrite(&maxY,sizeof(maxY),1,fp);
	fwrite(&maxZ,sizeof(maxZ),1,fp);
	//
	std::vector<Real> density_linearized(xn*yn*zn);
	// m_density.const_parallel_all([&](int i, int j, int k, auto &it) {
	// 	density_linearized[i+j*(xn)+k*(xn*yn)] = it();
	// });
	for (size_t n=0; n<xn*yn*zn; ++n) {
		fwrite(&density_linearized[n],sizeof(value),1,fp);
	}
	//
	fclose(fp);
}
//
void cloud_oc::render_density( int frame ) const {
	//
	scoped_timer timer(this);
	global_timer::pause();
	//
	assert(console::get_root_path().size());
	//
	std::string mitsuba_path = console::get_root_path() + "/smoke_mitsuba";
	std::string copy_from_path = filesystem::find_resource_path("smoke","mitsuba");
	if( ! filesystem::is_exist(mitsuba_path)) {
		if( filesystem::is_exist(copy_from_path)) {
			console::run( "cp -r %s %s", copy_from_path.c_str(), mitsuba_path.c_str());
		} else {
			console::dump( "Could not lcoate mitsuba files (%s).\n", copy_from_path.c_str());
			exit(0);
		}
	}
	//
	std::string render_command = console::format_str("cd %s; /usr/bin/python render.py %d %d %g %s",
				mitsuba_path.c_str(),frame,m_param.render_sample_count,m_param.volume_scale,"img");
	//
	console::dump("Running command: %s\n", render_command.c_str());
	console::system(render_command.c_str());
	//
	global_timer::resume();
}
//
// void cloud_oc::do_inject_external_density( double dt, double time, unsigned step ) {
// 	//
// 	if( m_do_inject ) {
// 		//
// 		std::vector<double> total_injected (m_grid->parallel.get_thread_num());
// 		m_grid->iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
// 			//
// 			const vec3d p = m_grid->get_cell_position(cell_id);
// 			double value (0.0); vec3d u;
// 			//
// 			if( m_inject_func(p,m_dx,dt,time,step,value,u)) {
// 				const double v = m_grid->density[cell_id.index];
// 				m_grid->density[cell_id.index] = std::min(value,v);
// 				if( value < 0.0 && v > 0.0 ) {
// 					const double dx = m_grid->get_cell_dx(cell_id);
// 					total_injected[tid] += dx*dx*dx;
// 				}
// 			}
// 		});
// 		m_grid->iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
// 			const vec3d p = m_grid->get_face_position(face_id);
// 			const double dx = m_grid->get_face_dx(face_id);
// 			double value (0.0); vec3d u;
// 			if( m_inject_func(p,m_dx,dt,time,step,value,u)) {
// 				if( value < dx ) {
// 					m_grid->velocity[face_id.index] = u[face_id.dim];
// 				}
// 			}
// 		});
// 		m_injected_volume = std::accumulate(total_injected.begin(),total_injected.end(),0.0);
// 	}
// }
//
extern "C" module * create_instance() {
	return new cloud_oc;
}
//
extern "C" const char *license() {
	return "MIT";
}
//
