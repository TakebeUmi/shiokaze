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
#include <numeric>
#include <vector>

#include "../../../PerlinNoise/PerlinNoise.hpp"

//
SHKZ_USING_NAMESPACE
//

cloud_oc::cloud_oc () {
	//
	int coeff = 1;
	m_shape = shape3(64*coeff,64*coeff,64*coeff);
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
double cloud_oc::circle_source(const vec3d &p) const {
	vec2d center (0.5,0.5);
	double radius(0.2);
	double ratio(0.8);
	vec2d p_2d (p[0], p[2]);
	if ((p_2d - center*m_param.scale).len() < radius*ratio*m_param.scale ) {
		return 3.0;
	} else if ((p_2d - center*m_param.scale).len() < radius*m_param.scale ) {
		return -12.5*p_2d.len()/m_param.scale + 13.0;
	}
	else return 0.5;
}

double cloud_oc::perlin_noise_source(const vec3d &p) const {
	const siv::PerlinNoise::seed_type seed = 123456u;
	static siv::PerlinNoise perlin{seed}; // Seed for reproducibility
	vec2d p_2d(p[0], p[2]);
	const double noise = perlin.octave2D_01(p_2d[0], p_2d[1], 3, 0.3);
	return noise;
}

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

void cloud_oc::print_position(grid3 &grid) const {
	grid.iterate_active_cells([&]( const cell_id3 &cell_id, int thread_index ) {
		vec3d pos = grid.get_cell_position(cell_id);
		if (grid.get_cell_position(cell_id)[0] > 11000){
			console::dump("Cell Index: (%d, %d, %d) Position: (%.6f, %.6f, %.6f)\n",
				cell_id.pi[0], cell_id.pi[1], cell_id.pi[2],
				pos[0], pos[1], pos[2]);
		}
	});
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
	view_scale *= m_param.scale;
	//
	double resolution_scale (1.0);
	config.get_double("ResolutionScale",resolution_scale,"Resolution doubling scale");
	//
	//resolution_scale = 10.0;
	m_shape *= resolution_scale;
	m_dx = view_scale * m_shape.dx();
	console::dump("Computed dx: %e\n",m_dx);
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
			m_grid->assign_qc([&](const vec3d &p) {
				return 0.0;
			});
			m_grid->assign_qv([&](const vec3d &p) {
				return 0.0;
			});
			m_grid->assign_qr([&](const vec3d &p) {
				return 0.0;
			});
			m_grid->atmospheric_temperature();
			// m_grid->assign_theta([&](const vec3d &p) {
			// 	return 300.0;
			// });
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
	m_camera->look_at(m_param.target.v,m_param.origin.v,vec3d(0.0,-1.0,0.0).v,10000.0);
	//look_at関数でup, ターゲット位置、カメラ原点位置、視野（Field of View:FOV）を設定
	//https://ryichando.graphics/shiokaze/classcamera3__interface.html#a4bdebf6662eef4ab7295686768ebe381
	console::dump( "<<< Initialization finished. Took %s\n", timer.stock("initialization").c_str());
	//
	if( m_param.show_graph ) {
		m_graphplotter->clear();
		m_graph_id = m_graphplotter->create_entry("Kinetic Energy");
	}
	//print_position(*m_grid);
	//m_grid->check_buoyancy();
	std::vector<cell_id3> cell_ids(m_grid->cell_count);
	m_grid->iterate_active_cells([&]( const cell_id3 &cell_id, int thread_index ) {
		cell_ids[cell_id.index] = cell_id;
	});
	console::dump("Total active cells after initialization: %zu\n", cell_ids.size());
	uf.initialize(m_grid->cell_count, cell_ids);

}

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
		if (z[1]>=grid.param.scale - m_dx) newqv = 0.0;
		if (z[1]>=grid.param.scale - m_dx) newqc = 0.0;
		if (z[1]>=grid.param.scale - m_dx) newqr = 0.0;
		grid.qv[cell_id.index] = newqv;
		grid.qc[cell_id.index] = newqc;
		grid.qr[cell_id.index] = newqr;
		grid.theta[cell_id.index] = newtheta;
	});
}


void cloud_oc::add_source_cloud (double time, double dt , grid3 &grid) {
	//
	scoped_timer timer(this);
	//
	// auto add_func = reinterpret_cast<void(*)(const vec3d &, vec3d &, double &, double, double)>(m_dylib.load_symbol("add"));
	int octaves = 3;
	double amplitude = 1.0;
	double persistence = 0.3;
	auto source_func = [&](const vec3d &p, double time, double dt, double &theta, double &vapor) {
		vec2d p_2d (p[0], p[2]);
		vec2d center (0.5,0.5);
		if (/*(p_2d - center*m_param.scale).len() < 0.2*m_param.scale && */p[1] <= m_dx) {
				vec3d offset(1.0, 1.0, 1.0);
				double T0 = grid.param.T0;
				double p0 = grid.param.p0;
				double Gamma = grid.param.gamma;
				double z1 = grid.param.z1;
				double g = grid.param.g;
				double TISA = grid.atmospheric_temperature(T0, Gamma, z1, p[1]);
				double pISA = grid.atmospheric_pressure(T0, p0, Gamma, g, p[1]);
				double theta_val = grid.thermal_potential_temperature(TISA, p0, pISA, 0.0);
		// console::dump("theta[%f,%f,%f] = %f\n", xyz[0], xyz[1], xyz[2], theta_val);

				theta = theta_val +  perlin_noise_source(p)*circle_source(p)/m_param.source_coeff;
				vapor =  perlin_noise_source(p+offset)/m_param.source_coeff;
				// console::dump("Source added at (%.6f, %.6f, %.6f): vapor=%.6f, theta=%.6f\n",
				// 	p[0], p[1], p[2],
				// 	vapor, theta);
		}

	};
	
		timer.tick(); console::dump( "Adding sources..." );
		//
		// Velocity
		// double min_theta(1000.0); double min_vapor(1000.0);
		grid.iterate_active_cells([&](const cell_id3 &cell_id, int tid) {
			const vec3d p = grid.get_cell_position(cell_id);
			double theta(grid.theta[cell_id.index]); double vapor(grid.qv[cell_id.index]);
			// if (cell_id.index == 163455) console::dump("p:(%.6f, %.6f, %.6f), cell_id.index=%d\n",
			// 	p[0], p[1], p[2], cell_id.index);

			source_func (p,time,dt,theta,vapor);
			// min_theta = std::min(min_theta,theta);
			// min_vapor = std::min(min_vapor,vapor);
			// density.increment(i,j,k,d);
			// if (vapor > 0.0) {
			// 	console::dump("Source added at (%d, %d, %d): vapor=%.6f, theta=%.6f\n",
			// 		cell_id.pi[0], cell_id.pi[1], cell_id.pi[2],
			// 		vapor, theta);
			// }
			grid.qv[cell_id.index] = vapor;
			grid.theta[cell_id.index] = theta;
		});
		// console::dump( "Done. Min theta=%.6f, Min vapor=%.6f.\n", min_theta, min_vapor);
}
//
void cloud_oc::idle() {
	//
	console::dump( "=================== Step %d ===================\n", m_timestepper->get_step_count() );
	//write_to_txt("before_timestep"); //added

	// m_grid->serial_iterate_active_cells([&]( const cell_id3 &cell_id ) {
	// 	if (m_grid->get_cell_position(cell_id)[1] == 93.75 && 
	// 		m_grid->get_cell_position(cell_id)[2] == 93.75) console::dump("cell_id.index: %d p: (%f, %f, %f)\n", cell_id.index, m_grid->get_cell_position(cell_id)[0], m_grid->get_cell_position(cell_id)[1], m_grid->get_cell_position(cell_id)[2]);
	// 	// if (cell_id.index == 163455 ) console::dump("Before timestep - p: (%.6f, %.6f, %.6f)\n", p[0], p[1], p[2]);
	// });
	scoped_timer timer(this);
	m_grid->iterate_active_cells([&]( const cell_id3 &cell_id, int thread_index ) {
		vec3d p = m_grid->get_cell_position(cell_id);
	if (cell_id.index < 0 ) console::dump("p: (%.6f, %.6f, %.6f)\n", p[0], p[1], p[2]);
	});
	//
	// Add to graph
	//add_to_graph();
	//
	// Compute the timestep size
	timer.tick(); console::dump( "Computing time step...");
	double max_u_per_unit (0.0);
	//get_finest_dxが0しか返していないのでdtがほぼ0になってしまう→たぶんm_gridまわりでlayerを追加する処理をしてないから
	for( size_t n=0; n<m_grid->velocity.size(); ++n ) max_u_per_unit = std::max(max_u_per_unit,(double)std::abs(m_grid->velocity[n]));
	m_param.Time_coeff * m_timestepper->advance(max_u_per_unit,m_grid->get_finest_dx());
	const double dt = 5.0;
	const double time = m_timestepper->get_current_time();
	const double CFL = m_timestepper->get_current_CFL();
	console::dump( "Done. time=%.2e. max_u_perunit=%.2e, dx=%2.e, dt=%.2e,CFL=%.2f. Took %s\n", time, max_u_per_unit, m_grid->get_finest_dx(), dt, CFL, timer.stock("compute_timestep").c_str());
	//timer.tick(); console::dump( ">>> %s step started (dt=%.2e,CFL=%.2f)...\n", CFL, console::nth(step).c_str());
	//
    m_accumulated_CFL += CFL;
	m_grid_prev->copy(*m_grid);
	//advect velocity
	m_grid->set_velocity([&]( const vec3d &p, char dim ) {
		vec3d u (m_grid_prev->sample_velocity(p));
		// if (u[1] <= 1.0) u[1] = 1.0; // prevent going down too fast
		return m_grid_prev->sample_velocity(p-dt*u,dim);
	});
	m_grid->vorticity_confinement(dt);
	m_grid->source_func(m_dx);
	m_grid->add_buoyancy(dt);
	//added
	m_grid->serial_iterate_active_cells([&]( const cell_id3 &cell_id ) {
		m_grid->right_cell_neighbor(cell_id, [&]( char dim, const cell_id3 &neighbor_id ) {
			if ((m_grid->qc[cell_id.index] == 0.0 && m_grid->qc[neighbor_id.index] == 0.0) || (m_grid->qc[cell_id.index] > 0.0 && m_grid->qc[neighbor_id.index] > 0.0 )) {
				uf.unite(*m_grid, cell_id.index, neighbor_id.index);
			}
		});
	});
	std::vector<cell_id3_and_is_merged> cell_ids_included_merged_cells;
	int merged_cell_count = 0;

	for (int n = 0; n < uf.par.size(); n++) {
		if (uf.par[n] == -1) {
			cell_ids_included_merged_cells.push_back({uf.cell_ids[n], uf.top[n].index!=uf.bottom[n].index});
			if (uf.top[n].index!=uf.bottom[n].index) {
				merged_cell_count++;
				//console::dump("Merged cell at index %d\n", n);
			}
			//cell_idsのindexを保存。この値をcell_ids[i]にアクセスすることでcell_idを得られる。その際、is_mergedがtrueならば、merged cellである。
		}
	}
	int all_cell_count = cell_ids_included_merged_cells.size();

	m_macoctreeproject.assemble_matrix_qc(*m_grid);
	auto solid_velocity_func = [&]( const vec3d &p ) {
		return m_moving_solid_func ? m_moving_solid_func(time,p).second : vec3d();
	};
	// // // Project(added)
	m_macoctreeproject.project_cloud(*m_grid,dt,solid_velocity_func);
	//m_grid->check_buoyancy();
	std::swap(m_grid,m_grid_prev);
    if( m_accumulated_CFL >= m_param.maximal_CFL_accumulation ) {
        m_accumulated_CFL = 0.0;

		timer.tick(); console::dump( ">>> Remeshing...\n");
		if( m_param.use_sizing_func ) {
				m_macoctreesizingfunc.activate_cells(*m_grid_prev,*m_grid,m_combined_solid_func,nullptr);
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
	


	// Advect scalar quantities
	double dt_advect = 5.0;
	m_grid->assign_qv([&]( const vec3d &p ) {
		vec3d u (m_grid_prev->sample_velocity(p));
		vec3d advected_point = p - dt_advect * u;
		if (u[1] <= 1.0) u[1] = 1.0; 
		//console::dump("Advecting at position (%.6f, %.6f, %.6f)\n", u[0], u[1], u[2]);
		// vec3d u (50.0,50.0,50.0); // test
		double d =  m_grid_prev->sample_qv(p - dt_advect * u);
		if (d > 1.0) return 0.0;
		else return d;
	});
	m_grid->assign_qc([&]( const vec3d &p ) {
		vec3d u (m_grid_prev->sample_velocity(p));
		if (u[1] <= 1.0) u[1] = 1.0; // prevent going down too fast
		double d =  m_grid_prev->sample_qc(p - dt_advect * u);
		if (d > 1.0) return 0.0;
		else return d;
	});
	m_grid->assign_qr([&]( const vec3d &p ) {
		vec3d u (m_grid_prev->sample_velocity(p));
		if (u[1] <= 1.0) u[1] = 1.0; // prevent going down too fast
		double d =  m_grid_prev->sample_qr(p - dt_advect * u);
		if (d > 1.0) return 0.0;
		else return d;
	});
	m_grid->assign_theta([&]( const vec3d &p ) {
		vec3d u (m_grid_prev->sample_velocity(p));
		double d =  m_grid_prev->sample_theta(p - dt_advect * u);
		return d;
	});
	//write_to_txt("after_density_advection"); //added

	//added
	




	// //added
	// m_grid->extrapolate_toward_solid(m_combined_solid_func);
	// m_grid->extrapolate(m_combined_solid_func);
	//added
	
	//add_source_cloud(m_timestepper->get_current_time(),dt, *m_grid);
	
	//
	
	//m_grid->check_buoyancy();
	microphysics_cloud(*m_grid, dt);
	// Add external force
	//inject_external_force(m_velocity);
	//
	// Projection
	//m_grid->check_buoyancy();
	//print_position(*m_grid);
}
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
	
	// Draw concentration
	m_gridvisualizer_oc->draw_qc(g,*m_grid);
	// m_gridvisualizer_oc->draw_qv(g,*m_grid);
	// m_gridvisualizer_oc->draw_qr(g,*m_grid);
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
extern "C" module * create_instance() {
	return new cloud_oc;
}
//
extern "C" const char *license() {
	return "MIT";
}
//
