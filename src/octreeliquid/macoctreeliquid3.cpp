/*
**	macoctreeliquid3.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on Februrary 8, 2019. All rights reserved.
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
#include "macoctreeliquid3.h"
#include "macoctreemesher3.h"
#include <shiokaze/array/array_interpolator3.h>
#include <shiokaze/array/array_derivative3.h>
#include <shiokaze/array/shared_array3.h>
#include <shiokaze/array/bitmacarray3.h>
#include <shiokaze/utility/utility.h>
#include <shiokaze/array/array_utility3.h>
#include <shiokaze/array/array_extrapolator3.h>
#include <shiokaze/array/macarray_interpolator3.h>
#include <shiokaze/array/macarray_extrapolator3.h>
#include <shiokaze/array/shared_array_core3.h>
#include <shiokaze/core/filesystem.h>
#include <shiokaze/core/console.h>
#include <shiokaze/core/timer.h>
#include <cmath>
#include <numeric>
//
SHKZ_USING_NAMESPACE
using namespace macotreeliquid3_namespace;
//
void macoctreeliquid3::load( configuration &config ) {
	//
	std::string name("waterdrop3"); config.get_string("Name",name,"Scene file name");
	m_dylib.open_library(filesystem::resolve_libname(name));
	m_dylib.load(config);
	m_dylib.overwrite(config);
}
//
void macoctreeliquid3::configure( configuration &config ) {
	//
	if( console::get_root_path().size()) {
		m_export_path = console::get_root_path() + "/mesh";
		if( ! filesystem::is_exist(m_export_path)) {
			filesystem::create_directory(m_export_path);
		}
	}
	//
	m_dylib.configure(config);
	//
	m_param.render_mesh = console::system("mitsuba > /dev/null 2>&1") == 0;
	//
	m_shape = shape3(64,32,64);
	//
	config.set_default_string("LinSolver","amg");
	//
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
	//
	double view_scale (1.0);
	config.get_double("ViewScale",view_scale,"View scale");
	//
	double resolution_scale (1.0);
	config.get_double("ResolutionScale",resolution_scale,"Resolution doubling scale");
	//
	m_shape *= resolution_scale;
	m_dx = view_scale * m_shape.dx();
	//
	unsigned solid_max_resolution (512);
	config.get_unsigned("SolidMaxResolution",solid_max_resolution,"Max resolution for solid level set for visualization");
	m_solid_shape = m_shape;
	m_solid_dx = m_dx;
	while( m_solid_shape.max() > solid_max_resolution ) {
		m_solid_shape = m_solid_shape/2;
		m_solid_dx = 2.0 * m_solid_dx;
	}
	m_solid_gridvisualizer->set_environment("shape",&m_solid_shape);
	m_solid_gridvisualizer->set_environment("dx",&m_solid_dx);
	//
	m_solid_mesher->set_environment("shape",&m_solid_shape);
	m_solid_mesher->set_environment("dx",&m_solid_dx);
}
//
void macoctreeliquid3::setup_window( std::string &name, int &width, int &height ) const {
	height = width;
}
//
bool macoctreeliquid3::should_quit() const {
	return m_should_quit_on_save || m_timestepper->should_quit();
}
//
void macoctreeliquid3::post_initialize( bool initialized_from_file ) {
	//
	scoped_timer timer(this);
	//
	// Initialize scene
	auto initialize_func = reinterpret_cast<void(*)(const shape3 &m_shape, double m_dx)>(m_dylib.load_symbol("initialize"));
	if( initialize_func ) {
		timer.tick(); console::dump( "Initializing scene..." );
		initialize_func(m_shape,m_dx);
		console::dump( "Done. Took %s.\n", timer.stock("initialize_scene").c_str());
	}
	//
	// Get function pointers
	auto fluid_func = reinterpret_cast<double(*)(const vec3d &)>(m_dylib.load_symbol("fluid"));
	m_solid_func = reinterpret_cast<double(*)(const vec3d &)>(m_dylib.load_symbol("solid"));
	m_draw_func = reinterpret_cast<void(*)(graphics_engine &,double)>(m_dylib.load_symbol("draw"));
	m_moving_solid_func = reinterpret_cast<std::pair<double,vec3d>(*)(double time, const vec3d &p)>(m_dylib.load_symbol("moving_solid"));
	m_check_inject_func = reinterpret_cast<bool(*)(double, double, double, unsigned)>(m_dylib.load_symbol("check_inject"));
	m_inject_func = reinterpret_cast<bool(*)(const vec3d &, double, double, double, unsigned, double &, vec3d &)>(m_dylib.load_symbol("inject"));
	m_post_inject_func = reinterpret_cast<void(*)(double, double, double, unsigned, double&)>(m_dylib.load_symbol("post_inject"));
	auto velocity_func = reinterpret_cast<vec3d(*)(const vec3d &)>(m_dylib.load_symbol("velocity"));
	m_gravity_func = reinterpret_cast<vec3d(*)(double)>(m_dylib.load_symbol("gravity"));
	m_export_moving_poygon_func = reinterpret_cast<void(*)(polygon_list3 &)>(m_dylib.load_symbol("export_moving_poygon"));
	m_get_moving_polygon_transforms_func = reinterpret_cast<void(*)(double time, std::vector<vec3d> &, std::vector<vec3d> &)>(m_dylib.load_symbol("get_moving_polygon_transforms"));
	auto solid_visualize_func = reinterpret_cast<double(*)(const vec3d &)>(m_dylib.load_symbol("solid_visualize"));
	m_set_boundary_flux = reinterpret_cast<void(*)( double, Real [DIM3][2] )>(m_dylib.load_symbol("set_boundary_flux"));
	m_combined_solid_func = [&]( const vec3d &p ) {
		double value (1.0);
		if( m_solid_func ) value = std::min(value,m_solid_func(p));
		if( m_moving_solid_func ) value = std::min(value,m_moving_solid_func(m_timestepper->get_current_time(),p).first);
		return value;
	};
	//
	if( m_moving_solid_func ) {
		m_macoctreeproject.set_moving_solid([&]( const vec3d &p) {
			return m_moving_solid_func(m_timestepper->get_current_time(),p).first;
		});
	}
	//
	if( initialized_from_file ) {
		//
		console::dump( "Initialized from a file...\n" );
		//
		// Export moving polygon if necessary
		if( ! m_export_path.empty()) {
			export_moving_polygon();
		}
	} else {
		//
		timer.tick(); console::dump( ">>> Started initialization (%dx%dx%d)\n", m_shape[0], m_shape[1], m_shape[2] );
		//
		timer.tick(); console::dump( ">>> Adding layers...\n" );
		//
		m_grid_0.clear();
		m_grid_1.clear();
		m_accumulated_CFL = 0.0;
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
			shape3 shape = m_shape/n;
			console::dump( "Adding a layer of (%dx%dx%d)\n", shape[0], shape[1], shape[2] );
			m_grid_0.add_layer(shape,n*m_dx);
			m_grid_1.add_layer(shape,n*m_dx);
			n *= 2;
		}
		console::dump( "<<< Done. Took %s.\n", timer.stock("add_layers").c_str());
		//
		if( m_set_boundary_flux ) {
			flux_boundary_condition3 boundary_cond;
			m_set_boundary_flux(0.0,boundary_cond.velocity);
			m_grid_0.set_flux_boundary_condition(boundary_cond);
			m_grid_1.set_flux_boundary_condition(boundary_cond);
		}
		//
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
			m_grid->assign_levelset(fluid_func,m_combined_solid_func);
			//
			if( m_param.use_sizing_func ) {
				m_macoctreesizingfunc.compute_sizing_function(*m_grid_prev,*m_grid,0.0,m_combined_solid_func,[&]( const vec3d &p ) {
					return m_moving_solid_func ? m_moving_solid_func(m_timestepper->get_current_time(),p).second : vec3d();
				});
			}
			console::dump( "<<< Done. Took %s.\n", timer.stock("refinement"+std::to_string(count)).c_str());
			count ++;
		}
		//
		if( velocity_func ) {
			m_grid->set_velocity([&]( const vec3d &p, char dim ) {
				return velocity_func(p)[dim];
			});
		}
		//
		if( m_param.use_FLIP ) {
			shared_array3<Real> highres_fluid(m_shape);
			shared_macarray3<Real> highres_velocity(m_shape);
			//
			if( m_grid->reconstruct_fluid(highres_fluid(),[&]( int i, int j, int k, const Real &levelset_value ) {
				const double value = fluid_func(m_dx*vec3i(i,j,k).cell());
				return value < m_dx && value > -m_narrowband_depth;
			})) {
				m_grid->reconstruct_velocity(highres_velocity(),[&]( int dim, int i, int j, int k ) {
					const double value = fluid_func(m_dx*vec3i(i,j,k).face(dim));
					return value < m_dx && value > -m_narrowband_depth;
				});
				//
				m_flip->resample(highres_fluid(),m_combined_solid_func,highres_velocity());
			}
		}
		//
		if( solid_visualize_func ) {
			//
			array3<Real> solid_visualize_array(m_solid_shape.nodal());
			//
			const double sqrt3 = sqrt(3.0);
			timer.tick(); console::dump( ">>> Assigning solid level set...\n" );
			timer.tick(); console::dump( "Collecting active nodes..." );
			std::vector<vec3i> valid_positions[m_parallel.get_thread_num()];
			m_parallel.for_each(m_solid_shape.nodal(),[&](int i, int j, int k, int tid) {
				double value = (*solid_visualize_func)(m_solid_dx*vec3i(i,j,k).nodal());
				if( std::abs(value) < sqrt3*m_solid_dx ) {
					valid_positions[tid].push_back(vec3i(i,j,k));
				}
			});
			size_t count (0);
			for( const auto &e : valid_positions ) count += e.size();
			console::dump( "Done. Count sum = %zu. Took %s.\n", count, timer.stock("collect_active_solid_nodes").c_str());
			//
			if( count ) {
				timer.tick(); console::dump( "Assigning on nodes..." );
				for( const auto &e : valid_positions ) for( const auto &pi : e ) {
					double value = (*solid_visualize_func)(m_solid_dx*pi.nodal());
					solid_visualize_array.set(pi,value);
				}
				console::dump( "Done. Took %s.\n", timer.stock("assign_active_solid_nodes").c_str());
			}
			console::dump( "<<< Done. Took %s.\n", timer.stock("solid_levelset").c_str());
			//
			if( UI_interface::has_graphical_interface() ) {
				m_solid_visualize.copy(solid_visualize_array);
			}
			if( m_export_path.size()) do_export_solid_mesh(solid_visualize_array);
		} else {
			if( m_export_path.size()) do_export_empty_solid_mesh();
		}
		//
		// Export moving polygon
		if( ! m_export_path.empty()) {
			export_moving_polygon();
		}
		//
		// Segregate region
		if( m_param.volume_correction && m_param.regional_volume_correction ) {
			timer.tick(); console::dump( "Segregating region..." );
			//
			m_region_count = m_macoctreesegregator.segregate(*m_grid,m_regions);
			m_y_list.resize(m_region_count);
			m_macoctreesegregator.compute(*m_grid,m_regions,m_volumes);
			const double initial_volume = std::accumulate(m_volumes.begin(),m_volumes.end(),0.0);
			console::dump( "Done. regions=%u. initial volume = %.2e. Took %s.\n", m_region_count, initial_volume, timer.stock("segregate_region").c_str());
			//
			console::write("num_region",m_region_count);
			console::write("volume_sum",initial_volume);
		}
		//
		console::dump( "<<< Initialization finished. Took %s\n", timer.stock("initialization").c_str());
		//
		if( ! m_export_path.empty()) {
			export_mesh(0);
			render_mesh(0);
		}
	}
	//
	m_camera->set_bounding_box(vec3d().v,m_shape.box(m_dx).v);
}
//
static double grid_kernel( const vec3d &r, double dx ) {
	//
	double x = std::abs(r[0]) / dx;
	double y = std::abs(r[1]) / dx;
	double z = std::abs(r[2]) / dx;
	return std::max(0.0,1.0-x) * std::max(0.0,1.0-y) * std::max(0.0,1.0-z);
}
//
void macoctreeliquid3::drag( double x, double y, double z, double u, double v, double w ) {
}
//
void macoctreeliquid3::idle() {
	//
	scoped_timer timer(this);
	//
	// High-res fluid level set and velocity
	shared_array3<Real> highres_fluid(m_shape);
	shared_macarray3<Real> highres_velocity(m_shape);
	shared_macarray3<Real> save_velocity(m_shape);
	//
	// Define reconstruct test functions
	auto reconstruct_test_func_cell = [&]( int i, int j, int k, const Real &levelset_value ) {
		return levelset_value < (m_param.erode_width+2)*m_dx && levelset_value > -m_narrowband_depth-(m_param.erode_width+2)*m_dx;
	};
	//
	auto reconstruct_test_func_face = [&]( int dim, int i, int j, int k) {
		return highres_fluid->active(m_shape.clamp(i,j,k)) || highres_fluid->active(m_shape.clamp(i-(dim==0),j-(dim==1),k-(dim==2)));
	};
	//
	const unsigned step = m_timestepper->get_step_count()+1;
	timer.tick(); console::dump( ">>> %s step started...\n", console::nth(step).c_str());
	//
	// Compute time step
	timer.tick(); console::dump( "Computing time step...");
	double max_u_per_unit (0.0);
	for( size_t n=0; n<m_grid->velocity.size(); ++n ) max_u_per_unit = std::max(max_u_per_unit,(double)std::abs(m_grid->velocity[n]));
	const double dt = m_timestepper->advance(max_u_per_unit,m_grid->get_finest_dx());
	const double time = m_timestepper->get_current_time();
	const double CFL = m_timestepper->get_current_CFL();
	console::dump( "Done. time=%.2e. dt=%.2e,CFL=%.2f. Took %s\n", time, dt, CFL, timer.stock("compute_timestep").c_str());
	//
	// Begin (check) injecting fluid
	begin_inject_external_fluid(dt,time,step);
	//
	if( m_param.use_sizing_func ) {
		timer.tick(); console::dump( "Evaluating sizing function...");
		m_macoctreesizingfunc.compute_sizing_function(*m_grid_prev,*m_grid,dt,m_combined_solid_func,[&]( const vec3d &p ) {
				return m_moving_solid_func ? m_moving_solid_func(m_timestepper->get_current_time(),p).second : vec3d();
		});
		console::dump( "Done. Took %s\n", timer.stock("sizing_value").c_str());
	}
	//
	// Swap grid and remesh
	m_accumulated_CFL += CFL;
	std::swap(m_grid,m_grid_prev);
	if( m_do_inject || m_accumulated_CFL >= m_param.maximal_CFL_accumulation ) {
		//
		m_accumulated_CFL = 0.0;
		timer.tick(); console::dump( ">>> Remeshing...\n");
		if( m_param.use_sizing_func ) {
			if( m_do_inject ) {
				m_macoctreesizingfunc.activate_cells(*m_grid_prev,*m_grid,m_combined_solid_func,
					[&]( const vec3d &p ) {
						vec3d u; double value;
						m_inject_func(p,m_dx,dt,time,step,value,u);
						return value;
					});
			} else {
				m_macoctreesizingfunc.activate_cells(*m_grid_prev,*m_grid,m_combined_solid_func,nullptr);
			}
		} else {
			m_grid->activate_cells([&]( const vec3d &p ){
				//
				double inject_levelset (std::numeric_limits<double>::max());
				if( m_do_inject ) {
					vec3d u; double value;
					if( m_inject_func(p,m_dx,dt,time,step,value,u)) {
						inject_levelset = value;
					}
				}
				vec3d u (m_grid_prev->sample_velocity(p));
				return std::min(inject_levelset,m_grid_prev->sample_levelset(p-dt*u));
			},m_combined_solid_func);
		}
		//
		m_grid->balance_layers();
		m_grid->assign_indices();
		console::dump( "<<< Done. Took %s.\n", timer.stock("remeshing").c_str());
		//
		// Update fluid levelset
		if( m_param.use_FLIP && m_flip->get_particle_count()) {
			//
			if( m_grid_prev->reconstruct_fluid(highres_fluid(),reconstruct_test_func_cell)) {
				m_flip->update(m_combined_solid_func,highres_fluid(),time,false);
				highres_fluid->parallel_actives([&]( int i, int j, int k, auto &it ) {
					const auto &highres_layer = m_grid_prev->layers[0];
					m_grid_prev->levelset[highres_layer->active_cells(i,j,k)] = it();
				});
			}
		}
	} else {
		console::dump( "Accumulated CFL=%.2f. Copying previous grid...", m_accumulated_CFL );
		timer.tick();
		m_grid->copy(*m_grid_prev);
		console::dump( "Done. Took %s.\n", timer.stock("copy_grid").c_str());
	}
	console::write("num_accumulated_CFL",m_accumulated_CFL);
	//
	// Set boundary flux
	if( m_set_boundary_flux ) {
		flux_boundary_condition3 boundary_cond;
		m_set_boundary_flux(time,boundary_cond.velocity);
		m_grid->set_flux_boundary_condition(boundary_cond);
		console::dump( "Boundary flux = (%.2e,%.2e),(%.2e,%2.e),(%.2e,%2.e).\n",
			boundary_cond.velocity[0][0],boundary_cond.velocity[0][1],
			boundary_cond.velocity[1][0],boundary_cond.velocity[1][1],
			boundary_cond.velocity[2][0],boundary_cond.velocity[2][1] );
		//
		console::write("boundary_flux_00",boundary_cond.velocity[0][0]);
		console::write("boundary_flux_01",boundary_cond.velocity[0][1]);
		console::write("boundary_flux_10",boundary_cond.velocity[1][0]);
		console::write("boundary_flux_11",boundary_cond.velocity[1][1]);
		console::write("boundary_flux_20",boundary_cond.velocity[2][0]);
		console::write("boundary_flux_21",boundary_cond.velocity[2][1]);
	}
	//
	// Advect level set
	m_grid->assign_levelset([&]( const vec3d &p ) {
		vec3d u (m_grid_prev->sample_velocity(p));
		return m_grid_prev->sample_levelset(p-dt*u);
	},m_combined_solid_func);
	//
	// Advect region
	std::vector<uint_type> regions01; // Region after the advection
	if( m_param.volume_correction && m_param.regional_volume_correction ) {
		timer.tick(); console::dump( "Advecting region...");
		m_macoctreesegregator.backtrace(*m_grid_prev,*m_grid,dt,
			[&]( const vec3d &p ) {
				vec3d u (m_grid_prev->sample_velocity(p));
				return m_grid_prev->sample_levelset(p-dt*u);
			},m_regions,regions01);
		m_macoctreesegregator.extrapolate_jacobi(*m_grid,regions01);
		m_macoctreesegregator.prune(*m_grid,regions01);
		console::dump( "Done. Took %s.\n", timer.stock("advect_region").c_str());
	}
	//
	// Advect FLIP particles
	if( m_param.use_FLIP ) {
		m_flip->advect(
			m_combined_solid_func,
			[&](const vec3d &p){
				return vec3d(m_grid_prev->sample_velocity(p));
			},
			m_timestepper->get_current_time(),dt);
	}
	//
	// Advect velocity
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
	//
	if( m_param.use_FLIP && m_flip->get_particle_count()) {
		//
		// Define interpolation functions
		auto interpolate_fluid = [&]( const vec3d &p) {
			return m_grid->sample_levelset(p);
		};
		auto interpolate_velocity = [&]( const vec3d &p) {
			return m_grid->sample_velocity(p);
		};
		//
		// Mark bullet particles
		m_flip->mark_bullet(
			[&](const vec3d &p){ return interpolate_fluid(p); },
			[&](const vec3d &p){ return interpolate_velocity(p); },
			m_timestepper->get_current_time()
		);
		//
		// Correct positions
		m_flip->correct([&](const vec3d &p){ return interpolate_fluid(p); },highres_velocity());
		//
		// Splat momentum and mass of FLIP particles onto grids
		shared_macarray3<macflip3_interface::mass_momentum3> mass_and_momentum(m_shape);
		m_flip->splat(time,mass_and_momentum());
		//
		// Velocity overwrite
		mass_and_momentum->parallel_actives([&](char dim, int i, int j, int k, auto &it, int tn ) {
			auto &highres_layer = m_grid->layers[0];
			if( highres_layer->active_faces[dim].active(i,j,k)) {
				const auto value = it();
				double grid_mass = std::max(0.0,1.0-value.mass);
				Real &grid_velocity = m_grid->velocity[highres_layer->active_faces[dim](i,j,k)];
				grid_velocity = (grid_mass*grid_velocity+value.momentum) / (grid_mass+value.mass);
			}
		});
		//
		// Save the current velocity
		m_grid->reconstruct_velocity(save_velocity(),reconstruct_test_func_face);
	}
	//
	// Add gravity force
	timer.tick(); console::dump( "Adding gravity..." );
	const vec3d gravity = m_gravity_func ? m_gravity_func(time) : m_param.gravity;
	m_grid->iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
		m_grid->velocity[face_id.index] += dt*gravity[face_id.dim];
	});
	console::dump( "Done. Took %s\n", timer.stock("add_gravity").c_str());
	//
	std::vector<Real> target_volumes1;
	//
	if( m_param.volume_correction && m_param.regional_volume_correction ) {
		//
		timer.tick(); console::dump( "Computing new segregation..." );
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
		const double volume_sum = m_macoctreesegregator.compute_target_volume(*m_grid,regions01,m_volumes,regions11,volumes1,y_lists0,m_current_volumes,target_volumes1,m_y_list);
		m_regions = regions11;
		m_volumes = target_volumes1;
		m_macoctreesegregator.compute(*m_grid,m_regions,volumes1);
		//
		console::dump( "Done. regions=%u. volume_sum=%.2e. Took %s\n", m_region_count, volume_sum, timer.stock("compute_new_regional_volume").c_str());
		//
		console::write("num_region",m_region_count);
		console::write("volume_sum",volume_sum);
	}
	//
	// Inject external fluid
	do_inject_external_fluid(dt,time,step);
	//
	// Ending injecting fluid
	end_inject_external_fluid(dt,time,step);
	//
	// Add surface tension force
	if( m_param.surftens_k ) {
		timer.tick(); console::dump( "Adding surface tension force..." );
		m_grid->add_surfacetense_force(m_param.surftens_k,dt);
		console::dump( "Done. Took %s\n", timer.stock("add_surftens_force").c_str());
	}
	//
	// Assemble matrix
	m_macoctreeproject.assemble_matrix(*m_grid);
	//
	// Define solid velocity
	auto solid_velocity_func = [&]( const vec3d &p ) {
		return m_moving_solid_func ? m_moving_solid_func(time,p).second : vec3d();
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
		if(m_grid->reconstruct_fluid(highres_fluid(),reconstruct_test_func_cell)) {
			//
			m_grid->reconstruct_velocity(highres_velocity(),reconstruct_test_func_face);
			if( m_param.erode_width ) {
				timer.tick(); console::dump( "Eroding cells (%u cells wide)...", m_param.erode_width );
				highres_fluid->erode(m_param.erode_width);
				console::dump( "Done. Took %s\n", timer.stock("highres_fluid_erosion").c_str());
			}
			//
			// Reseed particles
			m_flip->resample(highres_fluid(),m_combined_solid_func,highres_velocity());
			//
			// Remove particles
			m_flip->remove([&](const vec3r &p, bool bullet) {
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
			m_flip->remove([](const vec3r &p, bool bullet){ return true; });
		}
	}
	//
	console::dump( "<<< %s step done. Took %s\n", console::nth(step).c_str(), timer.stock("simstep").c_str());
	//
	if( ! m_export_path.empty()) {
		int export_frame = m_timestepper->should_export_frame();
		if( m_param.debug_mode == 1 ) export_frame = step;
		if( export_frame ) {
			export_mesh(export_frame);
			render_mesh(export_frame);
			//
			if( m_param.transfer_file && console::get_root_path().size() && filesystem::is_exist("sync.sh")) {
				global_timer::pause();
				console::dump("Tranfering files...");
				std::string root_path (console::get_root_path());
				console::run( "./sync.sh %s",root_path.c_str());
				console::dump("Done.\n");
				global_timer::resume();
			}
		}
	}
	// Save state
	save_state();
}
//
void macoctreeliquid3::begin_inject_external_fluid( double dt, double time, unsigned step ) {
	m_do_inject = m_check_inject_func && m_inject_func && m_check_inject_func(m_dx,dt,time,step);
}
//
void macoctreeliquid3::do_inject_external_fluid( double dt, double time, unsigned step ) {
	//
	if( m_do_inject ) {
		//
		std::vector<double> total_injected (m_grid->parallel.get_thread_num());
		m_grid->iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
			//
			const vec3d p = m_grid->get_cell_position(cell_id);
			double value (0.0); vec3d u;
			//
			if( m_inject_func(p,m_dx,dt,time,step,value,u)) {
				const double v = m_grid->levelset[cell_id.index];
				m_grid->levelset[cell_id.index] = std::min(value,v);
				if( value < 0.0 && v > 0.0 ) {
					const double dx = m_grid->get_cell_dx(cell_id);
					total_injected[tid] += dx*dx*dx;
				}
			}
		});
		m_grid->iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
			const vec3d p = m_grid->get_face_position(face_id);
			const double dx = m_grid->get_face_dx(face_id);
			double value (0.0); vec3d u;
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
void macoctreeliquid3::end_inject_external_fluid( double dt, double time, unsigned step ) {
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
			m_grid->iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
				bool found (false);
				if( m_grid->levelset[cell_id.index] < 0.0 && ! regions_save[cell_id.index]) {
					uint_type new_index (0);
					m_grid->iterate_cell_neighbors(cell_id,[&]( char dim, const cell_id3& neighbor_cell_id ){
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
		m_grid->serial_iterate_active_cells([&]( const cell_id3 &cell_id ) {
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
void macoctreeliquid3::do_export_solid_mesh( const array3<Real> &solid ) {
	//
	scoped_timer timer(this);
	auto uv_coordinate_func = [&](const vec3d &p) {
		return vec2d(p[0],p[2]);
	};
	//
	std::string static_solids_directory_path = console::format_str("%s/static_solids",m_export_path.c_str());
	std::string path_wo_suffix = console::format_str("%s/levelset_solid",static_solids_directory_path.c_str());
	//
	if( ! filesystem::is_exist(static_solids_directory_path)) {
		//
		filesystem::create_directory(static_solids_directory_path);
		//
		if( array_utility3::levelset_exist(solid)) {
			//
			timer.tick(); console::dump( "Generating solid mesh..." );
			//
			std::vector<vec3d> vertices;
			std::vector<std::vector<size_t> > faces;
			m_solid_mesher->generate_mesh(solid,vertices,faces);
			//
			m_mesh_exporter->set_mesh(vertices,faces);
			std::vector<vec2d> uv_coordinates (vertices.size());
			for( unsigned n=0; n<vertices.size(); ++n ) {
				uv_coordinates[n] = uv_coordinate_func(vertices[n]);
			}
			//
			m_mesh_exporter->set_texture_coordinates(uv_coordinates);
			m_mesh_exporter->export_ply(console::format_str("%s.ply",path_wo_suffix.c_str()));
			m_mesh_exporter->export_mitsuba(console::format_str("%s.serialized",path_wo_suffix.c_str()));
			m_mesh_exporter->clear();
			//
			console::dump("Done. Took %s.\n", timer.stock("export_solid_mesh").c_str());
			//
		} else {
			//
			// Export dummy solid file
			do_export_empty_solid_mesh(true);
		}
	}
}
//
void macoctreeliquid3::export_moving_polygon() {
	//
	if( m_export_moving_poygon_func ) {
		//
		std::string moving_solids_directory_path = console::format_str("%s/moving_solids",m_export_path.c_str());
		if( ! filesystem::is_exist(moving_solids_directory_path)) {
			filesystem::create_directory(moving_solids_directory_path);
			//
			polygon_list3 polygons;
			m_export_moving_poygon_func(polygons);
			//
			for( unsigned n=0; n<polygons.size(); ++n ) {
				//
				const std::vector<vec3d> &vertices = polygons[n].first;
				const std::vector<std::vector<size_t> > &faces = polygons[n].second;
				std::string path_wo_suffix = console::format_str("%s/polygon_%d",moving_solids_directory_path.c_str(),n);
				//
				m_mesh_exporter->set_mesh(vertices,faces);
				m_mesh_exporter->export_ply(console::format_str("%s.ply",path_wo_suffix.c_str()));
				m_mesh_exporter->export_mitsuba(console::format_str("%s.serialized",path_wo_suffix.c_str()));
				m_mesh_exporter->clear();
			}
		}
	}
}
//
void macoctreeliquid3::do_export_empty_solid_mesh( bool force ) {
	//
	std::string static_solids_directory_path = console::format_str("%s/static_solids",m_export_path.c_str());
	std::string path_wo_suffix = console::format_str("%s/levelset_solid",static_solids_directory_path.c_str());
	//
	if( force || ! filesystem::is_exist(static_solids_directory_path)) {
		//
		if( ! filesystem::is_exist(static_solids_directory_path)) filesystem::create_directory(static_solids_directory_path);
		//
		// Export dummy solid file
		std::vector<vec3d> vertices;
		std::vector<std::vector<size_t> > faces(1);
		vertices.push_back(vec3d(1e3,1e3,1e3));
		vertices.push_back(vec3d(1e3+1.0,1e3,1e3));
		vertices.push_back(vec3d(1e3,1e3,1e3+1.0));
		faces[0].push_back(0);
		faces[0].push_back(1);
		faces[0].push_back(2);
		m_mesh_exporter->set_mesh(vertices,faces);
		m_mesh_exporter->export_ply(console::format_str("%s.ply",path_wo_suffix.c_str()));
		m_mesh_exporter->export_mitsuba(console::format_str("%s.serialized",path_wo_suffix.c_str()));
		m_mesh_exporter->clear();
	}
}
//
void macoctreeliquid3::export_mesh( int frame ) {
	//
	scoped_timer timer(this);
	timer.tick(); console::dump( ">>> Exporting %s mesh...\n", console::nth(frame).c_str() );
	//
	std::vector<vec3d> vertices;
	std::vector<std::vector<size_t> > faces;
	//
	timer.tick(); console::dump( "Generating surface mesh..." );
	m_macoctreemesher.generate_mesh(*m_grid,m_grid->levelset,0.0,m_solid_func,vertices,faces,false);
	//
	if( m_param.remove_quater ) {
		std::vector<std::vector<size_t> > faces_save(faces);
		faces.clear();
		for( const auto &f : faces_save ) {
			vec3d v_sum;
			for( const auto &idx : f ) {
				v_sum += vertices[idx];
			}
			v_sum /= f.size();
			if( v_sum[0] > 0.5 || v_sum[2] < m_param.z-0.25*m_dx ) faces.push_back(f);
		}
	}
	//
	for( unsigned n=frame; n<1000; ++n ) {
		//
		std::string next_mesh_ply = console::format_str("%s/%d_mesh.ply",m_export_path.c_str(),n);
		std::string next_mesh_serialized = console::format_str("%s/%d_mesh.serialized",m_export_path.c_str(),n);
		std::string next_particle_path = console::format_str("%s/%d_particles.dat",m_export_path.c_str(),n);
		std::string next_translation_path = console::format_str("%s/%d_transforms.dat",m_export_path.c_str(),n);
		std::string svg_export_path = console::get_root_path() + "/SVG";
		std::string next_svg_path = svg_export_path + "/" + "output_" + std::to_string(n) + ".svg";
		std::string next_jpg_path = svg_export_path + "/" + "output_" + std::to_string(n) + ".jpg";
		//
		if( filesystem::is_exist(next_mesh_ply)) filesystem::remove_file(next_mesh_ply);
		if( filesystem::is_exist(next_mesh_serialized)) filesystem::remove_file(next_mesh_serialized);
		if( filesystem::is_exist(next_particle_path)) filesystem::remove_file(next_particle_path);
		if( filesystem::is_exist(next_translation_path)) filesystem::remove_file(next_translation_path);
		if( filesystem::is_exist(next_svg_path)) filesystem::remove_file(next_svg_path);
		if( filesystem::is_exist(next_jpg_path)) filesystem::remove_file(next_jpg_path);
	}
	//
	m_mesh_exporter->set_mesh(vertices,faces);
	m_mesh_exporter->export_ply(console::format_str("%s/%d_mesh.ply",m_export_path.c_str(),frame));
	m_mesh_exporter->export_mitsuba(console::format_str("%s/%d_mesh.serialized",m_export_path.c_str(),frame));
	console::dump( "Done. Took %s.\n", timer.stock("surface_mesh_gen").c_str());
	//
	timer.tick(); console::dump( "Generating enclosed mesh..." );
	m_macoctreemesher.generate_mesh(*m_grid,m_grid->levelset,0.0,m_solid_func,vertices,faces,true);
	m_mesh_exporter->set_mesh(vertices,faces);
	m_mesh_exporter->export_ply(console::format_str("%s/%d_mesh_enclosed.ply",m_export_path.c_str(),frame));
	m_mesh_exporter->export_mitsuba(console::format_str("%s/%d_mesh_enclosed.serialized",m_export_path.c_str(),frame));
	console::dump( "Done. Took %s.\n", timer.stock("enclosed_surface_mesh_gen").c_str());
	m_mesh_exporter->clear();
	//
	timer.tick(); console::dump( "Generating particles..." );
	std::vector<macflip3_interface::particle3> particles = m_flip->get_particles();
	std::string particle_path = console::format_str("%s/%d_particles.dat",m_export_path.c_str(),frame);
	FILE *fp = fopen(particle_path.c_str(),"wb");
	//
	auto particle_test = [&]( const macflip3_interface::particle3 &particle ) {
		if( m_param.remove_quater ) {
			return particle.p[0] > 0.5 || particle.p[2] < m_param.z-0.25*m_dx;
		} else {
			return true;
		}
	};
	//
	size_t size (0);
	for( size_t n=0; n<particles.size(); ++n ) {
		if( particles[n].bullet ) {
			if( particle_test(particles[n]) && m_grid->sample_levelset(particles[n].p) > particles[n].r ) {
				++ size;
			} else {
				particles[n].bullet = false;
			}
		}
	}
	fwrite(&size,1,sizeof(unsigned),fp);
	for( size_t n=0; n<particles.size(); ++n ) {
		if( particles[n].bullet ) {
			float position[3] = { (float)particles[n].p.v[0],
								  (float)particles[n].p.v[1],
								  (float)particles[n].p.v[2] };
			float radius = particles[n].r;
			assert( ! utility::is_nan(radius));
			fwrite(position,3,sizeof(float),fp);
			fwrite(&radius,1,sizeof(float),fp);
		}
	}
	fclose(fp);
	//
	if( m_get_moving_polygon_transforms_func ) {
		//
		const double time = m_timestepper->get_current_time();
		std::vector<vec3d> translations, rotations;
		m_get_moving_polygon_transforms_func(time,translations,rotations);
		//
		std::string translation_path = console::format_str("%s/%d_transforms.dat",m_export_path.c_str(),frame);
		FILE *fp = fopen(translation_path.c_str(),"wb");
		const unsigned number = translations.size();
		fwrite(&number,1,sizeof(unsigned),fp);
		for( unsigned n=0; n<translations.size(); ++n ) {
			const vec3d &t = translations[n];
			vec3d r = rotations[n];
			float len = r.len();
			if( ! len ) {
				r[0] = 1.0;
				len = 1.0;
			}
			const float t_float[3] = { (float)t.v[0],
	 								   (float)t.v[1],
								 	   (float)t.v[2] };
			const float r_float[3] = { (float)r.v[0] / len,
									   (float)r.v[1] / len,
									   (float)r.v[2] / len };
			fwrite(t_float,3,sizeof(float),fp);
			fwrite(&len,1,sizeof(float),fp);
			fwrite(r_float,3,sizeof(float),fp);
		}
		fclose(fp);
	}
	console::dump( "Done. Generated %u particles. Took %s.\n", size, timer.stock("particle_gen").c_str());
	//
	console::dump( "<<< Done. Took %s.\n", timer.stock("export_mesh").c_str());
	//
	if( m_param.export_svg ) {
		//
		std::string svg_export_path = console::get_root_path() + "/SVG";
		if( ! filesystem::is_exist(svg_export_path)) {
			filesystem::create_directory(svg_export_path);
		}
		//
		const double y = m_shape[1] / (double) m_shape[0];
		const double image_width (10240.0);
		m_svg_writer->setup_graphics();
		m_svg_writer->set_2D_coordinate(0.0,m_dx*m_shape[0],0.0,m_dx*m_shape[1]);
		m_svg_writer->set_viewport(0,0,image_width,image_width*y);
		m_svg_writer->clear();
		//
		const double z = m_param.z-0.5*m_dx;
		m_grid->draw_grid(*m_svg_writer.get(),z,true);
		//
		std::string final_path = svg_export_path + "/" + "output_" + std::to_string(frame) + ".svg";
		m_svg_writer->const_send_message("write",(char *)final_path.c_str());
		std::string jpg_path = svg_export_path + "/" + "output_" + std::to_string(frame) + ".jpg";
		console::system("convert %s %s",final_path.c_str(),jpg_path.c_str());
		//
		std::vector<vec3d> vertices;
		std::vector<vec2d> uv_coord;
		std::vector<std::vector<size_t> > faces;
		//
		const double scale_x = m_dx * m_shape[0];
		const double scale_y = m_dx * m_shape[1];
		const vec3d v0 (0.0,0.0,z);
		const vec3d v1 (scale_x,scale_y,z);
		//
		vertices.push_back(vec3d(v0[0],v0[1],z));
		vertices.push_back(vec3d(v0[0],v1[1],z));
		vertices.push_back(vec3d(v1[0],v1[1],z));
		vertices.push_back(vec3d(v1[0],v0[1],z));
		//
		const vec2d p0 (0.0,0.0);
		const vec2d p1 (1.0,1.0);
		uv_coord.push_back(vec2d(p0[0],1.0-p0[1]));
		uv_coord.push_back(vec2d(p0[0],1.0-p1[1]));
		uv_coord.push_back(vec2d(p1[0],1.0-p1[1]));
		uv_coord.push_back(vec2d(p1[0],1.0-p0[1]));
		//
		faces.push_back({0,1,2,3});
		//
		m_mesh_exporter->set_mesh(vertices,faces);
		m_mesh_exporter->set_texture_coordinates(uv_coord);
		m_mesh_exporter->export_ply(console::format_str("%s/xy_plane.ply",m_export_path.c_str()));
		m_mesh_exporter->clear();
	}
}
//
void macoctreeliquid3::render_mesh( unsigned frame ) const {
	//
	scoped_timer timer(this);
	global_timer::pause();
	//
	if(console::get_root_path().size()) {
		//
		std::string mitsuba_path = console::get_root_path() + "/mitsuba";
		std::string copy_from_path = filesystem::find_resource_path("octreeliquid","mitsuba");
		if( ! filesystem::is_exist(mitsuba_path)) {
			if( filesystem::is_exist(copy_from_path)) {
				console::run( "cp -r %s %s", copy_from_path.c_str(), mitsuba_path.c_str());
			} else {
				console::dump( "Could not lcoate mitsuba files (%s).\n", copy_from_path.c_str());
				exit(0);
			}
		}
		//
		auto remove_future_fles = [&]( std::string name ) {
			std::string image_dir_path = console::get_root_path() + "/" + name + "_img";
			for( unsigned n=frame; n<1000; ++n ) {
				std::string next_exr_path = console::format_str("%s/%d_%s.exr",image_dir_path.c_str(),n,name.c_str());
				std::string next_png_path = console::format_str("%s/%d_%s.png",image_dir_path.c_str(),n,name.c_str());
				if( filesystem::is_exist(next_exr_path)) filesystem::remove_file(next_exr_path);
				if( filesystem::is_exist(next_png_path)) filesystem::remove_file(next_png_path);
				if( name == "grid") {
					std::string next_composite_path = console::format_str("%s/%d_%s_composite.jpg",image_dir_path.c_str(),n,name.c_str());
					if( filesystem::is_exist(next_composite_path)) filesystem::remove_file(next_composite_path);
				}
			}
		};
		//
		remove_future_fles("mesh");
		remove_future_fles("transparent");
		remove_future_fles("wireframe");
		remove_future_fles("grid");
		//
		if( m_param.render_mesh || m_param.render_wireframe ) {
			for( int n=0; n<2; ++n ) {
				if( m_param.render_wireframe == false && n == 0 ) continue;
				if( m_param.render_mesh == false && n == 1 ) continue;
				//
				double liquid_color[3];
				if( n == 0 ) {
					liquid_color[0] = 181/255.0;
					liquid_color[1] = 209/255.0;
					liquid_color[2] = 253/255.0;
				} else {
					liquid_color[0] = 0.5;
					liquid_color[1] = 0.5;
					liquid_color[2] = 1.0;
				}
				std::string render_command = console::format_str("cd %s; /usr/bin/python render.py %d %s %g %g %g %d %g %g %g %g %g %g %g",
						mitsuba_path.c_str(),
						frame,
						n == 0 ? "wireframe" : "mesh",
						liquid_color[0], liquid_color[1], liquid_color[2],
						m_param.render_sample_count,
						m_param.target[0], m_param.target[1], m_param.target[2],
						m_param.origin[0], m_param.origin[1], m_param.origin[2],
						0.1 * m_dx );
				//
				console::dump("Running command: %s\n", render_command.c_str());
				console::system(render_command.c_str());
			}
		}
		//
		if( m_param.render_transparent ) {
			//
			std::string render_command = console::format_str("cd %s; /usr/bin/python render.py %d transparent %g %g %g %d %g %g %g %g %g %g %g",
					mitsuba_path.c_str(),
					frame,
					0.5, 0.5, 1.0,
					m_param.render_transparent_sample_count,
					m_param.target[0], m_param.target[1], m_param.target[2],
					m_param.origin[0], m_param.origin[1], m_param.origin[2],
					0.1 * m_dx );
			//
			console::dump("Running command: %s\n", render_command.c_str());
			console::system(render_command.c_str());
		}
		//
		if( m_param.render_grid && m_param.export_svg ) {
			//
			std::string render_command = console::format_str("cd %s; /usr/bin/python render.py %d grid %g %g %g %d %g %g %g %g %g %g %g",
					mitsuba_path.c_str(),
					frame,
					0.5, 0.5, 1.0,
					m_param.render_sample_count,
					m_param.target[0], m_param.target[1], m_param.target[2],
					m_param.origin[0], m_param.origin[1], m_param.origin[2],
					0.1 * m_dx );
			console::dump("Running command: %s\n", render_command.c_str());
			console::system(render_command.c_str());
			//
			render_command = console::format_str("cd %s; /usr/bin/python render.py %d grid_without %g %g %g %d %g %g %g %g %g %g %g",
					mitsuba_path.c_str(),
					frame,
					0.5, 0.5, 1.0,
					m_param.render_sample_count,
					m_param.target[0], m_param.target[1], m_param.target[2],
					m_param.origin[0], m_param.origin[1], m_param.origin[2],
					0.1 * m_dx );
			console::dump("Running command: %s\n", render_command.c_str());
			console::system(render_command.c_str());
			//
			std::string exp = console::get_root_path();
			for( unsigned frame1=0; frame1<=frame; ++frame1 ) {
				std::string final_path = console::format_str("%s/grid_img/%d_grid_composite.jpg",exp.c_str(),frame1);
				if( ! filesystem::is_exist(final_path)) {
					std::string composite_command = console::format_str("composite -blend 50 %s/grid_img/%d_grid.png %s/grid_without_img/%d_grid_without.png %s",
							exp.c_str(),frame1,exp.c_str(),frame1,final_path.c_str());
					console::dump("Running command: %s\n", composite_command.c_str());
					console::system(composite_command);
				}
				//
				std::string svg_export_path = console::get_root_path() + "/SVG";
				std::string mip_path = svg_export_path + "/" + "output_" + std::to_string(frame) + ".mip";
				if( filesystem::is_exist(mip_path)) {
					console::system("rm %s",mip_path.c_str());
				}
			}
			//
			if( frame && frame % 10 == 0 ) {
				console::system("rm -f %s/grid_img/grid_composite.mp4; avconv -r 60 -i %s/grid_img/%%d_grid_composite.jpg -b:v 12000k -pix_fmt yuv420p %s/grid_img/grid_composite.mp4",exp.c_str(),exp.c_str(),exp.c_str());
			}
		}
	}
	//
	global_timer::resume();
}
//
void macoctreeliquid3::save_state() {
	//
	global_timer::pause();
	scoped_timer timer(this);
	//
	if( filesystem::is_exist("quit")) {
		console::dump( "Quit on save found. Quitting...\n" );
		m_should_quit_on_save = true;
	}
	if( m_param.save_interval && console::get_root_path().size()) {
		unsigned step = m_timestepper->get_step_count();
		if( m_should_quit_on_save || (step % m_param.save_interval == 0) ) {
			//
			for( unsigned n=step; n<10000; ++n ) {
				std::string next_state = console::get_root_path()+"/state/"+std::to_string(n)+"_state.dat";
				if( filesystem::is_exist(next_state)) filesystem::remove_file(next_state);
			}
			//
			std::string path = console::get_root_path()+"/state";
			if( ! filesystem::is_exist(path)) filesystem::create_directory(path);
			path += "/"+std::to_string(step)+"_state.dat";
			filestream file(path,filestream::WRITE);
			recursive_serialize(file);
			timer.tick(); console::dump( "Saving %s state (%s)...", console::nth(step).c_str(), path.c_str());
			console::dump("Done.\n");
		}
	}
	global_timer::resume();
}
//
void macoctreeliquid3::draw( graphics_engine &g ) const {
	//
	const double time = m_timestepper->get_current_time();
	//
	// Draw from scene library
	if( m_draw_func ) m_draw_func(g,time);
	//
	g.color4(0.9,0.6,0.3,0.5);
	m_solid_gridvisualizer->draw_levelset(g,m_solid_visualize);
	//
	// Draw FLIP
	if( m_param.use_FLIP ) m_flip->draw(g,time);
	//
	// Draw surface
	std::vector<vec3d> vertices;
	std::vector<std::vector<size_t> > faces;
	m_macoctreemesher.generate_mesh(*m_grid,m_grid->levelset,0.0,m_solid_func,vertices,faces,false);
	//
	g.color4(1.0,1.0,1.0,0.5);
	for( unsigned i=0; i<faces.size(); i++ ) {
		g.begin(graphics_engine::MODE::LINE_LOOP);
		for( unsigned j=0; j<faces[i].size(); j++ ) g.vertex3v(vertices[faces[i][j]].v);
		g.end();
	}
}
//
extern "C" module * create_instance() {
	return new macoctreeliquid3;
}