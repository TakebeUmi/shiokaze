/*
**	macoctreeproject2.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on November 11, 2019. All rights reserved.
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
#include "macoctreeproject2.h"
#include <shiokaze/math/RCMatrix_utility.h>
#include <pcgsolver/pcg_solver.h>
#include <numeric>
#include <cassert>
//
SHKZ_USING_NAMESPACE
using namespace macotreeliquid2_namespace;
//
void macoctreeproject2::configure( configuration &config ) {
	//
	configuration::auto_group group(config,*this);
	//
	config.get_bool("VolumeCorrection",m_param.volume_correction,"Whether to perform volume correction");
	config.get_bool("VolumeCorrectSkipSurface",m_param.volume_correct_skip_surfaces,"Skip near surfaces for volume correct");
	config.set_default_bool("ForceGlobalResidual",true);
	config.get_bool("CheckSymmetric",m_param.check_symmetric,"Check matrix symmetircity");
	config.get_bool("CheckPositiveDiag",m_param.check_positive_diag,"Check the all the diag elements are positive");
	config.get_bool("FixDivergence",m_param.fix_divergence,"Fix divergence due to the moving solid");
	config.get_bool("DebugAssemble",m_param.debug_assemble,"Debug mode for assembling matrix");
	config.get_double("VolumeRecoverRatio",m_param.volume_recover_ratio,"Volume recover ratio for volume correction");
	config.get_double("MinimalVolumeCorrectTimeStep",m_param.minimal_volume_correct_timestep,"Minimal time step size for volume correction");
	config.get_double("NP",m_param.np,"Constant np for volume correction");
	config.get_bool("RemoveOneDegreesOfFreedom",m_param.remove_one_degrees_of_freedom,"Remove one degrees of freedom");
	config.get_integer("SpecialBoundaryCondition",m_param.special_boundary_condition,"Type of special boundary condition");
}
//
void macoctreeproject2::post_initialize( bool initialized_from_file ) {
	//
	m_initial_volume = 0.0;
	m_y_prev = 0.0;
}
//
void macoctreeproject2::assemble_matrix( grid2 &grid ) {
	//
	grid.compute_cell_map();
	grid.compute_face_map();
	//
	uint_type cell_count = grid.valid_cell_count;
	uint_type face_count = grid.valid_face_count;
	//
	if( m_matrix.allocated ) {
		m_matrix.Lhs->initialize(cell_count,cell_count);
		if( m_param.debug_assemble ) {
			m_matrix.G->initialize(face_count,cell_count);
			m_matrix.D->initialize(cell_count,face_count);
		}
	} else {
		m_matrix.Lhs = m_factory->allocate_matrix(cell_count,cell_count);
		if( m_param.debug_assemble ) {
			m_matrix.G = m_factory->allocate_matrix(face_count,cell_count);
			m_matrix.D = m_factory->allocate_matrix(cell_count,face_count);
		}
		m_matrix.allocated = true;
	}
	//
	if( m_param.debug_assemble ) {
		//
		// Scaled gradient matrix
		grid.iterate_active_faces([&]( const face_id2 &face_id, int tid ) {
			if( grid.face_map[face_id.index] ) {
				uint_type row = grid.face_map[face_id.index]-1;
				grid.get_scaled_gradient(face_id,[&]( const cell_id2 &cell_id, double value, const grid2::gradient_info2 &info ) {
					if( value ) {
						assert( grid.cell_map[cell_id.index] );
						uint_type column = grid.cell_map[cell_id.index]-1;
						m_matrix.G->add_to_element(row,column,value);
					}
				});
			}
		});
		//
		// Divergence matrix
		grid.iterate_active_cells([&]( const cell_id2 &cell_id, int tid ) {
			if( grid.cell_map[cell_id.index] ) {
				uint_type row = grid.cell_map[cell_id.index]-1;
				grid.get_divergence(cell_id,[&]( const face_id2 &face_id, double value0, double value1 ) {
					if( value0 && grid.face_map[face_id.index] ) {
						uint_type column = grid.face_map[face_id.index]-1;
						m_matrix.D->add_to_element(row,column,value0);
					}
				});
			}
		});
		//
		m_matrix.D->multiply(m_matrix.G.get(),m_matrix.Lhs.get());
		RCMatrix_utility<size_t,double>::report(m_matrix.Lhs.get(),"Lhs");
		//
	} else {
		//
		// Directily assemble Lhs matrix
		grid.iterate_active_cells([&]( const cell_id2 &cell_id, int tid ) {
			if( grid.cell_map[cell_id.index] ) {
				uint_type row0 = grid.cell_map[cell_id.index]-1;
				double diag (0.0);
				grid.get_divergence(cell_id,[&]( const face_id2 &face_id, double value0, double value1 ) {
					if( grid.face_map[face_id.index] ) {
						if( value0 ) {
							grid.get_scaled_gradient(face_id,[&]( const cell_id2 &cell_neigh_id, double value, const grid2::gradient_info2 &info ) {
								if( value ) {
									uint_type row1 = grid.cell_map[cell_neigh_id.index]-1;
									const double v = value0*value;
									if( row0 == row1 ) diag += v;
									else m_matrix.Lhs->add_to_element(row0,row1,v);
								}
							});
						}
					}
				});
				m_matrix.Lhs->add_to_element(row0,row0,diag);
			}
		});
		RCMatrix_utility<size_t,double>::report(m_matrix.Lhs.get(),"Lhs");
	}
	//
	if( m_param.check_symmetric ) {
		const double symm_error = RCMatrix_utility<size_t,double>::symmetricity_error(m_matrix.Lhs.get());
		assert( symm_error == 0.0 );
	}
	if( m_param.check_positive_diag ) {
		const double min_diag = RCMatrix_utility<size_t,double>::min_diag(m_matrix.Lhs.get());
		assert( min_diag > 0.0 );
	}
	m_matrix.assembled = true;
}
//
void macoctreeproject2::clear_matrix() {
	//
	if( m_matrix.allocated ) {
		m_matrix.Lhs.reset();
		if( m_param.debug_assemble ) {
			m_matrix.G.reset();
			m_matrix.D.reset();
		}
		m_matrix.allocated = false;
	}
	m_matrix.assembled = false;
}
//
void macoctreeproject2::set_moving_solid( std::function<double(const vec2d &p)> moving_solid_func ) {
	m_moving_solid_func = moving_solid_func;
}
//
void macoctreeproject2::project( grid2 &grid, double dt, std::function<vec2d(const vec2d &p)> solid_velocity, std::vector<Real> *pressure ) {
	//
	std::vector<uint_type> regions;
	std::vector<Real> current_volumes;
	std::vector<Real> target_volumes;
	std::vector<Real> y_list;
	//
	project(grid,dt,0,regions,current_volumes,target_volumes,y_list,solid_velocity,pressure);
}
//
void macoctreeproject2::project( grid2 &grid, double dt,
								 size_t region_count,
								 const std::vector<uint_type> &regions,
								 const std::vector<Real> &current_volumes,
								 const std::vector<Real> &target_volumes,
								 std::vector<Real> &y_list,
								 std::function<vec2d(const vec2d &p)> solid_velocity,
								 std::vector<Real> *pressure_vector ) {
	//
	double rhs_correct (0.0);
	std::vector<Real> rhs_corrects(region_count);
	//
	if( m_param.volume_correction ) {
		const double dt0 = std::max(m_param.minimal_volume_correct_timestep,dt);
		if( region_count ) {
			for( unsigned n=0; n<region_count; ++n ) {
				double x = (current_volumes[n]-target_volumes[n])/target_volumes[n];
				double y = y_list[n] + x*dt0; y_list[n] = y;
				double kp = -log(1.0-m_param.volume_recover_ratio)/(m_param.np*dt0);
				double ki = kp*kp/16.0;
				rhs_corrects[n] = -(kp*x+ki*y)/(x+1.0);
				if( rhs_corrects[n] ) rhs_correct = 1.0;
			}
		} else {
			//
			double current_volume = grid.get_volume();
			if( ! m_initial_volume ) m_initial_volume = current_volume;
			double x = (current_volume-m_initial_volume)/m_initial_volume;
			double y = m_y_prev + x*dt0; m_y_prev = y;
			double kp = -log(1.0-m_param.volume_recover_ratio)/(m_param.np*dt0);
			double ki = kp*kp/16.0;
			rhs_correct = -(kp*x+ki*y)/(x+1.0);
		}
	}
	//
	auto is_surface_cell = [&]( const cell_id2 &cell_id ) {
		//
		if( m_param.volume_correct_skip_surfaces ) {
			bool touching_air (false);
			grid.iterate_cell_neighbors(cell_id,[&]( char dim, const cell_id2 &cell_n_id ) {
				if( grid.levelset[cell_n_id.index] > 0.0 ) touching_air = true;
			});
			return touching_air;
		} else {
			return false;
		}
	};
	//
	assert( m_matrix.assembled );
	auto rhs = m_factory->allocate_vector(grid.valid_cell_count);
	auto result = m_factory->allocate_vector(grid.valid_cell_count);
	//
	grid.iterate_active_cells([&]( const cell_id2 &cell_id, int tid ) {
		if( grid.cell_map[cell_id.index] ) {
			uint_type row = grid.cell_map[cell_id.index]-1;
			if( rhs_correct ) {
				bool solid_cell (false);
				for( char dim : DIMS2 ) {
					grid.iterate_face_neighbors(cell_id,dim,[&]( const face_id2 &face_id ){
						if( grid.area[face_id.index] < 1.0 ) solid_cell = true;
					});
					if( solid_cell ) break;
				}
				if( rhs_correct && ! solid_cell && ! is_surface_cell(cell_id)) {
					const double dx = grid.get_cell_dx(cell_id);
					if( region_count ) {
						rhs->add(row,(dx*dx)*rhs_corrects[regions[cell_id.index]-1]);
					} else {
						rhs->add(row,(dx*dx)*rhs_correct);
					}
				}
			}
			//
			if( m_param.fix_divergence ) {
				double err (0.0);
				if( solid_velocity ) {
					grid.get_unmofidied_divergence(cell_id,[&]( const face_id2 &face_id, double value ) {
						err += value*grid.sample_solid_face_velocity(face_id,solid_velocity);
					});
				}
				bool boundary_flag (false);
				for( int dim : DIMS2 ) {
					boundary_flag = cell_id.pi[dim] == 0 || cell_id.pi[dim] == grid.layers[cell_id.depth]->shape[dim]-1;
					if( boundary_flag ) break;
				}
				if( ! boundary_flag ) rhs->add(row,-err);
			}
			//
			grid.get_divergence(cell_id,[&]( const face_id2 &face_id, double value0, double value1 ) {
				bool add_divergence (true);
				if( m_param.special_boundary_condition == 1 && cell_id.pi[0] == grid.layers[cell_id.depth]->shape[0]-1 && face_id.dim == 0 ) {
					add_divergence = false;
				}
				if( add_divergence ) {
					if( solid_velocity && value1 ) {
						rhs->add(row,value1*grid.sample_solid_face_velocity(face_id,solid_velocity));
					}
					if( value0 ) {
						rhs->add(row,value0*grid.velocity[face_id.index]);
					}
				}
			});
			//
			if( grid.flux_boundary_condition.has_flux()) {
				const double dx = grid.get_cell_dx(cell_id);
				const vec2d p = grid.get_cell_position(cell_id);
				if( m_moving_solid_func(p) > dx ) for( int dim : DIMS2 ) {
					if( cell_id.pi[dim] == 0 ) {
						rhs->add(row,(dx)*grid.flux_boundary_condition.velocity[dim][0]);
					} else if( cell_id.pi[dim] == grid.layers[cell_id.depth]->shape[dim]-1 ) {
						rhs->add(row,-(dx)*grid.flux_boundary_condition.velocity[dim][1]);
					}
				}
			}
		}
	});
	//
	auto compute_vector_kind = [&]( std::vector<unsigned char> &result ) {
		unsigned total_kinds (0);
		result.resize(grid.valid_cell_count);
		grid.serial_iterate_active_cells([&]( const cell_id2 &cell_id ) {
			if( grid.cell_map[cell_id.index] ) {
				uint_type row = grid.cell_map[cell_id.index]-1;
				result[row] = cell_id.depth;
				total_kinds = std::max(total_kinds,(unsigned)result[row]);
			}
		});
		return total_kinds;
	};
	//
	std::vector<unsigned char> vector_kind; compute_vector_kind(vector_kind);
	m_solver->register_vector_norm_kind(vector_kind);
	//
	if( m_param.remove_one_degrees_of_freedom ) {
		m_matrix.Lhs->clear(0);
		m_matrix.Lhs->add_to_element(0,0,1.0);
		rhs->set(0,0);
		m_parallel.for_each(grid.valid_cell_count,[&]( size_t row ) {
			m_matrix.Lhs->for_each(row,[&]( size_t column, double& value ) {
				if( column == 0 ) value = 0.0;
			});
		});
	}
	auto status = m_solver->solve(m_matrix.Lhs.get(),rhs.get(),result.get());
	//
	std::string str;
	if( status.vector_reresid.empty()) {
		str = console::format_str("%.2e",status.reresid);
	} else if( status.vector_reresid.size() <= 4 ) {
		for( const auto &e : status.vector_reresid ) {
			str += " " + console::format_str("%.2e",e);
		}
	} else {
		double max_reresid (0.0);
		double min_reresid (-1.0);
		unsigned slot_max (0), slot_min(0);
		for( int n=0; n<status.vector_reresid.size(); ++n ) {
			const auto &e = status.vector_reresid[n];
			if( e > max_reresid ) {
				max_reresid = e;
				slot_max = n;
			}
			if( min_reresid < 0.0 ) min_reresid = max_reresid;
			else if( e < min_reresid ) {
				min_reresid = e;
				slot_min = n;
			}
		}
		str = console::format_str("%.2e (num=%u, min[%u]=%.2e, max[%u]=%.2e)",status.reresid,status.vector_reresid.size(),min_reresid,slot_min,max_reresid,slot_max);
	}
	console::dump( "Iteration = %d. Reresid=%s\n", status.count, str.c_str() );
	//
	for( uint_type n=0; n<grid.face_count; ++n ) {
		if( ! grid.face_map[n] ) grid.velocity[n] = 0.0;
	}
	grid.iterate_active_faces([&]( const face_id2 &face_id, int tid ) {
		grid.get_scaled_gradient(face_id,[&]( const cell_id2 &cell_id, double value, const grid2::gradient_info2 &info ) {
			if( grid.face_map[face_id.index] && grid.cell_map[cell_id.index] && value ) {
				grid.velocity[face_id.index] -= value * result->at(grid.cell_map[cell_id.index]-1);
			}
		});
	});
	//
	if( pressure_vector ) {
		pressure_vector->resize(grid.cell_count);
		for( uint_type n=0; n<grid.cell_count; ++n ) {
			if( grid.cell_map[n] ) {
				(*pressure_vector)[n] = result->at(grid.cell_map[n]-1);
			} else {
				(*pressure_vector)[n] = 0.0;
			}
		}
	}
	//
	// Should be commneted out because it will be used in macoctreephelper functions
	// grid.clear_map();
}
//