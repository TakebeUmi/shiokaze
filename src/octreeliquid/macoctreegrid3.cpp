/*
**	macoctreegrid3.cpp
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
#include "macoctreegrid3.h"
#include <Eigen/Dense>
#include <shiokaze/utility/utility.h>
#include <shiokaze/array/shared_array3.h>
#include <shiokaze/array/array_interpolator3.h>
#include <shiokaze/core/console.h>
#include <shiokaze/core/timer.h>
#include "unstructured_extrapolator3.h"
#include "../../src/redistancer/unstructured_fastmarch3.h"
//
SHKZ_USING_NAMESPACE
using namespace macotreeliquid3_namespace;
//
void grid3::configure( configuration &config ) {
	//
	configuration::auto_group group(config,*this);
	config.get_unsigned("ChangeSparseGridResolution",param.change_sparse_array_resolution, "Resolution to change to tiled-based grid");
	config.get_string("SparseGrid",param.sparse_grid_module,"Sparse grid module name");
	config.get_bool("FirstOrder",param.first_order,"Whether to force 1st-order accuracy near T junctions");
	config.get_bool("ClampOrder",param.clamp_order,"Clamp order of accuracy");
	config.get_double("ClampFluidEps",param.clamp_fluid_eps,"Clamping eps for fluid faces");
	config.get_double("ClampSolidEps",param.clamp_solid_eps,"Clamping eps for solid faces");
	config.get_double("AccuracyEps",param.eps,"Second order accuracy clamp eps");
	config.get_double("PrecisionEps",param.precision_eps,"Precision eps");
	config.get_double("MLSEps",param.MLS_eps,"Eps for moving least square solve");
	config.get_bool("MLSConstDiag",param.MLS_constant_diag,"Use constant diagonal for MLS");
	config.get_bool("AccurateInterpolation",param.accurate_interpolation,"Accurate MLS interpolation");
	config.get_bool("UseInverse",param.use_inverse,"Use inverse for MLS solve");
	config.get_integer("AdaptivityType",param.adaptivity_type,"Adaptivity type");
	if( param.adaptivity_type == 2 || param.adaptivity_type == 3 ) {
		param.steep_adapvitiy = false;
	}
	config.get_bool("SteepAdaptivity",param.steep_adapvitiy,"Steep adaptivity");
	config.get_unsigned("DilateCount",param.dilate_count,"Dilation count for octree grading");
	config.get_double("PaddingForRemeshing",param.padding_for_remeshing,"Number of padding cells for remeshing");
	config.get_unsigned("SurfaceTensionSmoothCount",param.surftens_smooth_count,"Surface tension smoothing count");
	config.get_unsigned("PDEUpdateCount",param.pde_update_count,"PDE-based levelset redistancing count");
	config.get_bool("SimpleRedistance",param.simple_redistance,"Simple redistancing");
	config.get_bool("NormalizeLaplacian",param.normalize_laplacian,"Normalize laplacian" );
	config.get_unsigned("ErosionCount",param.erosion_count,"Eroson count for highresolution fluid reconstruct");
	config.get_bool("Debug",param.debug,"Debug mode");
}
//
void grid3::copy( const grid3 &grid ) {
	//
	levelset = grid.levelset;
	velocity = grid.velocity;
	area = grid.area;
	solid_cell = grid.solid_cell;
	ghost_cells = grid.ghost_cells;
	ghost_faces = grid.ghost_faces;
	cell_count = grid.cell_count;
	face_count = grid.face_count;
	ghost_cell_reference_count = grid.ghost_cell_reference_count;
	ghost_face_reference_count = grid.ghost_face_reference_count;
	flux_boundary_condition = grid.flux_boundary_condition;
	cell_map = grid.cell_map;
	face_map = grid.face_map;
	valid_cell_count = grid.valid_cell_count;
	valid_face_count = grid.valid_face_count;
	param = grid.param;
	//
	assert(layers.size()==grid.layers.size());
	for( unsigned depth=0; depth<layers.size(); ++depth ) {
		assert(layers[depth]->shape == grid.layers[depth]->shape );
		assert(layers[depth]->dx == grid.layers[depth]->dx );
	}
	//
	for( unsigned depth=0; depth<layers.size(); ++depth ) {
		//
		auto &layer_dst = *layers[depth];
		const auto &layer_src = *grid.layers[depth];
		//
		layer_dst.fill_flags.copy(layer_src.fill_flags);
		layer_dst.active_cells.copy(layer_src.active_cells);
		layer_dst.active_faces.copy(layer_src.active_faces);
		layer_dst.ghost_cell_indices.copy(layer_src.ghost_cell_indices);
		layer_dst.ghost_face_indices.copy(layer_src.ghost_face_indices);
	}
}
//
void grid3::clear() {
	//
	layers.clear();
	//
	levelset.clear();
	velocity.clear();
	area.clear();
	//
	levelset.shrink_to_fit();
	velocity.shrink_to_fit();
	area.shrink_to_fit();
	//
	cell_count = 0;
	face_count = 0;
	//
	clear_map();
}
//
void grid3::add_layer( const shape3 &shape, double dx ) {
	//
	auto layer = std::shared_ptr<layer3>(new layer3());
	layer->shape = shape;
	layer->dx = dx;
	//
	layer->set_environment("shape",&layer->shape);
	layer->set_environment("dx",&layer->dx);
	//
	std::string core_name ("lineararray3");
	if( shape.max() >= param.change_sparse_array_resolution ) {
		core_name = param.sparse_grid_module;
	}
	layer->fill_flags.set_core_name(core_name);
	layer->active_cells.set_core_name(core_name);
	layer->active_faces.set_core_name(core_name);
	//
	layer->setup_now();
	layers.push_back(layer);
}
//
void grid3::activate_cells( std::function<bool(char depth, const vec3d &p)> func ) {
	//
	scoped_timer timer(this);
	timer.tick(); console::dump( "Setting active fluid cells..." );
	//
	for( char depth=0; depth < layers.size(); ++depth ) {
		auto &layer = *layers[depth];
		layer.active_cells.clear();
		layer.fill_flags.clear();
	}
	//
	for( char depth=layers.size()-1; depth > 0; --depth ) {
		//
		auto &layer = *layers[depth];
		double dx = layer.dx;
		std::vector<vec3i> evaluation_points;
		//
		int next_depth = depth-1;
		auto &prev_layer = *layers[depth];
		auto &next_layer = *layers[next_depth];
		double prev_dx = prev_layer.dx;
		double next_dx = next_layer.dx;
		//
		if( depth == layers.size()-1 ) {
			prev_layer.active_cells.const_serial_all([&](int i, int j, int k, const auto &it) {
				evaluation_points.push_back(vec3i(i,j,k));
			});
		} else {
			prev_layer.active_cells.const_serial_actives([&](int i, int j, int k, const auto &it) {
				evaluation_points.push_back(vec3i(i,j,k));
			});
		}
		//
		parallel.for_each(evaluation_points.size(),[&]( size_t n ) {
			vec3d p = prev_dx*evaluation_points[n].cell();
			if( ! func(depth,p)) evaluation_points[n] = vec3i(-1,-1,-1);
		});
		//
		for( const auto &pi : evaluation_points ) {
			if( pi != vec3i(-1,-1,-1)) {
				for( int ii=0; ii<2; ++ii ) for( int jj=0; jj<2; ++jj ) for( int kk=0; kk<2; ++kk ) {
					next_layer.active_cells.set(2*pi+vec3i(ii,jj,kk),0);
				}
			}
		}
		//
		next_layer.active_cells.parallel_actives([&]( int i, int j, int k, auto &it ) {
			if( ! func(next_depth,next_dx*vec3i(i,j,k).cell())) it.set_off();
		});
		//
		next_layer.active_cells.dilate(1);
	}
	//
	for( char depth=0; depth < layers.size(); ++depth ) {
		auto &layer = *layers[depth];
		layer.active_cells.parallel_actives([&]( int i, int j, int k, auto &it ) {
			if( ! func(depth,layer.dx*vec3i(i,j,k).cell())) it.set_off();
		});
	}
	//
	if( param.steep_adapvitiy ) {
		for( char depth=1; depth < layers.size(); ++depth ) {
			layers[depth]->active_cells.clear();
		}
	}
	//
	size_t count (0);
	for( char depth=0; depth < layers.size(); ++depth ) {
		count += layers[depth]->active_cells.count();
	}
	console::dump( "Done. Count sum = %u. Took %s.\n", count, timer.stock("set_active").c_str());
	console::write("num_unbalanced_active_cells",count);
}
//
void grid3::activate_cells( std::function<double(const vec3d &p)> fluid, std::function<double(const vec3d &p)> solid ) {
	//
	auto still_active = [&]( const vec3d &p, char depth ) -> bool {
		if( param.adaptivity_type == 0 ) {
			return true;
		} else if( param.adaptivity_type == 1 ) {
			auto &layer = *layers[depth];
			const double d = (p-vec3d(0.5,p[1],0.5)).len();
			return std::max(0.0,d-0.2) < layer.dx;
		} else if( param.adaptivity_type == 2 ) {
			return 3*p[0] < depth+1;
		} else if( param.adaptivity_type == 3 ) {
			const double shift_x = (3+param.dilate_count) * layers[0]->dx;
			if( p[0] < 1/4.0+1/6.0-shift_x ) return depth >= 0;
			else if( p[0] < 1/4.0+2/6.0-shift_x ) return depth >= 1;
			else if( p[0] < 1/4.0+3/6.0-shift_x ) return depth >= 2;
			else return depth >= 3;
		} else {
			return true;
		}
	};
	//
	auto func = [&]( char depth, const vec3d &p ) -> bool {
		const double dx = layers[depth]->dx;
		const double d = std::max(0.0,std::abs(fluid(p))-param.padding_for_remeshing * dx);
		return d <= dx && still_active(p,depth);
	};
	//
	activate_cells(func);
}
//
void grid3::assign_levelset( std::function<double( const vec3d &p )> fluid, std::function<double( const vec3d &p )> solid ) {
	//
	scoped_timer timer(this);
	timer.tick(); console::dump( "Assigning level set..." );
	//
	// Assigning fluid levelset
	iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
		levelset[cell_id.index] = fluid(get_cell_position(cell_id));
	});
	//
	if( solid ) {
		for( size_t depth=0; depth < layers.size(); ++depth ) {
			//
			auto &layer = *layers[depth];
			double dx = layer.dx;
			//
			shared_array3<Real> solid_array(layer.shape.nodal());
			shared_macarray3<Real> areas(layer.shape);
			//
			// Compute nodal solid level set
			layer.active_cells.const_serial_actives([&]( int i, int j, int k, const auto &it ) {
				for( int ii=0; ii<2; ++ii ) for( int jj=0; jj<2; ++jj ) for( int kk=0; kk<2; ++kk ) {
					solid_array->set(i+ii,j+jj,k+kk,0.0);
				}
			});
			//
			solid_array->parallel_actives([&]( int i, int j, int k, auto &it ) {
				it.set(solid(dx*vec3i(i,j,k).nodal()));
			});
			//
			// Compute area fraction on each layer
			layer.macutility->compute_area_fraction(solid_array(),areas());
			layer.active_faces.const_parallel_actives([&]( char dim, int i, int j, int k, const auto &it, int tid ) {
				area[it()] = areas()[dim](i,j,k);
			});
		}
		//
		std::vector<char> solid_cell_tmp(cell_count);
		iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
			bool has_solid (false);
			for( int dim : DIMS3 ) {
				iterate_face_neighbors(cell_id,dim,[&]( const face_id3 &face_id ){
					if( area[face_id.index] < 1.0 ) {
						has_solid = true;
					}
				});
			}
			solid_cell_tmp[cell_id.index] = has_solid;
		});
		for( uint_type n=0; n<cell_count; ++n ) {
			solid_cell[n] = solid_cell_tmp[n];
		}
	} else {
		std::fill(area.begin(),area.end(),1.0);
	}
	console::dump( "Done. Took %s.\n", timer.stock("assign_levelset").c_str());
}
//
void grid3::set_flux_boundary_condition( const flux_boundary_condition3 &boundary_cond ) {
	//
	flux_boundary_condition = boundary_cond;
}
//
void grid3::set_velocity( std::function<double( const vec3d &p, char dim )> func ) {
	//
	scoped_timer timer(this);
	timer.tick(); console::dump( "Set velocity..." );
	//
	iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
		velocity[face_id.index] = func(get_face_position(face_id),face_id.dim);
	});
	//
	console::dump( "Done. Took %s.\n", timer.stock("set_velocity").c_str());
}
//
void grid3::balance_layers() {
	//
	scoped_timer timer(this);
	timer.tick(); console::dump( "Balancing layers..." );
	//
	layers[0]->active_cells.dilate(param.dilate_count);
	//
	for( unsigned depth=0; depth < layers.size()-1; ++depth ) {
		//
		const auto &coarse_layer = layers[depth+1];
		auto &high_layer = layers[depth];
		//
		const unsigned factor_x = high_layer->shape[0] / coarse_layer->shape[0];
		const unsigned factor_y = high_layer->shape[1] / coarse_layer->shape[1];
		const unsigned factor_z = high_layer->shape[2] / coarse_layer->shape[2];
		//
		high_layer->active_cells.serial_actives([&]( int i, int j, int k, auto &it ) {
			coarse_layer->active_cells.set(vec3i(i/factor_x,j/factor_y,k/factor_z),0);
		});
		//
		coarse_layer->active_cells.dilate(param.dilate_count);
	}
	//
	for( unsigned depth=0; depth < layers.size()-1; ++depth ) {
		//
		const auto &coarse_layer = layers[depth+1];
		auto &high_layer = layers[depth];
		//
		const unsigned factor_x = high_layer->shape[0] / coarse_layer->shape[0];
		const unsigned factor_y = high_layer->shape[1] / coarse_layer->shape[1];
		const unsigned factor_z = high_layer->shape[2] / coarse_layer->shape[2];
		//
		coarse_layer->active_cells.serial_actives([&]( int i, int j, int k, auto &it ) {
			//
			bool found (false);
			for( int fk=factor_z*k; fk<factor_z*(k+1); ++fk ) for( int fj=factor_y*j; fj<factor_y*(j+1); ++fj ) for( int fi=factor_x*i; fi<factor_x*(i+1); ++fi ) {
				if(high_layer->active_cells.active(fi,fj,fk)) {
					found = true;
					break;
				}
			}
			if( found ) {
				for( int fk=factor_z*k; fk<factor_z*(k+1); ++fk ) for( int fj=factor_y*j; fj<factor_y*(j+1); ++fj ) for( int fi=factor_x*i; fi<factor_x*(i+1); ++fi ) {
					high_layer->active_cells.set(fi,fj,fk,0);
				}
				coarse_layer->fill_flags.set(i,j,k);
			}
		});
		//
	}
	//
	for( auto layer_ptr : layers ) {
		layer_ptr->active_cells.parallel_actives([&]( int i, int j, int k, auto &it ) {
			if( layer_ptr->fill_flags(i,j,k)) it.set_off();
		});
	}
	//
	layers.back()->active_cells.activate_all();
	layers.back()->active_cells.parallel_actives([&]( int i, int j, int k, auto &it ) {
		if( layers.back()->fill_flags(i,j,k)) it.set_off();
	});
	//
	console::dump( "Done. Took %s.\n", timer.stock("balance_layers").c_str());
}
//
void grid3::assign_indices() {
	//
	scoped_timer timer(this);
	timer.tick(); console::dump( "Assigning indices..." );
	//
	size_t pressure_index (0);
	size_t velocity_index (0);
	//
	for( char depth=layers.size()-1; depth >= 0; --depth ) {
		//
		auto &layer = *layers[depth];
		layer.active_faces.clear();
		layer.active_cells.serial_actives([&]( int i, int j, int k, auto &it ) {
			//
			it.set(pressure_index++);
			//
			for( char dim : DIMS3 ) for( int dir=-1; dir<=1; dir+=2 ) {
				vec3i pj = vec3i(i+dir*(dim==0),j+dir*(dim==1),k+dir*(dim==2));
				if( ! layer.shape.out_of_bounds(pj)) {
					if( layer.active_cells.active(pj) || layer.fill_flags(pj)) {
						int offset = dir == -1 ? 0 : 1;
						layer.active_faces[dim].set(i+offset*(dim==0),j+offset*(dim==1),k+offset*(dim==2),0);
					}
				}
			}
		});
		//
		layer.active_faces.serial_actives([&]( char dim, int i, int j, int k, auto &it ) {
			it.set(velocity_index++);
		});
	}
	//
	assert( pressure_index <= std::numeric_limits<uint_type>::max() );
	assert( velocity_index <= std::numeric_limits<uint_type>::max() );
	//
	levelset.clear();
	velocity.clear();
	area.clear();
	solid_cell.clear();
	//
	levelset.resize(pressure_index);
	velocity.resize(velocity_index);
	area.resize(velocity_index);
	solid_cell.resize(pressure_index);
	//
	levelset.shrink_to_fit();
	velocity.shrink_to_fit();
	area.shrink_to_fit();
	solid_cell.shrink_to_fit();
	//
	cell_count = pressure_index;
	face_count = velocity_index;
	//
	console::dump( "Done. Num cells = %u, Num face = %u. Took %s.\n", pressure_index, velocity_index, timer.stock("assign_indices").c_str());
	console::write("num_active_cells",pressure_index);
	console::write("num_active_faces",velocity_index);
	//
	if( param.accurate_interpolation ) {
		//
		/********************************** Ghost cell preprocessing **********************************/
		//
		size_t ghost_cell_reference_index (0);
		//
		for( char depth=1; depth < layers.size(); depth++ ) {
			auto &fine_layer = *layers[depth-1];
			const auto &coarse_layer = *layers[depth];
			fine_layer.ghost_cell_indices.clear();
			coarse_layer.active_faces.const_serial_actives([&]( int dim, int i, int j, int k, const auto &it ) {
				for( int dir=-1; dir<=0; dir++ ) {
					const vec3i neighbor_pi = vec3i(i,j,k)+dir*vec3i(dim==0,dim==1,dim==2);
					if( coarse_layer.fill_flags(neighbor_pi)) {
						vec3i ivec, jvec;
						if( dim == 0 ) {
							ivec = vec3i(0,1,0);
							jvec = vec3i(0,0,1);
						} else if( dim == 1 ) {
							ivec = vec3i(1,0,0);
							jvec = vec3i(0,0,1);
						} else if( dim == 2 ) {
							ivec = vec3i(1,0,0);
							jvec = vec3i(0,1,0);
						}
						const vec3i nvec = vec3i(dim==0,dim==1,dim==2);
						const vec3i bottom = 2*vec3i(i,j,k);
						vec3i small_cells[4] = { bottom, bottom+ivec, bottom+ivec+jvec, bottom+jvec };
						for( auto &pi : small_cells ) {
							if( dir == -1 ) {
								fine_layer.ghost_cell_indices.set(pi,0);
							} else {
								pi -= nvec;
								fine_layer.ghost_cell_indices.set(pi,0);
							}
						}
					}
				}
			});
			//
			fine_layer.ghost_cell_indices.serial_actives([&]( int i, int j, int k, auto &it ) {
				it.set(ghost_cell_reference_index++);
			});
		}
		//
		ghost_cell_reference_count = ghost_cell_reference_index;
		ghost_cells.clear();
		ghost_cells.resize(ghost_cell_reference_count);
		ghost_cells.shrink_to_fit();
		//
		for( char depth=0; depth < layers.size()-1; depth++ ) {
			auto &layer = *layers[depth];
			const auto &coarse_layer = *layers[depth+1];
			layer.ghost_cell_indices.const_parallel_actives([&]( int i, int j, int k, const auto &it ) {
				bool valid (false);
				ghost_cell3 ghost_cell;
				for( int dim : DIMS3 ) {
					vec3i neighbor_pi;
					int dir;
					for( dir=-1; dir<=1; dir+=2 ) {
						neighbor_pi = vec3i(i,j,k)+dir*vec3i(dim==0,dim==1,dim==2);
						if( layer.active_cells.safe_active(neighbor_pi)) {
							//
							// Found active cell, look for coarse cell with the opposite
							interpolation_data3 interpolation_data;
							vec3d coarse_p = (vec3i(i,j,k).cell()-0.5*dir*vec3d(dim==0,dim==1,dim==2)) / 2.0;
							vec3i indices[8];
							double coef[8];
							array_interpolator3::interpolate_coef(coarse_layer.shape,coarse_p-vec3d(0.5,0.5,0.5),indices,coef);
							int slot (0);
							double sum (0.0);
							for( int n=0; n<8; ++n ) {
								if( coef[n] > param.precision_eps && coarse_layer.active_cells.active(indices[n])) {
									assert(slot<4);
									sum += coef[n];
									interpolation_data.coef[slot] = coef[n];
									interpolation_data.indices[slot] = coarse_layer.active_cells(indices[n]);
									slot ++;
								}
							}
							if( slot && std::abs(sum-1.0) < param.precision_eps ) {
								interpolation_data.count = slot;
								interpolation_data.p = coarse_layer.dx*coarse_p;
								ghost_cell.data.push_back(interpolation_data);
								valid = true;
								break;
							}
						}
					}
					if( valid ) {
						interpolation_data3 interpolation_data;
						interpolation_data.count = 1;
						interpolation_data.p = layer.dx*neighbor_pi.cell();
						interpolation_data.coef[0] = 1.0;
						interpolation_data.indices[0] = layer.active_cells(neighbor_pi);
						ghost_cell.data.push_back(interpolation_data);
						ghost_cell.p = layer.dx*(vec3i(i,j,k).cell()+0.5*dir*vec3d(dim==0,dim==1,dim==2));
						//
						double sum (0.0);
						assert(ghost_cell.data.size()==2);
						for( int n=0; n<ghost_cell.data.size(); ++n ) {
							double scale;
							if( n == 0 ) scale = 0.5/1.5;
							else if( n == 1 ) scale = 1.0/1.5;
							const auto &d = ghost_cell.data[n];
							for( int m=0; m<d.count; ++m ) {
								const double w = scale*d.coef[m];
								sum += w;
								ghost_cell.combined.push_back(std::make_pair(d.indices[m],w));
							}
						}
						if( std::abs(sum-1.0) > param.precision_eps ) {
							console::dump( "C0: sum = %e\n", sum );
							exit(0);
						}
						if( dim == 0 ) {
							if( dir == -1 ) ghost_cell.shrink_info |= shrink_left;
							if( dir == 1 ) ghost_cell.shrink_info |= shrink_right;
						} else if( dim == 1 ) {
							if( dir == -1 ) ghost_cell.shrink_info |= shrink_bottom;
							if( dir == 1 ) ghost_cell.shrink_info |= shrink_top;
						} else if( dim == 2 ) {
							if( dir == -1 ) ghost_cell.shrink_info |= shrink_front;
							if( dir == 1 ) ghost_cell.shrink_info |= shrink_back;
						}
						break;
					}
				}
				if( valid ) {
					ghost_cells[it()].push_back(ghost_cell);
				} else {
					// Diagonal test
					ghost_cell3 ghost_cell;
					vec3i shift_dir;
					int num_adjacent_count (0);
					for( int dim : DIMS3 ) {
						for( int dir=-1; dir<=1; dir+=2 ) {
							const vec3i neighbor_pi = vec3i(i,j,k)+dir*vec3i(dim==0,dim==1,dim==2);
							if( layer.active_cells.safe_active(neighbor_pi)) {
								shift_dir[dim] += dir;
								num_adjacent_count ++;
								interpolation_data3 interpolation_data;
								interpolation_data.count = 1;
								interpolation_data.p = layer.dx*neighbor_pi.cell();
								interpolation_data.coef[0] = 1.0;
								interpolation_data.indices[0] = layer.active_cells(neighbor_pi);
								ghost_cell.data.push_back(interpolation_data);
							}
						}
					}
					if( num_adjacent_count >= 2 ) {
						interpolation_data3 interpolation_data;
						vec3d coarse_p = (vec3i(i,j,k).cell()-0.5*vec3d(shift_dir)) / 2.0;
						vec3i indices[8];
						double coef[8];
						array_interpolator3::interpolate_coef(coarse_layer.shape,coarse_p-vec3d(0.5,0.5,0.5),indices,coef);
						int slot (0);
						for( int n=0; n<8; ++n ) {
							if( coef[n] > param.precision_eps && coarse_layer.active_cells.active(indices[n])) {
								assert(slot<2);
								interpolation_data.coef[slot] = coef[n];
								interpolation_data.indices[slot] = coarse_layer.active_cells(indices[n]);
								slot ++;
							}
						}
						if( slot && slot <= 2 ) {
							interpolation_data.count = slot;
							interpolation_data.p = coarse_layer.dx*coarse_p;
							ghost_cell.data.push_back(interpolation_data);
							valid = true;
						}
					}
					if( valid ) {
						//
						// Split
						for( int dim : DIMS3 ) {
							const vec3i neighbor_pi = vec3i(i+(dim==0)*shift_dir[0],j+(dim==1)*shift_dir[1],k+(dim==2)*shift_dir[2]);
							if( layer.active_cells.active(neighbor_pi)) {
								//
								ghost_cell3 ghost_cell_split(ghost_cell);
								ghost_cell_split.p = layer.dx*vec3d(i+0.5*(dim==0)*shift_dir[0],j+0.5*(dim==1)*shift_dir[1],k+0.5*(dim==2)*shift_dir[2]).cell();
								interpolation_data3 interpolation_data;
								interpolation_data.count = 1;
								interpolation_data.p = layer.dx*neighbor_pi.cell();
								interpolation_data.coef[0] = 1.0;
								interpolation_data.indices[0] = layer.active_cells(neighbor_pi);
								ghost_cell_split.data.push_back(interpolation_data);
								//
								double sum (0.0);
								assert(ghost_cell_split.data.size()==num_adjacent_count+2);
								for( int n=0; n<ghost_cell_split.data.size(); ++n ) {
									double scale;
									if( n < num_adjacent_count ) scale = 0.25 / num_adjacent_count;
									else if( n == num_adjacent_count ) scale = 0.25;
									else scale = 0.5;
									const auto &d = ghost_cell_split.data[n];
									double sum_check (0.0);
									for( int m=0; m<d.count; ++m ) {
										sum_check += d.coef[m];
										const double w = scale*d.coef[m];
										sum += w;
										ghost_cell_split.combined.push_back(std::make_pair(d.indices[m],w));
									}
									assert( std::abs(sum_check-1.0) < param.precision_eps );
								}
								if( std::abs(sum-1.0) > param.precision_eps ) {
									console::dump( "C1: sum = %e\n", sum );
									exit(0);
								}
								ghost_cell_split.shrink_info = 0;
								if( dim == 0 ) {
									if( shift_dir[0] == -1 ) ghost_cell_split.shrink_info |= shrink_left;
									if( shift_dir[0] == 1 ) ghost_cell_split.shrink_info |= shrink_right;
								} else if( dim == 1 ) {
									if( shift_dir[1] == -1 ) ghost_cell_split.shrink_info |= shrink_bottom;
									if( shift_dir[1] == 1 ) ghost_cell_split.shrink_info |= shrink_top;
								} else if( dim == 2 ) {
									if( shift_dir[2] == -1 ) ghost_cell_split.shrink_info |= shrink_front;
									if( shift_dir[2] == 1 ) ghost_cell_split.shrink_info |= shrink_back;
								} 
								ghost_cells[it()].push_back(ghost_cell_split);
							}
						}
					}
				}
			});
		}
		//
		/********************************** Ghost face preprocessing **********************************/
		//
		size_t ghost_face_reference_index (0);
		//
		for( char depth=1; depth < layers.size(); depth++ ) {
			auto &fine_layer = *layers[depth-1];
			const auto &coarse_layer = *layers[depth];
			fine_layer.ghost_face_indices.clear();
			coarse_layer.active_faces.const_serial_actives([&]( int dim, int i, int j, int k, const auto &it ) {
				for( int dir=-1; dir<=0; dir++ ) {
					const vec3i neighbor_pi = vec3i(i,j,k)+dir*vec3i(dim==0,dim==1,dim==2);
					if( coarse_layer.fill_flags(neighbor_pi)) {
						vec3i ivec, jvec;
						int ivec_dim;
						int jvec_dim;
						if( dim == 0 ) {
							ivec = vec3i(0,1,0); ivec_dim=1;
							jvec = vec3i(0,0,1); jvec_dim=2;
						} else if( dim == 1 ) {
							ivec = vec3i(1,0,0); ivec_dim=0;
							jvec = vec3i(0,0,1); jvec_dim=2;
						} else if( dim == 2 ) {
							ivec = vec3i(1,0,0); ivec_dim=0;
							jvec = vec3i(0,1,0); jvec_dim=1;
						}
						const vec3i nvec = vec3i(dim==0,dim==1,dim==2);
						const vec3i bottom = 2*vec3i(i,j,k);
						for( int dim2 : DIMS3 ) {
							if( dim == dim2 ) continue;
							int ii_len = ivec_dim == dim2 ? 3 : 2;
							int jj_len = jvec_dim == dim2 ? 3 : 2;
							for( int ii=0; ii<ii_len; ++ii ) for( int jj=0; jj<jj_len; ++jj ) {
								vec3i pi = bottom+ii*ivec+jj*jvec;
								if( dir == 0 ) {
									pi -= nvec;
									fine_layer.ghost_face_indices[dim2].set(pi,0);
								} else {
									fine_layer.ghost_face_indices[dim2].set(pi,0);
								}
							}
						}
					}
				}
			});
			//
			fine_layer.ghost_face_indices.serial_actives([&]( int dim, int i, int j, int k, auto &it ) {
				it.set(ghost_face_reference_index++);
			});
		}
		//
		ghost_face_reference_count = ghost_face_reference_index;
		ghost_faces.clear();
		ghost_faces.resize(ghost_face_reference_count);
		ghost_faces.shrink_to_fit();
		//
		auto fix_interpolation = []( int dim, const shape3 &shape, vec3i indices[8] ) {
			for( int n=0; n<8; ++n ) {
				if( indices[n][dim] == 0 ) for( int m=0; m<8; ++m ) {
					if( indices[n]+vec3i(dim==0,dim==1,dim==2) == indices[m] ) {
						indices[n] = indices[m];
						break;
					}
				}
				if( indices[n][dim] == shape[dim]-1 ) for( int m=0; m<8; ++m ) {
					if( indices[n]-vec3i(dim==0,dim==1,dim==2) == indices[m] ) {
						indices[n] = indices[m];
						break;
					}
				}
			}
		};
		//
		auto check_interpolation = [&]( const interpolation_data3 &data ) {
			double sum (0.0);
			for( int m=0; m<data.count; ++m ) {
				sum += data.coef[m];
			}
			return std::abs(sum-1.0) < param.precision_eps;
		};
		//
		auto combine_with_scale = [&]( ghost_face3 &ghost_face, const std::vector<double> &scales ) {
			double sum (0.0);
			assert(scales.size() == ghost_face.data.size());
			for( int n=0; n<ghost_face.data.size(); ++n ) {
				const double &scale = scales[n];
				const auto &d = ghost_face.data[n];
				assert(check_interpolation(d));
				double sub_check (0.0);
				for( int m=0; m<d.count; ++m ) {
					const double w = scale*d.coef[m];
					sub_check += d.coef[m];
					sum += w;
					ghost_face.combined.push_back(std::make_pair(d.indices[m],w));
				}
				assert(std::abs(sub_check-1.0) < param.precision_eps);
			}
			assert(std::abs(sum-1.0) < param.precision_eps);
		};
		//
		for( char depth=0; depth < layers.size()-1; depth++ ) {
			auto &layer = *layers[depth];
			const auto &coarse_layer = *layers[depth+1];
			layer.ghost_face_indices.const_parallel_actives([&]( int dim, int i, int j, int k, const auto &it ) {
				//
				vec3i edge_shifted_pi(i,j,k);
				if( edge_shifted_pi[dim] == 0 ) edge_shifted_pi[dim] ++;
				if( edge_shifted_pi[dim] == layer.shape[dim] ) edge_shifted_pi[dim] --;
				//
				const vec3d &face_position = layer.dx * vec3i(i,j,k).face(dim);
				const vec3d &edge_shiftted_face_position = layer.dx * edge_shifted_pi.face(dim);
				const auto &coarse_layer = *layers[depth+1];
				const vec3i cell_pi = coarse_layer.shape.find_cell(face_position/coarse_layer.dx);
				const vec3d &cell_position = coarse_layer.dx * cell_pi.cell();
				//
				vec3i shift_dir;
				for( int dim2 : DIMS3 ) if( dim != dim2 ) {
					int dir = cell_position[dim2] < face_position[dim2] ? 1 : -1;
					const vec3i &neighbor_pi = edge_shifted_pi+dir*vec3i(dim2==0,dim2==1,dim2==2);
					if( layer.active_faces[dim].safe_active(neighbor_pi)) {
						shift_dir[dim2] = dir;
					}
				}
				//
				for( int dim2 : DIMS3 ) if( dim2 != dim ) if( shift_dir[dim2] ) {
					//
					ghost_face3 ghost_face;
					int sgn = shift_dir[dim2];
					const vec3d &half_shift = 0.5*sgn*layer.dx*vec3d(dim2==0,dim2==1,dim2==2);
					ghost_face.p = face_position+half_shift;
					//
					if( dim2 == 0 ) {
						if( sgn == 1 ) ghost_face.shrink_info |= shrink_right;
						if( sgn == -1 ) ghost_face.shrink_info |= shrink_left;
					} else if( dim2 == 1 ) {
						if( sgn == 1 ) ghost_face.shrink_info |= shrink_top;
						if( sgn == -1 ) ghost_face.shrink_info |= shrink_bottom;
					} else if( dim2 == 2 ) {
						if( sgn == 1 ) ghost_face.shrink_info |= shrink_back;
						if( sgn == -1 ) ghost_face.shrink_info |= shrink_front;
					}
					//
					interpolation_data3 neighbor_data;
					//
					const vec3i &neighbor_pi = edge_shifted_pi+sgn*vec3i(dim2==0,dim2==1,dim2==2);
					if( layer.active_faces[dim].safe_active(neighbor_pi)) {
						neighbor_data.count = 1;
						neighbor_data.p = layer.dx * neighbor_pi.face(dim);
						neighbor_data.coef[0] = 1.0;
						neighbor_data.indices[0] = layer.active_faces[dim](neighbor_pi);
					}
					//
					auto interpolate_coarse = [&]( const vec3d &p ) {
						//
						interpolation_data3 result, interpolation_data;
						vec3i indices[8];
						double coef[8];
						const vec3d pos = p/coarse_layer.dx-0.5*vec3d(dim!=0,dim!=1,dim!=2);
						array_interpolator3::interpolate_coef(coarse_layer.shape.face(dim),pos,indices,coef);
						fix_interpolation(dim,coarse_layer.shape.face(dim),indices);
						int slot (0);
						double sum (0.0);
						for( int n=0; n<8; ++n ) {
							if( coef[n] > param.precision_eps && coarse_layer.active_faces[dim].active(indices[n])) {
								assert(slot<4);
								sum += coef[n];
								interpolation_data.coef[slot] = coef[n];
								interpolation_data.indices[slot] = coarse_layer.active_faces[dim](indices[n]);
								slot ++;
							}
						}
						if( std::abs(sum-1.0) < param.precision_eps ) {
							interpolation_data.count = slot;
							interpolation_data.p = p;
							result = interpolation_data;
						}
						return result;
					};
					//
					if( neighbor_data.count ) {
						interpolation_data3 coarse_data = interpolate_coarse(edge_shiftted_face_position-half_shift);
						if( coarse_data.count ) {
							assert(check_interpolation(neighbor_data));
							assert(check_interpolation(coarse_data));
							ghost_face.data.push_back(neighbor_data);
							ghost_face.data.push_back(coarse_data);
							combine_with_scale(ghost_face,{1.0/1.5,0.5/1.5});
							ghost_faces[it()].push_back(ghost_face);
						} else {
							//
							// Split
							interpolation_data3 split_data;
							int slot (0);
							for( int dim3 : DIMS3 ) if( dim != dim3 ) {
								const vec3i &neighbor_pi = edge_shifted_pi+vec3i((dim3==0)*shift_dir[0],(dim3==1)*shift_dir[1],(dim3==2)*shift_dir[2]);
								split_data.indices[slot] = layer.active_faces[dim](neighbor_pi);
								split_data.coef[slot] = 0.5;
								slot ++;
							}
							split_data.count = slot;
							assert(slot==2);
							//
							coarse_data = interpolate_coarse(face_position-0.5*layer.dx*vec3d(shift_dir));
							assert(coarse_data.count);
							assert(check_interpolation(coarse_data));
							assert(check_interpolation(split_data));
							assert(check_interpolation(coarse_data));
							assert(check_interpolation(neighbor_data));
							//
							ghost_face.data.push_back(split_data);
							ghost_face.data.push_back(coarse_data);
							ghost_face.data.push_back(neighbor_data);
							//
							combine_with_scale(ghost_face,{0.25,0.25,0.5});
							ghost_faces[it()].push_back(ghost_face);
						}
					}
				}
			});
		}
	}
}
//
void grid3::iterate_cell_neighbors( const cell_id3 &cell_id, std::function<void( char dim, const cell_id3 &cell_id )> func ) const {
	//
	const auto &layer = *layers[cell_id.depth];
	for( char dim : DIMS3 ) for( int dir=-1; dir<=1; dir+=2 ) {
		//
		vec3i ivec, jvec;
		if( dim == 0 ) {
			ivec = vec3i(0,1,0);
			jvec = vec3i(0,0,1);
		} else if( dim == 1 ) {
			ivec = vec3i(1,0,0);
			jvec = vec3i(0,0,1);
		} else if( dim == 2 ) {
			ivec = vec3i(1,0,0);
			jvec = vec3i(0,1,0);
		}
		vec3i nvec = vec3i(dim==0,dim==1,dim==2);
		vec3i pi = cell_id.pi+dir*nvec;
		bool found (false);
		//
		if( cell_id.depth > 0 ) {
			const auto &high_layer = *layers[cell_id.depth-1];
			vec3i bottom = 2*pi;
			vec3i small_cells[4] = { bottom, bottom+ivec, bottom+jvec, bottom+ivec+jvec };
			if( dir == -1 ) for( auto &pi : small_cells ) pi += nvec;
			for( int n=0; n<4; ++n ) {
				if( high_layer.active_cells.safe_active(small_cells[n])) {
					uint_type column = high_layer.active_cells(small_cells[n]);
					func(dim,{(char)(cell_id.depth-1),small_cells[n],column});
					found = true;
				}
			}
		}
		if( ! found ) {
			if( layer.active_cells.safe_active(pi)) {
				uint_type column = layer.active_cells(pi);
				func(dim,{cell_id.depth,pi,column});
				found = true;
			}
		}
		if( ! found && cell_id.depth < layers.size()-1 ) {
			const auto &coarse_layer = *layers[cell_id.depth+1];
			vec3i cell = 0.5 * pi;
			if( coarse_layer.active_cells.safe_active(cell)) {
				uint_type column = coarse_layer.active_cells(cell);
				func(dim,{(char)(cell_id.depth+1),cell,column});
				found = true;
			}
		}
	}
}
//
void grid3::iterate_face_neighbors( const face_id3 &face_id, std::function<void( const face_id3 &face_id )> func ) const {
	//
	const auto &layer = *layers[face_id.depth];
	//
	vec3i ivec, jvec;
	if( face_id.dim == 0 ) {
		ivec = vec3i(0,1,0);
		jvec = vec3i(0,0,1);
	} else if( face_id.dim == 1 ) {
		ivec = vec3i(1,0,0);
		jvec = vec3i(0,0,1);
	} else if( face_id.dim == 2 ) {
		ivec = vec3i(1,0,0);
		jvec = vec3i(0,1,0);
	}
	vec3i nvec = vec3i(face_id.dim==0,face_id.dim==1,face_id.dim==2);
	//
	for( int dir=-1; dir<=1; dir+=2 ) {
		bool found (false);
		if( face_id.depth > 0 ) {
			const auto &high_layer = *layers[face_id.depth-1];
			vec3i bottom = 2*face_id.pi+dir*nvec;
			vec3i small_faces[4] = { bottom, bottom+ivec, bottom+jvec, bottom+ivec+jvec };
			for( int n=0; n<4; ++n ) {
				if( high_layer.active_faces[face_id.dim].safe_active(small_faces[n])) {
					uint_type column = high_layer.active_faces[face_id.dim](small_faces[n]);
					func({(char)(face_id.depth-1),face_id.dim,small_faces[n],column});
					found = true;
				}
			}
		}
		if( ! found ) {
			vec3i pi = face_id.pi+dir*nvec;
			if( layer.active_faces[face_id.dim].safe_active(pi)) {
				uint_type column = layer.active_faces[face_id.dim](pi);
				func({face_id.depth,face_id.dim,pi,column});
				found = true;
			}
		}
		if( ! found && face_id.depth < layers.size()-1 ) {
			const auto &coarse_layer = *layers[face_id.depth+1];
			vec3i cell = 0.5*face_id.pi+(dir+1)/2*nvec;
			if( coarse_layer.active_faces[face_id.dim].safe_active(cell)) {
				uint_type column = coarse_layer.active_faces[face_id.dim](cell);
				func({(char)(face_id.depth+1),face_id.dim,cell,column});
				found = true;
			}
		}
	}
	//
	for( int dir=-1; dir<=1; dir+=2 ) for( int n=0; n<2; ++n ) {
		//
		vec3i tvec = (n == 0) ? ivec : jvec;
		vec3i ttvec = (n == 0) ? jvec : ivec;
		bool found (false);
		//
		if( face_id.depth > 0 ) {
			const auto &high_layer = *layers[face_id.depth-1];
			for( int m=0; m<2; ++m ) {
				vec3i face = 2*face_id.pi+(dir==-1 ? -tvec : 2*tvec)+m*ttvec;
				if( high_layer.active_faces[face_id.dim].safe_active(face)) {
					uint_type column = high_layer.active_faces[face_id.dim](face);
					func({(char)(face_id.depth-1),face_id.dim,face,column});
					found = true;
				}
			}
		}
		if( ! found ) {
			vec3i pi = face_id.pi+dir*tvec;
			if( layer.active_faces[face_id.dim].safe_active(pi)) {
				uint_type column = layer.active_faces[face_id.dim](pi);
				func({face_id.depth,face_id.dim,pi,column});
				found = true;
			}
		}
		if( ! found ) {
			if( face_id.depth < layers.size()-1 ) {
				const auto &coarse_layer = *layers[face_id.depth+1];
				vec3i face_back = 0.5*face_id.pi+dir*tvec;
				vec3i coarse_faces[2] = { face_back, face_back+nvec };
				for( int m=0; m<2; ++m ) {
					if( face_id.pi[face_id.dim] % 2 == 0 && m==1 ) continue;
					if( coarse_layer.active_faces[face_id.dim].safe_active(coarse_faces[m])) {
						uint_type column = coarse_layer.active_faces[face_id.dim](coarse_faces[m]);
						func({(char)(face_id.depth+1),face_id.dim,coarse_faces[m],column});
						found = true;
					}
				}
			}
		}
	}
}
//
void grid3::iterate_face_neighbors( const cell_id3 &cell_id, char dim, std::function<void( const face_id3 &face_id )> func ) const {
	//
	const auto &layer = *layers[cell_id.depth];
	assert(layer.active_cells.active(cell_id.pi));
	//
	for( int dir=0; dir<=1; ++dir ) {
		vec3i fpi = cell_id.pi+dir*vec3i(dim==0,dim==1,dim==2);
		if( layer.active_faces[dim].active(fpi)) {
			func({cell_id.depth,dim,fpi,layer.active_faces[dim](fpi)});
		} else if( cell_id.depth < layers.size()-1 ) {
			const auto &coarse_layer = *layers[cell_id.depth+1];
			vec3i fpi = cell_id.pi/2+dir*vec3i(dim==0,dim==1,dim==2);
			if(coarse_layer.active_faces[dim].active(fpi)) {
				func({(char)(cell_id.depth+1),dim,fpi,coarse_layer.active_faces[dim](fpi)});
			}
		}
	}
}
//
void grid3::iterate_active_cells( const std::function<void( const cell_id3 &cell_id, int thread_index )> func ) const {
	//
	for( char depth=0; depth<layers.size(); ++depth ) {
		const auto &layer = *layers[depth];
		layer.active_cells.const_parallel_actives([&]( int i, int j, int k, const auto &it, int thread_index ) {
			func({depth,vec3i(i,j,k),it()},thread_index);
		});
	}
}
//
void grid3::iterate_active_faces( const std::function<void( const face_id3 &face_id, int thread_index )> func ) const {
	//
	for( char depth=0; depth<layers.size(); ++depth ) {
		const auto &layer = *layers[depth];
		layer.active_faces.const_parallel_actives([&]( char dim, int i, int j, int k, const auto &it, int thread_index ) {
			func({depth,dim,vec3i(i,j,k),it()},thread_index);
		});
	}
}
//
void grid3::serial_iterate_active_cells( const std::function<void( const cell_id3 &cell_id )> func ) const {
	//
	for( char depth=0; depth<layers.size(); ++depth ) {
		const auto &layer = *layers[depth];
		layer.active_cells.const_serial_actives([&]( int i, int j, int k, const auto &it ) {
			func({depth,vec3i(i,j,k),it()});
		});
	}
}
//
void grid3::serial_iterate_active_faces( const std::function<void( const face_id3 &face_id )> func ) const {
	//
	for( char depth=0; depth<layers.size(); ++depth ) {
		const auto &layer = *layers[depth];
		layer.active_faces.const_serial_actives([&]( char dim, int i, int j, int k, const auto &it ) {
			func({depth,dim,vec3i(i,j,k),it()});
		});
	}
}
//
void grid3::get_gradient( const face_id3 &face_id, std::function<void( const cell_id3 &cell_id, double value, const gradient_info3 &info )> func ) const {
	//
	const auto &layer = *layers[face_id.depth];
	const vec3i p_forward = face_id.pi;
	const vec3i p_backward = face_id.pi-vec3i(face_id.dim==0,face_id.dim==1,face_id.dim==2);
	const Real area = this->area[face_id.index];
	//
	const Real dx = layer.dx;
	const double scale = 1.0 / dx;
	//
	if( layer.active_cells.active(p_forward) && layer.active_cells.active(p_backward)) {
		//
		const uint_type foward_index = layer.active_cells(p_forward);
		const uint_type backward_index = layer.active_cells(p_backward);
		//
		const cell_id3 &forward_id = {face_id.depth,p_forward,foward_index};
		const cell_id3 &backward_id = {face_id.depth,p_backward,backward_index};
		//
		const Real levelsets[] = { levelset[foward_index], levelset[backward_index] };
		const Real rho = utility::fraction(levelsets[0],levelsets[1]);
		const bool cross_interface = rho > 0.0 && rho < 1.0;
		//
		func(forward_id,scale,{levelsets[0],dx,rho,area,false,cross_interface,false});
		func(backward_id,-scale,{levelsets[1],dx,rho,area,false,cross_interface,false});
		//
	} else if( face_id.depth > 0 ) {
		//
		const auto &high_layer = *layers[face_id.depth-1];
		const double T_scale = 1.0 / (dx * 0.75);
		//
		if( layer.fill_flags(p_forward) && layer.active_cells.active(p_backward)) {
			//
			std::vector<cell_id3> indices;
			vec3i ivec, jvec;
			if( face_id.dim == 0 ) {
				ivec = vec3i(0,1,0);
				jvec = vec3i(0,0,1);
			} else if( face_id.dim == 1 ) {
				ivec = vec3i(1,0,0);
				jvec = vec3i(0,0,1);
			} else if( face_id.dim == 2 ) {
				ivec = vec3i(1,0,0);
				jvec = vec3i(0,1,0);
			}
			//
			vec3i forward_bottom = 2*(p_backward+vec3i(face_id.dim==0,face_id.dim==1,face_id.dim==2));
			uint_type backward_index = layer.active_cells(p_backward);
			indices.push_back({face_id.depth,p_backward,backward_index});
			for( int ii=0; ii<2; ++ii ) for( int jj=0; jj<2; ++jj ) {
				vec3i hpi = forward_bottom+ii*ivec+jj*jvec;
				indices.push_back({(char)(face_id.depth-1),hpi,high_layer.active_cells(hpi)});
			}
			double min_levelset (0.0), max_levelset (0.0);
			for( unsigned n=0; n<indices.size(); ++n ) {
				min_levelset = std::min(min_levelset,(double)levelset[indices[n].index]);
				max_levelset = std::max(max_levelset,(double)levelset[indices[n].index]);
			}
			const bool cross_interface = min_levelset * max_levelset < 0.0;
			const Real rho = utility::fraction(min_levelset,max_levelset);
			//
			func(indices[0],-T_scale,{levelset[indices[0].index],dx,rho,area,true,cross_interface,false});
			for( unsigned n=1; n<indices.size(); ++n ) {
				func(indices[n],0.25*T_scale,{levelset[indices[n].index],dx,rho,area,true,cross_interface,false});
			}
			//
		} else if( layer.active_cells.active(p_forward) && layer.fill_flags(p_backward)) {
			//
			std::vector<cell_id3> indices;
			vec3i ivec, jvec;
			if( face_id.dim == 0 ) {
				ivec = vec3i(0,1,0);
				jvec = vec3i(0,0,1);
			} else if( face_id.dim == 1 ) {
				ivec = vec3i(1,0,0);
				jvec = vec3i(0,0,1);
			} else if( face_id.dim == 2 ) {
				ivec = vec3i(1,0,0);
				jvec = vec3i(0,1,0);
			}
			vec3i backward_bottom = 2*p_forward-vec3i(face_id.dim==0,face_id.dim==1,face_id.dim==2);
			uint_type forward_index = layer.active_cells(p_forward);
			indices.push_back({face_id.depth,p_forward,forward_index});
			for( int ii=0; ii<2; ++ii ) for( int jj=0; jj<2; ++jj ) {
				vec3i hpi = backward_bottom+ii*ivec+jj*jvec;
				indices.push_back({(char)(face_id.depth-1),hpi,high_layer.active_cells(hpi)});
			}
			//
			double min_levelset (0.0), max_levelset (0.0);
			for( unsigned n=0; n<indices.size(); ++n ) {
				min_levelset = std::min(min_levelset,(double)levelset[indices[n].index]);
				max_levelset = std::max(max_levelset,(double)levelset[indices[n].index]);
			}
			const bool cross_interface = min_levelset * max_levelset < 0.0;
			const Real rho = utility::fraction(min_levelset,max_levelset);
			//
			func(indices[0],T_scale,{levelset[indices[0].index],dx,rho,area,true,cross_interface,false});
			for( unsigned n=1; n<indices.size(); ++n ) {
				func(indices[n],-0.25*T_scale,{levelset[indices[n].index],dx,rho,area,true,cross_interface,false});
			}
		}
	}
}
//
size_t grid3::get_compromised_gradient_count() const {
	//
	std::vector<size_t> count_bucket(parallel.get_thread_num());
	iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
		get_scaled_gradient(face_id,[&]( const cell_id3 &cell_id, double value, const gradient_info3 &info ) {
			if( info.compromised ) {
				count_bucket[tid] ++;
			}
		});
	});
	return std::accumulate(count_bucket.begin(),count_bucket.end(),0);
}
//
size_t grid3::get_surface_T_junction_count() const {
	//
	std::vector<size_t> count_bucket(parallel.get_thread_num());
	iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
		get_scaled_gradient(face_id,[&]( const cell_id3 &cell_id, double value, const gradient_info3 &info ) {
			if( info.t_junction && info.cross_interface ) {
				count_bucket[tid] ++;
			}
		});
	});
	return std::accumulate(count_bucket.begin(),count_bucket.end(),0);
}
//
void grid3::get_scaled_gradient( const face_id3 &face_id, std::function<void( const cell_id3 &cell_id, double value, const gradient_info3 &info )> func ) const {
	//
	const auto &layer = *layers[face_id.depth];
	const Real area = this->area[face_id.index];
	if( area ) {
		//
		using gradient_info_DIM = gradient_info3;
		using cell_id_DIM = cell_id3;
		//
		double fluid_average (0.0);
		double W1 (0.0), W2 (0.0);
		bool has_solid (false);
		//
		auto compute_fluid_average = [&]() {
			get_gradient(face_id,[&]( const cell_id_DIM &cell_id, double value, const gradient_info_DIM &info ) {
				if( info.levelset < 0.0 ) {
					W1 += info.levelset;
					W2 += value * info.levelset;
				}
				fluid_average += value * info.levelset;
				has_solid = has_solid || solid_cell[cell_id.index];
			});
			if( std::abs(W1) < param.eps ) W1 = std::copysign(param.eps,W1);
			if( std::abs(W2) < param.eps ) W2 = std::copysign(param.eps,W2);
		};
		//
		bool rescale_computed (false), t_junction_compute (false);
		double scale (1.0);
		//
		get_gradient(face_id,[&]( const cell_id_DIM &cell_id, double value, const gradient_info_DIM &info ) {
			//
			if( info.cross_interface && ! rescale_computed ) {
				compute_fluid_average();
				if( info.t_junction ) {
					if( ! param.first_order ) t_junction_compute = true;
				} else {
					scale = 1.0 / std::max(param.eps,(double)info.rho);
				}
				rescale_computed = true;
			}
			//
			if( info.levelset < 0.0 ) {
				if( t_junction_compute ) {
					gradient_info_DIM new_info (info);
					double a;
					if( fluid_average * W2 > 0.0 ) {
						a = fluid_average / W2 * value;
					} else {
						new_info.compromised = true;
						if( param.clamp_order ) {
							if( has_solid ) a = param.clamp_solid_eps * value;
							else a = param.clamp_fluid_eps * value;
						} else {
							a = fluid_average / W1;
							if( has_solid ) {
								if( std::abs(a) < param.clamp_solid_eps ) a = std::copysign(param.clamp_solid_eps,a);
							} else {
								if( std::abs(a) < param.clamp_fluid_eps ) a = std::copysign(param.clamp_fluid_eps,a);
							}
						}
					}
					func(cell_id,a,new_info);
				} else {
					func(cell_id,scale*value,info);
				}
			}
		});
	}
}
//
void grid3::get_divergence( const cell_id3 &cell_id, std::function<void( const face_id3 &face_id, double value0, double value1 )> func ) const {
	//
	for( char dim : DIMS3 ) {
		iterate_face_neighbors(cell_id,dim,[&]( const face_id3 &face_id ){
			//
			const auto &face_layer = *layers[face_id.depth];
			const Real area = this->area[face_id.index];
			uint_type column = face_layer.active_faces[dim](face_id.pi);
			get_gradient(face_id,[&]( const cell_id3 &neighbor_cell_id, double value, const gradient_info3 &info ) {
				if( info.levelset < 0.0 && cell_id.index == neighbor_cell_id.index ) {
					const double t = (info.t_junction ? 0.75 : 1.0) * info.dx * info.dx * info.dx;
					const double scale0 = t * area;
					const double scale1 = t * (1.0-area);
					func(face_id,scale0*value,scale1*value);
				}
			});
		});
	}
}
//
void grid3::get_unmofidied_divergence( const cell_id3 &cell_id, std::function<void( const face_id3 &face_id, double value )> func ) const {
	//
	for( char dim : DIMS3 ) {
		iterate_face_neighbors(cell_id,dim,[&]( const face_id3 &face_id ){
			this->get_gradient(face_id,[&]( const cell_id3 &neighbor_cell_id, double value, const gradient_info3 &info ) {
				if( cell_id.index == neighbor_cell_id.index ) {
					const double scale = (info.t_junction ? 0.75 : 1.0) * info.dx * info.dx * info.dx;
					func(face_id,scale*value);
				}
			});
		});
	}
}
//
void grid3::add_surfacetense_force( double coeff, double dt ) {
	//
	// Compute the gradient of levelset
	std::vector<Real> grad_levelset(face_count);
	iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
		bool valid (false);
		double sum (0.0);
		this->get_gradient(face_id,[&]( const cell_id3 &cell_id, double value, const grid3::gradient_info3 &info ) {
			sum += levelset[cell_id.index] * value;
		});
		grad_levelset[face_id.index] = sum;
	});
	//
	if( param.normalize_laplacian ) {
		//
		// Compute the full gradient maginitude per cell
		std::vector<Real> cell_grad_levelset_magnitude(cell_count);
		iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
			double dx = get_cell_dx(cell_id);
			vec3d p = get_cell_position(cell_id);
			cell_grad_levelset_magnitude[cell_id.index] = 0.5*vec3d(
				sample_face(p+0.5*vec3d(dx,0.0,0.0),0,grad_levelset)+sample_face(p-0.5*vec3d(dx,0.0,0.0),0,grad_levelset),
				sample_face(p+0.5*vec3d(0.0,dx,0.0),1,grad_levelset)+sample_face(p-0.5*vec3d(0.0,dx,0.0),1,grad_levelset),
				sample_face(p+0.5*vec3d(0.0,0.0,dx),2,grad_levelset)+sample_face(p-0.5*vec3d(0.0,0.0,dx),2,grad_levelset)
			).len();
		});
		//
		// Scale gradient by this
		iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
			vec3d p = get_face_position(face_id);
			grad_levelset[face_id.index] /= sample_cell(p,cell_grad_levelset_magnitude);
		});
	}
	//
	// Compute the divergence of the levelset gradient
	std::vector<Real> lap_levelset(cell_count);
	iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
		double dx = get_cell_dx(cell_id);
		vec3d p = get_cell_position(cell_id);
		double value = 
			(sample_face(p+0.5*vec3d(dx,0.0,0.0),0,grad_levelset)-sample_face(p-0.5*vec3d(dx,0.0,0.0),0,grad_levelset)+
			 sample_face(p+0.5*vec3d(0.0,dx,0.0),1,grad_levelset)-sample_face(p-0.5*vec3d(0.0,dx,0.0),1,grad_levelset)+
			 sample_face(p+0.5*vec3d(0.0,0.0,dx),2,grad_levelset)-sample_face(p-0.5*vec3d(0.0,0.0,dx),2,grad_levelset)) / dx;
		lap_levelset[cell_id.index] = value;
	});
	//
	// Add one step smoothing
	if( param.surftens_smooth_count ) {
		unsigned count = param.surftens_smooth_count;
		while( count -- ) {
			std::vector<Real> lap_levelset_save(lap_levelset);
			iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
				unsigned wsum (0);
				unsigned w (2);
				Real new_value (w*lap_levelset_save[cell_id.index]);
				wsum += w;
				iterate_cell_neighbors(cell_id,[&]( char dim, const cell_id3 &neighbor_cell_id ) {
					new_value += lap_levelset_save[neighbor_cell_id.index];
					wsum += 1;
				});
				lap_levelset[cell_id.index] = new_value / (Real)wsum;
			});
		}
	}
	//
	// Add surface tension force
	iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
		//
		Real avg_levelest_inside (0.0);
		Real avg_laplacian_inside (0.0);
		unsigned count_inside (0);
		//
		// Compute the average position and the levelset of inside cells
		this->get_scaled_gradient(face_id,[&]( const cell_id3 &cell_id, double value, const grid3::gradient_info3 &info ) {
			if( info.cross_interface ) {
				avg_levelest_inside += info.levelset;
				avg_laplacian_inside += lap_levelset[cell_id.index];
				count_inside ++;
			}
		});
		//
		if( count_inside ) {
			avg_levelest_inside /= count_inside;
			avg_laplacian_inside /= count_inside;
			assert( avg_levelest_inside <= 0.0 );
			//
			// Compute the average position and the levelset of outside cells
			Real avg_levelest_outside (0.0);
			Real avg_laplacian_outside (0.0);
			unsigned count_outside (0);
			//
			this->get_gradient(face_id,[&]( const cell_id3 &cell_id, double value, const grid3::gradient_info3 &info ) {
				if( info.levelset > 0.0 ) {
					avg_levelest_outside += info.levelset;
					avg_laplacian_outside += lap_levelset[cell_id.index];
					count_outside ++;
				}
			});
			assert( count_outside );
			//
			avg_levelest_outside /= count_outside;
			avg_laplacian_outside /= count_outside;
			assert( avg_levelest_outside >= 0.0 );
			//
			// Compute the curvature
			const double curvature = (avg_levelest_outside * avg_laplacian_inside - avg_levelest_inside * avg_laplacian_outside) / (avg_levelest_outside - avg_levelest_inside);
			//
			double sum (0.0);
			this->get_scaled_gradient(face_id,[&]( const cell_id3 &cell_id, double value, const grid3::gradient_info3 &info ) {
				sum += value * curvature;
			});
			//
			// Add surface tension force
			velocity[face_id.index] += dt * coeff * sum;
		}
	});
}
//
void grid3::extrapolate( std::function<double(const vec3d &p)> solid_func ) {
	//
	scoped_timer timer(this);
	assert(layers.size());
	shape3 shape = layers[0]->shape;
	double dx = layers[0]->dx;
	const double max_dist = dx * shape.max();
	//
	std::vector<cell_id3> cell_ids(cell_count);
	iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
		cell_ids[cell_id.index] = cell_id;
	});
	//
	std::vector<Real> levelset_new(levelset);
	auto velocity_extrapolate = [&]() {
		//
		std::vector<char> face_fixed(face_count);
		iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
			bool fix (false);
			get_scaled_gradient(face_id,[&]( const cell_id3 &cell_id, double value, const gradient_info3 &info ) {
				if( info.area && info.rho ) fix = true;
			});
			face_fixed[face_id.index] = fix;
			if( ! fix ) velocity[face_id.index] = 0.0;
		});
		//
		std::vector<char> cell_fixed(face_count);
		std::vector<vec3d> cell_velocity(cell_count);
		//
		iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
			vec3d wsum;
			for( char dim : DIMS3 ) {
				double value (0.0);
				iterate_face_neighbors(cell_id,dim,[&]( const face_id3 &face_id ) {
					if( face_fixed[face_id.index] ) {
						value += velocity[face_id.index];
						wsum[dim] += 1.0;
					}
				});
				if( wsum[dim] ) {
					cell_velocity[cell_id.index][dim] = value / wsum[dim];
				}
			}
			if( wsum[0] && wsum[1] && wsum[2] ) {
				cell_fixed[cell_id.index] = true;
			} else {
				cell_fixed[cell_id.index] = false;
				cell_velocity[cell_id.index] = vec3d();
			}
		});
		//
		std::vector<Real> levelset_combined(levelset.size());
		iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
			levelset_combined[cell_id.index] = std::max(levelset[cell_id.index],-(Real)solid_func(get_cell_position(cell_id)));
		});
		unstructured_extrapolator3<vec3d>::extrapolate(
			[&]( size_t index ) {
				return get_cell_position(cell_ids[index]);
			},
			[&]( size_t index, std::function<void( size_t j )> func ) {
				iterate_cell_neighbors(cell_ids[index],[&]( char dim, const cell_id3& cell_id){ func(cell_id.index); });
			},levelset_combined,cell_fixed,cell_velocity,max_dist);
		//
		iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
			if( ! face_fixed[face_id.index] ) {
				vec3d vel;
				double wsum (0.0);
				get_gradient(face_id,[&]( const cell_id3 &cell_id, double value, const gradient_info3 &info ) {
					if( cell_fixed[cell_id.index] ) {
						double w = std::abs(value);
						vel += w*cell_velocity[cell_id.index];
						wsum += w;
					}
				});
				if( wsum ) {
					velocity[face_id.index] = vel[face_id.dim] / wsum;
				}
			}
		});
	};
	//
	auto levelset_extrapolate = [&]() {
		//
		std::vector<char> fixed(cell_count);
		iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
			const uint_type n = cell_id.index;
			const double levelset0 = levelset[n];
			iterate_cell_neighbors(cell_ids[n],[&]( char dim, const cell_id3& neighbor_cell_id ){
				const double levelset1 = levelset[neighbor_cell_id.index];
				if( levelset0 * levelset1 < 0.0 ) {
					fixed[n] = true;
				}
			});
		});
		const double sqrt3 = sqrt(3.0);
		iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
			const uint_type n = cell_id.index;
			const double levelset0 = levelset[n];
			const double dx = get_cell_dx(cell_id);
			if( fixed[n] ) {
				if( param.simple_redistance ) {
					if( std::abs(levelset[n]) > sqrt3*dx ) {
						levelset_new[n] = std::copysign(sqrt3*dx,levelset[n]);
					}
				} else {
					const double glen = get_upwind_gradient(cell_id,levelset).len();
					if( glen ) levelset_new[cell_id.index] /= glen;
				}
			} else {
				levelset_new[n] = std::copysign(sqrt3*dx,levelset[n]);
			}
		});
		//
		if( param.pde_update_count ) {
			for( unsigned n=0; n<param.pde_update_count; ++n ) {
				std::vector<Real> levelset_prev(levelset_new);
				iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
					if( ! fixed[cell_id.index] ) {
						const double glen = get_upwind_gradient(cell_id,levelset_prev).len();
						const double dx = get_cell_dx(cell_id);
						const double v = levelset_prev[cell_id.index];
						const double sgn = v / sqrt(v*v+dx*dx);
						levelset_new[cell_id.index] += 0.5*sgn*dx*(1.0-glen);
					}
				});
			}
		} else {
			unstructured_fastmarch3::fastmarch(
				[&]( size_t index ) {
					return get_cell_position(cell_ids[index]);
				},
				[&]( size_t index, std::function<void( size_t j )> func ) {
					iterate_cell_neighbors(cell_ids[index],[&]( char dim, const cell_id3& cell_id ){ func(cell_id.index); });
				},levelset_new,fixed,max_dist,parallel,meshutility.get());
		}
	};
	//
	timer.tick(); console::dump( "Extrapolating velocity and levelset..." );
	//
	std::vector<std::function<void()> > operations;
	operations.push_back(levelset_extrapolate);
	operations.push_back(velocity_extrapolate);
	parallel.run(operations);
	levelset = levelset_new;
	//
	console::dump( "Done. Took %s.\n", timer.stock("extrapolate").c_str());
}
//
void grid3::extrapolate_toward_solid( std::function<double(const vec3d &p)> solid_func ) {
	//
	scoped_timer timer(this);
	timer.tick(); console::dump( "Extrapolating levelset towards solid..." );
	//
	auto solid_derivative = [&]( const vec3d &p, double dx ) {
		vec3d derivative;
		for( char dim : DIMS3 ) {
			derivative[dim] = (solid_func(p+0.5*dx*vec3d(dim==0,dim==1,dim==2))-solid_func(p-0.5*dx*vec3d(dim==0,dim==1,dim==2))) / dx;
		}
		return derivative;
	};
	//
	std::vector<Real> new_fluid (levelset);
	parallel.for_each( layers.size(), [&]( size_t depth ) {
		//
		auto &layer = *layers[depth];
		const shape3 &shape = layer.shape;
		double dx = layer.dx;
		//
		layer.active_cells.const_parallel_actives([&]( int i, int j, int k, const auto &it ) {
			vec3d p = dx*vec3i(i,j,k).cell();
			double solid = solid_func(p);
			if( solid < 0.0 ) {
				vec3d derivative = solid_derivative(p,dx);
				if( derivative[1] > 0.0 ) {
					derivative[1] = 0.0;
					derivative.normalize();
				}
				new_fluid[it()] = this->sample_levelset(p - solid * derivative);
			}
		});
	});
	//
	levelset = new_fluid;
	console::dump( "Done. Took %s.\n", timer.stock("extrapolate_toward_solid").c_str());
}
//
static double grid_kernel( const vec3d &r, double dx, double dy, double dz ) {
	//
	const double x = std::abs(r[0]) / dx;
	const double y = std::abs(r[1]) / dy;
	const double z = std::abs(r[2]) / dz;
	return std::max(0.0,1.0-x) * std::max(0.0,1.0-y) * std::max(0.0,1.0-z);
}
//
static double grid_kernel( const vec3d &r, double dx ) {
	//
	const double x = std::abs(r[0]) / dx;
	const double y = std::abs(r[1]) / dx;
	const double z = std::abs(r[2]) / dx;
	return std::max(0.0,1.0-x) * std::max(0.0,1.0-y) * std::max(0.0,1.0-z);
}
//
static double MLS_interpolate( const vec3d &p, const std::vector<grid3::point_info3> &points, unsigned char type, bool use_const_diag=false ) {
	//
	double max_dx (0.0);
	unsigned size = points.size();
	for( unsigned n=0; n<points.size(); ++n ) {
		max_dx = std::max(max_dx,points[n].dx);
	}
	const double scale = 1.0/max_dx;
	Eigen::MatrixXd A(size,4);
	for( unsigned row=0; row<size; ++row ) {
		const vec3d &q = points[row].position;
		A(row,0) = scale * q[0];
		A(row,1) = scale * q[1];
		A(row,2) = scale * q[2];
		A(row,3) = 1.0;
	}
	Eigen::MatrixXd x(1,4);
	x << scale * p[0], scale * p[1], scale * p[2], 1.0;
	Eigen::MatrixXd a(size,1);
	for( unsigned n=0; n<size; ++n ) a(n,0) = points[n].value;
	//
	Eigen::VectorXd W_vec(size);
	for( unsigned row=0; row<size; ++row ) {
		W_vec(row) = use_const_diag ? 1.0 : points[row].weight;
	}
	auto W = W_vec.asDiagonal();
	Eigen::MatrixXd AtW = A.transpose() * W;
	Eigen::MatrixXd AtWA = AtW * A;
	const double result = (x * AtWA.llt().solve(AtW * a))(0,0);
	if( utility::is_nan(result)) {
		console::dump("type = %s, size = %u\n",type==0 ? "cell" : "face",size);
		for( unsigned row=0; row<size; ++row ) {
			const vec3d &q = points[row].position;
			console::dump("max_dx=%.2e, p=(%.2e,%.2e,%.2e),P=(%.2e,%.2e,%.2e), W=%.2e, V=%.2e\n",
				max_dx,p[0],p[1],p[2],q[0],q[1],q[2],W_vec(row),a(row,0));
		}
		exit(0);
		return 0.0;
	} else {
		return result;
	}
}
//
static unsigned mirror_samples( const grid3 &grid, const vec3d &p, std::vector<grid3::point_info3> &points, double eps ) {
	//
	thread_local std::vector<grid3::point_info3> mirror_points;
	mirror_points.clear();
	//
	const double &dx = grid.layers[0]->dx;
	const shape3 &shape = grid.layers[0]->shape;
	vec3d domain_size = dx*vec3d(shape[0],shape[1],shape[2]);
	//
	for( unsigned n=0; n<points.size(); ++n ) {
		const double cellsize = points[n].dx;
		const double h = 0.1*cellsize;
		for( int dim : DIMS3 ) {
			if( points[n].position[dim] > h && points[n].position[dim]-cellsize < h ) {
				vec3d mirror_pos = points[n].position;
				mirror_pos[dim] = -mirror_pos[dim];
				double w = grid_kernel(mirror_pos-p,cellsize);
				if( eps ) w = std::max(w,eps);
				if( w ) {
					grid3::point_info3 info;
					info.dx = cellsize;
					info.value = points[n].value;
					info.weight = w;
					info.position = mirror_pos;
					points.push_back(info);
				}
			} else if( points[n].position[dim] < domain_size[dim]-h &&
					   points[n].position[dim]+cellsize > domain_size[dim]-h ) {
				vec3d mirror_pos = points[n].position;
				mirror_pos[dim] = 2.0*domain_size[dim]-mirror_pos[dim];
				double w = grid_kernel(mirror_pos-p,cellsize);
				if( eps ) w = std::max(w,eps);
				if( w ) {
					grid3::point_info3 info;
					info.dx = cellsize;
					info.value = points[n].value;
					info.weight = w;
					info.position = mirror_pos;
					points.push_back(info);
				}
			}
		}
	}
	//
	if( mirror_points.size()) {
		points.insert(points.end(),mirror_points.begin(),mirror_points.end());
	}
	return mirror_points.size();
}
//
bool grid3::find_cell( const vec3d &p, cell_id3 &cell_id ) const {
	//
	for( char depth=0; depth < layers.size(); ++depth ) {
		//
		const auto &layer = *layers[depth];
		const vec3i pi = layer.shape.find_cell(p/layer.dx);
		if( layer.active_cells.active(pi)) {
			cell_id = {depth,pi,layer.active_cells(pi)};
			return true;
		}
	}
	return false;
}
//
bool grid3::find_face( const vec3d &p, char dim, face_id3 &face_id ) const {
	//
	for( char depth=0; depth < layers.size(); ++depth ) {
		//
		const auto &layer = *layers[depth];
		vec3i pi = layer.shape.find_face(p/layer.dx,dim);
		if( pi[dim] == 0 ) pi[dim] ++;
		if( pi[dim] == layer.shape[dim] ) pi[dim] --;
		//
		if( layer.active_faces[dim].active(pi)) {
			face_id = {depth,dim,pi,layer.active_faces[dim](pi)};
			return true;
		}
	}
	return false;
}
//
double grid3::sample_cell( const vec3d &p_in, const std::vector<Real> &cell_values, Real *min_max_values ) const {
	//
	double value (0.0);
	if( try_interp_cell(p_in,cell_values,value,min_max_values)) {
		return value;
	}
	//
	const shape3 shape = layers[0]->shape;
	const double dx = layers[0]->dx;
	vec3d p;
	for( char dim : DIMS3 ) {
		p[dim] = std::min(dx*(shape[dim]),std::max(0.0,p_in[dim]));
	}
	//
	if( min_max_values ) {
		min_max_values[0] = std::numeric_limits<Real>::max();
		min_max_values[1] = std::numeric_limits<Real>::min();
	}
	//
	thread_local std::vector<point_info3> points; points.clear();
	gather_neighbor_cells(p,points,cell_values);
	//
	if( min_max_values ) {
		for( unsigned n=0; n<points.size(); ++n ) {
			min_max_values[0] = std::min(min_max_values[0],points[n].value);
			min_max_values[1] = std::max(min_max_values[1],points[n].value);
		}
	}
	//
	return MLS_interpolate(p,points,0,param.MLS_constant_diag);
}
//
double grid3::sample_levelset( const vec3d &p, Real *min_max_values ) const {
	//
	const shape3 shape = layers[0]->shape;
	const double dx = layers[0]->dx;
	auto clamp_ceil = [&]( double x ) {
		return std::max(x,p[1]-dx*(shape[1]-0.5));
	};
	//
	Real local_min_max_values[2] = { 0.0, 0.0 };
	double result = sample_cell(p,levelset,local_min_max_values);
	result = clamp_ceil(std::max(std::min(result,(double)local_min_max_values[1]),(double)local_min_max_values[0]));
	//local_min_max_values
	if( min_max_values ) {
		min_max_values[0] = local_min_max_values[0];
		min_max_values[1] = local_min_max_values[1];
	}
	//
	return result;
}
//
double grid3::sample_face( const vec3d &p_in, char dim, const std::vector<Real> &face_values, Real *min_max_values ) const {
	//
	double value (0.0);
	if( try_interp_face(p_in,dim,face_values,value,min_max_values)) {
		return value;
	}
	//
	const shape3 shape = layers[0]->shape;
	const double dx = layers[0]->dx;
	vec3d p;
	for( char dim : DIMS3 ) {
		p[dim] = std::min(dx*(shape[dim]),std::max(0.0,p_in[dim]));
	}
	//
	if( min_max_values ) {
		min_max_values[0] = std::numeric_limits<Real>::max();
		min_max_values[1] = std::numeric_limits<Real>::min();
	}
	//
	thread_local std::vector<point_info3> points; points.clear();
	gather_neighbor_faces(p,dim,points,face_values);
	//
	if( min_max_values ) {
		for( unsigned n=0; n<points.size(); ++n ) {
			min_max_values[0] = std::min(min_max_values[0],points[n].value);
			min_max_values[1] = std::max(min_max_values[1],points[n].value);
		}
	}
	//
	return MLS_interpolate(p,points,1,param.MLS_constant_diag);
}
//
double grid3::sample_velocity( const vec3d &p, char dim, Real *min_max_values ) const {
	return sample_face(p,dim,velocity,min_max_values);
}
//
vec3d grid3::sample_velocity( const vec3d &p ) const {
	return vec3d(sample_velocity(p,0),sample_velocity(p,1),sample_velocity(p,2));
}
//
double grid3::sample_solid_face_velocity( const face_id3 &face_id, std::function<vec3d( const vec3d &p )> func ) const {
	//
	const int dim = face_id.dim;
	const double dx = get_face_dx(face_id);
	const vec3d p = get_face_position(face_id);
	vec3i ivec, jvec;
	if( dim == 0 ) {
		ivec = vec3i(0,1,0);
		jvec = vec3i(0,0,1);
	} else if( dim == 1 ) {
		ivec = vec3i(1,0,0);
		jvec = vec3i(0,0,1);
	} else if( dim == 2 ) {
		ivec = vec3i(1,0,0);
		jvec = vec3i(0,1,0);
	}
	//
	double wsum (0.0);
	double usum (0.0);
	const vec3d o = p-0.5*dx*(ivec+jvec);
	//
	for( int dir0=0; dir0<=1; ++dir0 ) for( int dir1=0; dir1<=1; ++dir1 ) {
		const double u = func(o+dx*(dir0*ivec+dir1*jvec))[dim];
		if( u ) {
			usum += u;
			wsum += 1.0;
		}
	}
	//
	if( wsum ) {
		return usum / wsum;
	} else {
		return 0.0;
	}
}
//
bool grid3::try_interp_cell( const vec3d &p_in, const std::vector<Real> &cell_values, double &value, Real *min_max_values ) const {
	//
	if( min_max_values ) {
		min_max_values[0] = std::numeric_limits<Real>::max();
		min_max_values[1] = std::numeric_limits<Real>::min();
	}
	//
	assert(layers.size());
	shape3 shape = layers[0]->shape;
	double dx = layers[0]->dx;
	//
	vec3d p;
	for( char dim2 : DIMS3 ) {
		p[dim2] = std::min(dx*(shape[dim2]),std::max(0.0,p_in[dim2]));
	}
	//
	bool prev_hit (false);
	bool invalid (false);
	int hit_count (0);
	double sum (0.0);
	//
	for( char depth=0; depth < layers.size(); ++depth ) {
		//
		const auto &layer = *layers[depth];
		const double dx = layer.dx;
		//
		vec3i indices[8];
		double coef[8];
		array_interpolator3::interpolate_coef(layer.shape,p/layer.dx-vec3d(0.5,0.5,0.5),indices,coef);
		//
		bool exit_after_loop = prev_hit;
		for( int n=0; n<8; ++n ) {
			if( layer.active_cells.active(indices[n])) {
				prev_hit = true;
				hit_count ++;
				if( exit_after_loop ) {
					invalid = true;
					break;
				}
				const uint_type index = layer.active_cells(indices[n]);
				sum += coef[n] * cell_values[index];
				//
				if( min_max_values ) {
					min_max_values[0] = std::min(min_max_values[0],cell_values[index]);
					min_max_values[1] = std::max(min_max_values[1],cell_values[index]);
				}
			}
		}
		if( exit_after_loop ) break;
	}
	//
	if( ! invalid && hit_count == 8 ) {
		value = sum;
		return true;
	} else {
		value = 0.0;
		return false;
	}
	//
	return false;
}
//
bool grid3::try_interp_face( const vec3d &p_in, char dim, const std::vector<Real> &face_values, double &value, Real *min_max_values ) const {
	//
	if( min_max_values ) {
		min_max_values[0] = std::numeric_limits<Real>::max();
		min_max_values[1] = std::numeric_limits<Real>::min();
	}
	//
	assert(layers.size());
	shape3 shape = layers[0]->shape;
	double dx = layers[0]->dx;
	//
	vec3d p;
	for( char dim2 : DIMS3 ) {
		p[dim2] = std::min(dx*(shape[dim2]),std::max(0.0,p_in[dim2]));
	}
	//
	bool prev_hit (false);
	bool invalid (false);
	int hit_count (0);
	double wsum (0.0);
	double sum (0.0);
	//
	auto fix_interpolation = [this]( int dim, const shape3 &shape, vec3i indices[8], double coef[8] ) {
		for( int n=0; n<8; ++n ) {
			if( indices[n][dim] == 0 ) for( int m=0; m<8; ++m ) {
				if( indices[n]+vec3i(dim==0,dim==1,dim==2) == indices[m] ) {
					indices[n] = indices[m];
					break;
				}
			}
			if( indices[n][dim] == shape[dim]-1 ) for( int m=0; m<8; ++m ) {
				if( indices[n]-vec3i(dim==0,dim==1,dim==2) == indices[m] ) {
					indices[n] = indices[m];
					break;
				}
			}
		}
	};
	//
	for( char depth=0; depth < layers.size(); ++depth ) {
		//
		const auto &layer = *layers[depth];
		const double dx = layer.dx;
		//
		vec3i indices[8];
		double coef[8];
		vec3d pos = vec3d(p[0]/dx-0.5*(dim!=0),p[1]/dx-0.5*(dim!=1),p[2]/dx-0.5*(dim!=2));
		array_interpolator3::interpolate_coef(layer.shape.face(dim),pos,indices,coef);
		fix_interpolation(dim,layer.shape.face(dim),indices,coef);
		//
		bool exit_after_loop = prev_hit;
		for( int n=0; n<8; ++n ) {
			if( layer.active_faces[dim].active(indices[n])) {
				//
				prev_hit = true;
				hit_count ++;
				if( exit_after_loop ) {
					invalid = true;
					break;
				}
				const uint_type index = layer.active_faces[dim](indices[n]);
				sum += coef[n] * face_values[index];
				//
				if( min_max_values ) {
					min_max_values[0] = std::min(min_max_values[0],face_values[index]);
					min_max_values[1] = std::max(min_max_values[1],face_values[index]);
				}
			}
		}
		if( exit_after_loop ) break;
	}
	//
	if( ! invalid && hit_count == 8 ) {
		value = sum;
		return true;
	} else {
		value = 0.0;
		return false;
	}
}
//
void grid3::gather_neighbor_faces( const vec3d &p_in, char dim, std::vector<point_info3> &points, const std::vector<Real> &values ) const {
	//
	assert(layers.size());
	shape3 shape = layers[0]->shape;
	double dx = layers[0]->dx;
	//
	vec3d p;
	for( char dim2 : DIMS3 ) {
		p[dim2] = std::min(dx*(shape[dim2]),std::max(0.0,p_in[dim2]));
	}
	auto fix_interpolation = []( int dim, const shape3 &shape, vec3i indices[8] ) {
		for( int n=0; n<8; ++n ) {
			if( indices[n][dim] == 0 ) for( int m=0; m<8; ++m ) {
				if( indices[n]+vec3i(dim==0,dim==1,dim==2) == indices[m] ) {
					indices[n] = indices[m];
					break;
				}
			}
			if( indices[n][dim] == shape[dim]-1 ) for( int m=0; m<8; ++m ) {
				if( indices[n]-vec3i(dim==0,dim==1,dim==2) == indices[m] ) {
					indices[n] = indices[m];
					break;
				}
			}
		}
	};
	//
	bool prev_hit (false);
	char hit_depth (0);
	for( char depth=0; depth < layers.size(); ++depth ) {
		//
		const auto &layer = *layers[depth];
		const double dx = layer.dx;
		//
		vec3i indices[8];
		double coef[8];
		vec3d pos = vec3d(p[0]/dx-0.5*(dim!=0),p[1]/dx-0.5*(dim!=1),p[2]/dx-0.5*(dim!=2));
		array_interpolator3::interpolate_coef(layer.shape.face(dim),pos,indices,coef);
		fix_interpolation(dim,layer.shape.face(dim),indices);
		//
		bool exit_after_loop = prev_hit;
		for( int n=0; n<8; ++n ) if( coef[n] ) {
			if( layer.active_faces[dim].active(indices[n])) {
				prev_hit = true;
				hit_depth = depth;
				const uint_type index = layer.active_faces[dim](indices[n]);
				point_info3 info;
				info.dx = layer.dx;
				info.value = values[index];
				info.weight = coef[n];
				info.position = layer.dx*indices[n].face(dim);
				info.index = index;
				points.push_back(info);
			} else if( param.accurate_interpolation && layer.ghost_face_indices[dim].active(indices[n]) ) {
				//
				for( const auto &ghost_face : ghost_faces[layer.ghost_face_indices[dim](indices[n])]) {
					//
					point_info3 info;
					double dx(layer.dx), dy(layer.dx), dz(layer.dx);
					const auto &ghost_position = ghost_face.p;
					if( p[0] < ghost_position[0] && (ghost_face.shrink_info & shrink_left )) dx *= 0.5;
					if( p[0] > ghost_position[0] && (ghost_face.shrink_info & shrink_right )) dx *= 0.5;
					if( p[1] < ghost_position[1] && (ghost_face.shrink_info & shrink_bottom )) dy *= 0.5;
					if( p[1] > ghost_position[1] && (ghost_face.shrink_info & shrink_top )) dy *= 0.5;
					if( p[2] < ghost_position[2] && (ghost_face.shrink_info & shrink_front )) dz *= 0.5;
					if( p[2] > ghost_position[2] && (ghost_face.shrink_info & shrink_back )) dz *= 0.5;
					//
					const double w = grid_kernel(p-ghost_position,dx,dy,dz);
					if( w ) {
						info.shrink_info = ghost_face.shrink_info;
						info.weight = w;
						info.dx = layer.dx;
						double value (0.0);
						for( const auto &e : ghost_face.combined ) {
							value += e.second*values[e.first];
						}
						info.value = value;
						info.position = ghost_position;
						info.index = std::numeric_limits<uint_type>::max();
						points.push_back(info);
					}
				}
			}
		}
		if( exit_after_loop ) break;
	}
	//
	bool hit_face (false);
	for( char depth=std::max(0,hit_depth-1); depth < layers.size(); ++depth ) {
		//
		const auto &layer = *layers[depth];
		vec3i pi = layer.shape.find_face(p/layer.dx,dim);
		if( pi[dim] == 0 ) pi[dim] ++;
		if( pi[dim] == layer.shape[dim] ) pi[dim] --;
		//
		if( layer.active_faces[dim].active(pi)) {
			const face_id3 &f = {depth,dim,pi,layer.active_faces[dim](pi)};
			iterate_face_neighbors(f,[&]( const face_id3 &face_id ) {
				if( std::find_if(points.begin(),points.end(),[&]( const auto &it ) {
					return it.index == face_id.index;
				}) == points.end()) {
					//
					const vec3d &face_p = get_face_position(face_id);
					const double dx = get_face_dx(face_id);
					double w = grid_kernel(face_p-p,dx);
					if( ! param.accurate_interpolation ) w = std::max(w,param.MLS_eps);
					if( w ) {
						point_info3 info;
						info.dx = dx;
						info.value = values[face_id.index];
						info.weight = w;
						info.position = face_p;
						info.index = face_id.index;
						points.push_back(info);
					}
				}
			});
			hit_face = true;
			break;
		}
	}
	if( ! hit_face ) {
		console::dump( "Error: gather_neighbor_faces failed.\n" );
		exit(0);
	}
	//
	for( char depth=std::max(0,hit_depth-1); depth <= hit_depth; ++depth ) {
		const auto &layer = *layers[depth];
		const double &dx = layer.dx;
		vec3i indices[8];
		double coef[8];
		vec3d pos = vec3d(p[0]/dx-0.5*(dim!=0),p[1]/dx-0.5*(dim!=1),p[2]/dx-0.5*(dim!=2));
		array_interpolator3::interpolate_coef(layer.shape.face(dim),pos,indices,coef);
		for( unsigned n=0; n<8; n++ ) {
			if( indices[n][dim] == 0 || indices[n][dim] == layer.shape[dim] ) {
				const face_id3 &face_id = {depth,dim,indices[n],std::numeric_limits<uint_type>::max()};
				const vec3d &face_p = get_face_position(face_id);
				double w = grid_kernel(face_p-p,dx);
				if( ! param.accurate_interpolation ) w = std::max(w,param.MLS_eps);
				if( w ) {
					int dir = indices[n][dim] == 0 ? 1 : -1;
					const vec3i neighbor_pi = indices[n]+dir*vec3i(dim==0,dim==1,dim==2);
					if( layer.active_faces[dim].safe_active(neighbor_pi)) {
						point_info3 info;
						info.dx = layer.dx;
						info.value = values[layer.active_faces[dim](neighbor_pi)];
						info.weight = w;
						info.position = face_p;
						info.index = face_id.index;
						points.push_back(info);
					}
				}
			}
		}
	}
	//
	mirror_samples(*this,p,points,param.accurate_interpolation ? 0.0 : param.MLS_eps);
}
//
void grid3::gather_neighbor_cells( const vec3d &p_in, std::vector<point_info3> &points, const std::vector<Real> &values ) const {
	//
	assert(layers.size());
	shape3 shape = layers[0]->shape;
	double dx = layers[0]->dx;
	//
	vec3d p;
	for( char dim : DIMS3 ) {
		p[dim] = std::min(dx*(shape[dim]),std::max(0.0,p_in[dim]));
	}
	//
	bool prev_hit (false);
	char hit_depth (0);
	for( char depth=0; depth < layers.size(); ++depth ) {
		//
		const auto &layer = *layers[depth];
		vec3i indices[8];
		double coef[8];
		array_interpolator3::interpolate_coef(layer.shape,p/layer.dx-vec3d(0.5,0.5,0.5),indices,coef);
		//
		bool exit_after_loop = prev_hit;
		for( int n=0; n<8; ++n ) if( coef[n] ) {
			if( layer.active_cells.active(indices[n])) {
				prev_hit = true;
				hit_depth = depth;
				const uint_type index = layer.active_cells(indices[n]);
				point_info3 info;
				info.dx = layer.dx;
				info.value = values[index];
				info.weight = coef[n];
				info.position = layer.dx*indices[n].cell();
				info.index = index;
				points.push_back(info);
			} if( param.accurate_interpolation && layer.ghost_cell_indices.active(indices[n]) ) {
				//
				for( const auto &ghost_cell : ghost_cells[layer.ghost_cell_indices(indices[n])]) {
					double dx (layer.dx), dy (layer.dx), dz (layer.dx);
					const auto &ghost_position = ghost_cell.p;
					if( p[0] < ghost_position[0] && (ghost_cell.shrink_info & shrink_left )) dx *= 0.5;
					if( p[0] > ghost_position[0] && (ghost_cell.shrink_info & shrink_right )) dx *= 0.5;
					if( p[1] < ghost_position[1] && (ghost_cell.shrink_info & shrink_bottom )) dy *= 0.5;
					if( p[1] > ghost_position[1] && (ghost_cell.shrink_info & shrink_top )) dy *= 0.5;
					if( p[2] < ghost_position[2] && (ghost_cell.shrink_info & shrink_front )) dz *= 0.5;
					if( p[2] > ghost_position[2] && (ghost_cell.shrink_info & shrink_back )) dz *= 0.5;
					//
					const double w = grid_kernel(p-ghost_position,dx,dy,dz);
					if( w ) {
						point_info3 info;
						info.shrink_info = ghost_cell.shrink_info;
						info.weight = w;
						info.dx = layer.dx;
						double value (0.0);
						for( const auto &e : ghost_cell.combined ) {
							value += e.second*values[e.first];
						}
						info.value = value;
						info.position = ghost_position;
						info.index = std::numeric_limits<uint_type>::max();
						points.push_back(info);
					}
				}
			}
		}
		if( exit_after_loop ) break;
	}
	//
	bool hit_cell (false);
	for( char depth=std::max(0,hit_depth-1); depth < layers.size(); ++depth ) {
		const auto &layer = *layers[depth];
		const vec3i pi = layer.shape.find_cell(p/layer.dx);
		if( layer.active_cells.active(pi)) {
			//
			iterate_cell_neighbors({depth,pi,layer.active_cells(pi)},[&]( char dim, const cell_id3 &cell_id ) {
				if( std::find_if(points.begin(),points.end(),[&]( const auto &it ) {
					return it.index == cell_id.index;
				}) == points.end()) {
					//
					const vec3d &cell_p = get_cell_position(cell_id);
					const double dx = layers[cell_id.depth]->dx;
					double w = grid_kernel(cell_p-p,dx);
					if( ! param.accurate_interpolation ) w = std::max(w,param.MLS_eps);
					if( w ) {
						point_info3 info;
						info.dx = dx;
						info.value = values[cell_id.index];
						info.weight = w;
						info.position = cell_p;
						info.index = cell_id.index;
						points.push_back(info);
					}
				}
			});
			hit_cell = true;
			break;
		}
	}
	if( ! hit_cell ) {
		console::dump( "Error: gather_neighbor_cells failed.\n" );
		exit(0);
	}
	//
	mirror_samples(*this,p,points,param.accurate_interpolation ? 0.0 : param.MLS_eps);
}
//
bool grid3::reconstruct_fluid( array3<Real> &fluid, std::function<bool(int i, int j, int k, const Real &levelset_value )> test_func ) const {
	//
	scoped_timer timer(this);
	//
	assert(layers.size());
	double dx = layers[0]->dx;
	//
	std::vector<char> T_junction_surface_flags (cell_count);
	if( param.erosion_count ) {
		timer.tick(); console::dump( "Extrapolating erosion cells (%d cells)...", param.erosion_count );
		iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
			if( is_surface_cell(cell_id) && is_T_junction_cell(cell_id)) {
				T_junction_surface_flags[cell_id.index] = 1;
			}
		});
		for( unsigned n=0; n<param.erosion_count-1; ++n ) {
			std::vector<char> prev_T_junction_surface_flags(T_junction_surface_flags);
			layers[0]->active_cells.const_parallel_actives([&]( int i, int j, int k, const auto &it, int thread_index ) {
				const cell_id3 &cell_id = {0,vec3i(i,j,k),it()};
				if( ! prev_T_junction_surface_flags[cell_id.index] ) {
					this->iterate_cell_neighbors(cell_id,[&]( char dim, const cell_id3 &cell_neigh_id ) {
						if( prev_T_junction_surface_flags[cell_neigh_id.index]) {
							T_junction_surface_flags[cell_id.index] = 1;
						}
					});
				}
			});
		}
		console::dump("Done. Took %s\n", timer.stock("extrapolate_erode_cells").c_str());
	}
	//
	uint_type erosion_total_count = std::accumulate(T_junction_surface_flags.begin(),T_junction_surface_flags.end(),0);
	console::dump( "Erosion total count: %u\n", erosion_total_count );
	console::write("num_erosion_count",erosion_total_count);
	//
	timer.tick(); console::dump( "Reconstructing high-res fluid levelset..." );
	fluid.clear();
	fluid.set_as_levelset(3.0*dx);
	bool found_active (false);
	layers[0]->active_cells.const_serial_actives([&](int i, int j, int k, const auto &it) {
		const Real levelset_value = levelset[it()];
		if( test_func(i,j,k,levelset_value) && ! T_junction_surface_flags[it()]) {
			fluid.set(i,j,k,levelset_value);
			found_active = true;
		}
	});
	console::dump( "Done. active=%s, count=%u. Took %s\n", found_active ? "Yes" : "No", fluid.count(), timer.stock("reconstruct_fluid").c_str() );
	return found_active;
}
//
bool grid3::reconstruct_velocity( macarray3<Real> &velocity, std::function<bool(int dim, int i, int j, int k)> test_func ) const {
	//
	scoped_timer timer(this);
	timer.tick(); console::dump( "Reconstructing high-res velocity..." );
	//
	velocity.clear();
	unsigned char found_active [] = {0,0,0};
	parallel.for_each(DIM3,[&]( size_t dim ) {
		layers[0]->active_faces[dim].const_serial_actives([&](int i, int j, int k, const auto &it) {
			if( test_func(dim,i,j,k)) {
				const Real velocity_value = this->velocity[it()];
				velocity[dim].set(i,j,k,velocity_value);
				found_active[dim] = 1;
			}
		});
	});
	const bool result = found_active[0]+found_active[1]+found_active[2];
	console::dump( "Done. active=%s, count=%u. Took %s\n", result ? "Yes" : "No", velocity.count(), timer.stock("reconstruct_velocity").c_str() );
	return result;
}
//
void grid3::compute_cell_map () {
	//
	scoped_timer timer(this);
	timer.tick(); console::dump( "Gathering active liquid cells..." );
	//
	valid_cell_count = 0;
	cell_map.clear();
	cell_map.resize(cell_count);
	//
	iterate_active_cells([&](const cell_id3 &cell_id, int tid ) {
		bool valid (false);
		if( levelset[cell_id.index] < 0.0 ) {
			bool open (false);
			for( char dim : DIMS3 ) {
				this->iterate_face_neighbors(cell_id,dim,[&]( const face_id3 &face_id ){
					if( area[face_id.index] ) open = true;
				});
				if( open ) break;
			}
			if( open ) {
				cell_map[cell_id.index] = 1;
			}
		}
	});
	//
	for( uint_type n=0; n<cell_count; ++n ) {
		if( cell_map[n] ) {
			cell_map[n] = ++ valid_cell_count;
		}
	}
	//
	cell_map.shrink_to_fit();
	//
	console::dump( "Done. Found %u cells. Took %s.\n", valid_cell_count, timer.stock("gather_liquid_cells").c_str());
	console::write("num_gather_liquid_cells",valid_cell_count);
}
//
void grid3::compute_face_map () {
	//
	scoped_timer timer(this);
	timer.tick(); console::dump( "Gathering active velocity faces..." );
	//
	valid_face_count = 0;
	face_map.clear();
	face_map.resize(face_count);
	//
	iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
		bool valid (false);
		this->get_scaled_gradient(face_id,[&]( const cell_id3 &cell_id, double value, const grid3::gradient_info3 &info ) {
			valid = true;
		});
		if( valid ) {
			face_map[face_id.index] = 1;
		}
	});
	//
	for( uint_type n=0; n<face_count; ++n ) {
		if( face_map[n] ) {
			face_map[n] = ++ valid_face_count;
		}
	}
	//
	face_map.shrink_to_fit();
	//
	console::dump( "Done. Found %u faces. Took %s.\n", valid_face_count, timer.stock("gather_liquid_faces").c_str());
	console::write("num_gather_liquid_faces",valid_face_count);
}
//
void grid3::clear_map () {
	//
	valid_cell_count = 0;
	valid_face_count = 0;
	cell_map.clear();
	face_map.clear();
	cell_map.shrink_to_fit();
	face_map.shrink_to_fit();
}
//
double grid3::get_finest_dx() const {
	//
	for( char depth=0; depth<layers.size(); ++depth ) {
		const auto &layer = *layers[depth];
		if( layer.active_cells.count()) {
			return layer.dx;
		}
	}
	return 0.0;
}
//
bool grid3::is_surface_cell( const cell_id3 &cell_id ) const {
	//
	const Real levelset0 = levelset[cell_id.index];
	bool cross_surface (false);
	iterate_cell_neighbors(cell_id,[&]( char dim, const cell_id3 &cell_neigh_id ) {
		if( levelset0 * levelset[cell_neigh_id.index] < 0.0 ) {
			cross_surface = true;
		}
	});
	return cross_surface;
}
//
bool grid3::is_T_junction_cell( const cell_id3 &cell_id ) const {
	//
	char depth (cell_id.depth);
	bool result (false);
	iterate_cell_neighbors(cell_id,[&]( char dim, const cell_id3 &cell_neigh_id ) {
		if( cell_neigh_id.depth != depth ) result = true;
	});
	return result;
}
//
double grid3::get_volume() const {
	//
	// Compute volume
	std::vector<double> volume_threads(parallel.get_thread_num(),0.0);
	iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
		bool rho_added (false);
		get_scaled_gradient(face_id,[&]( const cell_id3 &cell_id, double value, const grid3::gradient_info3 &info ) {
			if( ! rho_added ) {
				volume_threads[tid] += info.area * (info.t_junction ? 0.75 : 1.0) * info.rho * (info.dx * info.dx * info.dx);
				rho_added = true;
			}
		});
	});
	double current_volume (0.0);
	for( auto &e : volume_threads ) current_volume += e / DIM3;
	return current_volume;
}
//
double grid3::get_cell_volume( const cell_id3 &cell_id ) const {
	//
	const uint_type index = cell_id.index;
	const vec3d p = get_cell_position(cell_id);
	const Real phi = levelset[index];
	//
	double volume (0.0);
	if( phi < 0.0 ) for( int dim : DIMS3 ) {
		//
		const auto &layer = *layers[cell_id.depth];
		const auto &index_array = layer.active_cells;
		const auto &fill_flags = layer.fill_flags;
		const int i = cell_id.pi[0];
		const int j = cell_id.pi[1];
		const int k = cell_id.pi[2];
		const double dx = get_cell_dx(cell_id);
		//
		vec3i ijk (i,j,k);
		double phi_backward (0.0), phi_forward (0.0);
		//
		if( ijk[dim] > 0 ) {
			double h (dx);
			vec3i qi(i-(dim==0),j-(dim==1),k-(dim==2));
			if( index_array.active(qi)) {
				phi_backward = levelset[index_array(qi)];
			} else {
				if(fill_flags(qi)) h = 0.75 * dx;
				else h = 1.25 * dx;
				vec3d q = p-h*vec3i(dim==0,dim==1,dim==2);
				phi_backward = sample_levelset(q);
			}
			if( phi_backward > 0.0 ) {
				double f = utility::fraction(phi,phi_backward);
				volume += (dx*dx*h) * f;
			} else {
				volume += 0.5 * (dx*dx*dx);
			}
		} else {
			// Same as above
			volume += 0.5 * (dx*dx*dx);
		}
		//
		if( ijk[dim] < index_array.shape()[dim]-1 ) {
			double h (dx);
			vec3i qi(i+(dim==0),j+(dim==1),k+(dim==2));
			if( index_array.active(qi)) {
				phi_forward = levelset[index_array(qi)];
			} else {
				if(fill_flags(qi)) h = 0.75 * dx;
				else h = 1.25 * dx;
				vec3d q = p+h*vec3i(dim==0,dim==1,dim==2);
				phi_forward = sample_levelset(q);
			}
			if( phi_forward > 0.0 ) {
				double f = utility::fraction(phi,phi_forward);
				volume += (dx*dx*h) * f;
			} else {
				volume += 0.5 * (dx*dx*dx);
			}
		} else {
			// Same as above
			volume += 0.5 * (dx*dx*dx);
		}
	}
	return volume / DIM3;
}
//
vec3d grid3::get_upwind_gradient( const cell_id3 &cell_id, const std::vector<Real> &levelset ) const {
	//
	const Real levelset0 (levelset[cell_id.index]);
	const Real sgn = levelset0 > 0.0 ? 1.0 : -1.0;
	vec3d result;
	//
	for( char dim : DIMS3 ) {
		//
		iterate_face_neighbors(cell_id,dim,[&]( const face_id3 &face_id ){
			//
			unsigned count (0);
			Real sum (0.0), grad (0.0);
			this->get_gradient(face_id,[&]( const cell_id3 &neighbor_cell_id, double value, const gradient_info3 &info ) {
				sum += info.levelset;
				count ++;
				grad += value * levelset[neighbor_cell_id.index];
			});
			if( count ) {
				sum = sum / count;
				if( sgn*sum < sgn*levelset0 ) {
					if( std::abs(grad) > std::abs(result[dim])) result[dim] = grad;
				}
			}
		});
	}
	return result;
}
//
void grid3::draw_grid( graphics_engine &g, double slice_z, bool fill_checkboard ) const {
	//
	if( fill_checkboard ) {
		//
		serial_iterate_active_cells([&]( const cell_id3 &cell_id ) {
			//
			const double dx = get_cell_dx(cell_id);
			const vec3d p = get_cell_position(cell_id);
			if( p[2]-0.5*dx < slice_z && slice_z <= p[2]+0.5*dx ) {
				//
				const int i = cell_id.pi[0];
				const int j = cell_id.pi[1];
				//
				if( (i+j) % 2 ) {
					g.color4(1.0,1.0,1.0,0.25);
				} else {
					g.color4(0.5,0.5,0.5,0.25);
				}
				//
				g.begin(graphics_engine::MODE::TRIANGLE_FAN);
				g.vertex2v((dx*vec2d(i,j)).v);
				g.vertex2v((dx*vec2d(i+1,j)).v);
				g.vertex2v((dx*vec2d(i+1,j+1)).v);
				g.vertex2v((dx*vec2d(i,j+1)).v);
				g.end();
			}
		});
	}
	//
	serial_iterate_active_cells([&]( const cell_id3 &cell_id ) {
		//
		const double dx = get_cell_dx(cell_id);
		const vec3d p = get_cell_position(cell_id);
		if( p[2]-0.5*dx < slice_z && slice_z <= p[2]+0.5*dx ) {
			//
			const int i = cell_id.pi[0];
			const int j = cell_id.pi[1];
			//
			double fgc[4];
			g.get_foreground_color(fgc);
			g.color4(fgc[0],fgc[1],fgc[2],0.25);
			g.begin(graphics_engine::MODE::LINE_LOOP);
			g.vertex2v((dx*vec2d(i,j)).v);
			g.vertex2v((dx*vec2d(i+1,j)).v);
			g.vertex2v((dx*vec2d(i+1,j+1)).v);
			g.vertex2v((dx*vec2d(i,j+1)).v);
			g.end();
		}
	});
}
//