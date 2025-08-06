/*
**	macoctreesegregator2.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on December 8, 2019. All rights reserved.
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
#include "macoctreesegregator2.h"
#include <shiokaze/utility/utility.h>
#include <shiokaze/image/color.h>
#include <shiokaze/core/console.h>
#include <shiokaze/graphics/graphics_utility.h>
#include <algorithm>
#include <numeric>
#include <map>
#include <set>
#include <cassert>
//
SHKZ_USING_NAMESPACE
using namespace macotreeliquid2_namespace;
//
void macoctreesegregator2::configure( configuration &config ) {
	//
	configuration::auto_group group(config,*this);
	config.get_double("RatioLimit",m_param.ratio_limit,"Volume change maximal ratio");
	config.get_bool("FormMatrix",m_param.form_matrix,"Form matrix for segregation computation");
}
//
size_t macoctreesegregator2::segregate( const grid2 &grid, std::vector<uint_type> &regions ) const {
	//
	regions.resize(grid.cell_count);
	grid.iterate_active_cells([&]( const cell_id2 &cell_id, int tid ) {
		if( grid.levelset[cell_id.index] < 0.0 ) {
			regions[cell_id.index] = cell_id.index+1;
		} else {
			regions[cell_id.index] = 0;
		}
	});
	//
	std::vector<std::vector<uint_type> > connections;
	if( m_param.form_matrix ) {
		connections.resize(grid.cell_count);
		grid.iterate_active_cells([&]( const cell_id2 &cell_id, int tid ) {
			grid.iterate_cell_neighbors(cell_id,[&]( char dim, const cell_id2& neighbor_cell_id ) {
				if( regions[neighbor_cell_id.index] ) {
					connections[cell_id.index].push_back(neighbor_cell_id.index);
				}
			});
		});
	}
	//
	auto iterate_neighbors = [&]( const cell_id2 &cell_id, std::function<void(const uint_type& neighbor_index)> func ) {
		if( m_param.form_matrix ) {
			for( const auto &index : connections[cell_id.index] ) {
				func(index);
			}
		} else {
			grid.iterate_cell_neighbors(cell_id,[&]( char dim, const cell_id2& neighbor_cell_id ) {
				func(neighbor_cell_id.index);
			});
		}
	};
	//
	size_t iterations (0);
	while(true) {
		//
		std::vector<uint_type> regions_save(regions);
		std::vector<uint_type> changed_bucket(grid.parallel.get_thread_num());
		grid.iterate_active_cells([&]( const cell_id2 &cell_id, int tid ) {
			uint_type index = regions_save[cell_id.index];
			if( index ) {
				iterate_neighbors(cell_id,[&]( const uint_type& neighbor_index ) {
					if( regions_save[neighbor_index] ) {
						index = std::min(index,regions_save[neighbor_index]);
					}
				});
				if( regions_save[cell_id.index] != index ) {
					regions[cell_id.index] = index;
					changed_bucket[tid] ++;
				}
			}
		});
		//
		++ iterations;
		if( ! std::accumulate(changed_bucket.begin(),changed_bucket.end(),0)) {
			break;
		}
	}
	//
	uint_type region_count (0);
	uint_type index_cache (0), map_cache (0);
	std::map<uint_type,uint_type> index_map;
	grid.serial_iterate_active_cells([&]( const cell_id2 &cell_id ) {
		uint_type index = regions[cell_id.index];
		if( index ) {
			if( index_cache == index ) {
				regions[cell_id.index] = map_cache;
			} else {
				const auto it = index_map.find(index);
				if( it != index_map.end()) {
					uint_type region_index = it->second;
					regions[cell_id.index] = region_index;
					index_cache = cell_id.index;
					map_cache = region_index;
				} else {
					index_map[index] = ++ region_count;
					regions[cell_id.index] = region_count;
				}
			}
		}
	});
	//
	return region_count;
}
//
void macoctreesegregator2::extrapolate( const grid2 &grid, std::vector<uint_type> &regions ) const {
	//
	// Sort by distance
	std::vector<uint_type> order_map(grid.cell_count);
	std::iota(order_map.begin(),order_map.end(),0);
	std::sort(order_map.begin(),order_map.end(),[&](uint_type a, uint_type b){
		return std::abs(grid.levelset[a]) < std::abs(grid.levelset[b]);
	});
	//
	// Record cell_id
	std::vector<cell_id2> cell_ids(grid.cell_count);
	grid.serial_iterate_active_cells([&]( const cell_id2 &cell_id ) {
		cell_ids[cell_id.index] = cell_id;
	});
	//
	// Propagate
	for( uint_type _n=0; _n<order_map.size(); ++_n ) {
		uint_type n = order_map[_n];
		if( ! regions[n] ) {
			uint_type new_index (0);
			grid.iterate_cell_neighbors(cell_ids[n],[&]( char dim, const cell_id2& neighbor_cell_id ){
				if( regions[neighbor_cell_id.index] ) {
					new_index = regions[neighbor_cell_id.index];
				}
			});
			regions[n] = new_index;
		}
	}
	//
	// Fill the gap
	while(true) {
		bool found (false);
		grid.serial_iterate_active_cells([&]( const cell_id2 &cell_id ) {
			if( ! regions[cell_id.index]) {
				uint_type new_index (0);
				grid.iterate_cell_neighbors(cell_id,[&]( char dim, const cell_id2& neighbor_cell_id ){
					if( regions[neighbor_cell_id.index] ) {
						new_index = regions[neighbor_cell_id.index];
					}
				});
				if( new_index ) regions[cell_id.index] = new_index;
				else found = true;
			}
		});
		if( ! found ) break;
	}
	//
	// Assertion check
	grid.serial_iterate_active_cells([&]( const cell_id2 &cell_id ) {
		assert(regions[cell_id.index]);
	});
}
//
void macoctreesegregator2::extrapolate_jacobi( const grid2 &grid, std::vector<uint_type> &regions ) const {
	//
	while(true) {
		//
		std::vector<unsigned char> found_bucket(grid.parallel.get_thread_num());
		std::vector<uint_type> regions_save(regions);
		grid.iterate_active_cells([&]( const cell_id2 &cell_id, int tid ) {
			bool found (false);
			if( ! regions_save[cell_id.index]) {
				uint_type new_index (0);
				grid.iterate_cell_neighbors(cell_id,[&]( char dim, const cell_id2& neighbor_cell_id ){
					if( regions_save[neighbor_cell_id.index] ) {
						new_index = regions_save[neighbor_cell_id.index];
					}
				});
				if( new_index ) {
					regions[cell_id.index] = new_index;
				} else {
					if( grid.levelset[cell_id.index] < 0.0 ) found = true;
				}
			}
			found_bucket[tid] = found_bucket[tid] || found;
		});
		if( ! std::accumulate(found_bucket.begin(),found_bucket.end(),0) ) break;
	}
	//
	// Assertion check
	grid.iterate_active_cells([&]( const cell_id2 &cell_id, int tid ) {
		if( grid.levelset[cell_id.index] < 0.0 ) {
			assert(regions[cell_id.index]);
		}
	});
}
//
void macoctreesegregator2::backtrace( const grid2 &grid0, const grid2 &grid1,
						double dt, std::function<vec2d( const vec2d &p )> velocity,
						const std::vector<uint_type> &regions0, std::vector<uint_type> &regions1 ) const {
	//
	regions1.resize(grid1.cell_count);
	grid1.iterate_active_cells([&]( const cell_id2 &cell_id, int tid ) {
		vec2d p = grid1.get_cell_position(cell_id);
		p -= dt * velocity(p);
		for( char depth=0; depth < grid0.layers.size(); ++depth ) {
			const auto &layer = *grid0.layers[depth];
			const auto pi = layer.shape.find_cell(p/layer.dx);
			if( layer.active_cells.active(pi)) {
				regions1[cell_id.index] = regions0[layer.active_cells(pi)];
				break;
			}
		}
	});
}
//
void macoctreesegregator2::prune( const grid2 &grid, std::vector<uint_type> &regions ) const {
	//
	grid.iterate_active_cells([&]( const cell_id2 &cell_id, int tid ) {
		if( grid.levelset[cell_id.index] > 0.0 ) regions[cell_id.index] = 0;
	});
}
//
static uint_type compute_max_index( const grid2 &grid, const std::vector<uint_type> &regions ) {
	uint_type max_index (0);
	grid.serial_iterate_active_cells([&]( const cell_id2 &cell_id ) {
		const uint_type index = cell_id.index;
		if( grid.levelset[index] < 0.0 ) {
			max_index = std::max(max_index,regions[index]);
		}
	});
	return max_index;
}
//
void macoctreesegregator2::compute( const grid2 &grid, const std::vector<uint_type> &regions, std::vector<Real> &volumes ) const {
	//
	uint_type max_index = compute_max_index(grid,regions);
	unsigned size = max_index;
	//
	std::vector<std::vector<Real> > volume_bucket(size);
	for( auto &e : volume_bucket ) e.resize(grid.parallel.get_thread_num());
	//
	grid.iterate_active_cells([&]( const cell_id2 &cell_id, int tid ) {
		if( grid.levelset[cell_id.index] < 0.0 ) {
			assert(regions[cell_id.index]);
			volume_bucket[regions[cell_id.index]-1][tid] += grid.get_cell_volume(cell_id);
		}
	});
	//
	volumes.resize(size);
	for( int n=0; n<size; ++n ) {
		volumes[n] = std::accumulate(volume_bucket[n].begin(),volume_bucket[n].end(),0.0);
	}
}
//
double macoctreesegregator2::compute_target_volume( const grid2 &grid1,
									const std::vector<uint_type> &regions01, const std::vector<Real> &volumes0,
									const std::vector<uint_type> &regions11, const std::vector<Real> &volumes1,
									const std::vector<Real> &y_lists0,
									std::vector<Real> &current_volumes1,
									std::vector<Real> &target_volumes1,
									std::vector<Real> &new_y_list1 ) const {
	//
	grid1.iterate_active_cells([&]( const cell_id2 &cell_id, int tid ) {
		if( grid1.levelset[cell_id.index] < 0.0 ) {
			assert(regions01[cell_id.index]);
			assert(regions11[cell_id.index]);
		}
	});
	//
	uint_type max_index0 = compute_max_index(grid1,regions01);
	uint_type max_index1 = compute_max_index(grid1,regions11);
	//
	// Compute the volume change ratio
	std::vector<Real> current_volumes0(max_index0);
	compute(grid1,regions01,current_volumes0);
	compute(grid1,regions11,current_volumes1);
	//
	std::vector<Real> volume_change_ratio(max_index0);
	for( unsigned n=0; n<max_index0; ++n ) {
		const double r = current_volumes0[n] ? volumes0[n] / current_volumes0[n] : 1.0;
		if( m_param.ratio_limit ) {
			volume_change_ratio[n] = std::min(m_param.ratio_limit,std::max(-m_param.ratio_limit,r));
		} else {
			volume_change_ratio[n] = r;
		}
	}
	//
	std::vector<std::vector<std::vector<Real> > > cell_volume_region1_bucket (grid1.parallel.get_thread_num());
	for( auto &e0 : cell_volume_region1_bucket ) {
		e0.resize(max_index1);
		for( auto &e1 : e0 ) {
			e1.resize(max_index0);
		}
	}
	//
	grid1.iterate_active_cells([&]( const cell_id2 &cell_id, int tid ) {
		//
		const uint_type index1 = cell_id.index;
		if( grid1.levelset[index1] < 0.0 ) {
			//
			const double ratio = volume_change_ratio[regions01[index1]-1];
			const double v = grid1.get_cell_volume(cell_id);
			//
			cell_volume_region1_bucket[tid][regions11[index1]-1][regions01[index1]-1] += v * ratio;
		}
	});
	//
	double volume_sum (0.0);
	std::vector<std::vector<Real> > cell_volume_region1 (max_index1);
	for( auto &e : cell_volume_region1 ) e.resize(max_index0);
	for( const auto &e : cell_volume_region1_bucket ) {
		for( unsigned i=0; i<max_index1; ++i ) for( unsigned j=0; j<max_index0; ++j ) {
			cell_volume_region1[i][j] += e[i][j];
			volume_sum += e[i][j];
		}
	}
	console::dump( "volume_sum = %.2e\n", volume_sum);
	//
	target_volumes1.resize(max_index1);
	for( unsigned n=0; n<max_index1; ++n ) {
		target_volumes1[n] = std::accumulate(
			cell_volume_region1[n].begin(),
			cell_volume_region1[n].end(),0.0);
	}
	//
	new_y_list1.resize(max_index1);
	for( unsigned n=0; n<max_index1; ++n ) {
		new_y_list1[n] = 0.0;
		double D = std::accumulate(cell_volume_region1[n].begin(),cell_volume_region1[n].end(),0.0);
		if( D ) {
			for( unsigned m=0; m<max_index0; ++m ) {
				new_y_list1[n] += y_lists0[m] * cell_volume_region1[n][m] / D;
			}
		}
	}
	return volume_sum;
}
//
void macoctreesegregator2::draw_region( graphics_engine &g, const grid2 &grid, const std::vector<uint_type> &regions ) const {
	//
	if( regions.size() == grid.cell_count) {
		//
		uint_type max_index = compute_max_index(grid,regions);
		grid.serial_iterate_active_cells([&]( const cell_id2 &cell_id ) {
			uint_type index = regions[cell_id.index];
			if( index ) {
				//
				double rgb[4];
				color::heatcolor(index/(double)max_index,rgb); rgb[3] = 0.5;
				const int i = cell_id.pi[0];
				const int j = cell_id.pi[1];
				const double dx = grid.get_cell_dx(cell_id);
				//
				g.color4v(rgb);
				g.begin(graphics_engine::MODE::TRIANGLE_FAN);
				g.vertex2v((dx*vec2d(i,j)).v);
				g.vertex2v((dx*vec2d(i+1,j)).v);
				g.vertex2v((dx*vec2d(i+1,j+1)).v);
				g.vertex2v((dx*vec2d(i,j+1)).v);
				g.end();
			}
		});
	}
}
//