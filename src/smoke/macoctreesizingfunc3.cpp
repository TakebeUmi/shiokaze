/*
**	macoctreesizingfunc3.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on December 31, 2019. All rights reserved.
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
#include <limits>
#include "macoctreesizingfunc3.h"
#include <shiokaze/core/console.h>
#include <numeric>
#include "fastapprox/fastexp.h"
#include "fastapprox/fastlog.h"
//
SHKZ_USING_NAMESPACE
using namespace macotreeliquid3_namespace;
//
void macoctreesizingfunc3::configure( configuration &config ) {
	//
	configuration::auto_group group(config,*this);
	//
	config.get_double("HalfBandWidth",m_param.halfband_width,"Half bandwidth size");
	config.get_integer("DiffuseCount",m_param.diffuse_count,"Diffusion count");
	config.get_double("CurvatureBase",m_param.curvature_base,"Curvature base");
	config.get_double("VelocityBase",m_param.velocity_base,"Velocity base");
	config.get_double("DecayRate",m_param.decay_rate,"Sizing function decay rate");
	config.get_double("DecayTime",m_param.decay_time,"Sizing function decay time");
	config.get_double("SizingStrength",m_param.sizing_strength,"Sizing function strength");
	config.get_integer("InjectionDepth",m_param.injection_depth,"Injection depth");
	config.get_double("VelocityScale",m_param.velocity_scale,"Velocity scale");
}
//
static void diffuse( const grid3 &grid, std::vector<Real> &result, std::vector<unsigned char> &active_flags, int diffuse_count ) {
	//
	bool has_value (false);
	for( const auto &e : active_flags ) {
		if( e ) has_value = true;
	}
	if( ! has_value ) return;
	//
	while( true ) {
		//
		std::vector<size_t> undone_count(grid.parallel.get_thread_num());
		std::vector<Real> save_values(result);
		std::vector<unsigned char> save_flags(active_flags);
		//
		grid.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
			//
			const double dx = grid.get_cell_dx(cell_id);
			if( save_flags[cell_id.index] < diffuse_count ) {
				//
				double wsum (0.0);
				double value (0.0);
				const Real value_center = save_values[cell_id.index];
				//
				if( save_flags[cell_id.index] ) {
					const double w = dx*dx;
					wsum += w;
					value += w * value_center;
				}
				//
				grid.iterate_cell_neighbors(cell_id,[&]( char dim, const cell_id3 &cell_neigh_id ) {
					if( save_flags[cell_neigh_id.index] ) {
						const double dx = grid.get_cell_dx(cell_neigh_id);
						const double w = dx*dx;
						wsum += w;
						value += w * std::max(value_center,save_values[cell_neigh_id.index]);
					}
				});
				//
				if( wsum ) {
					active_flags[cell_id.index] ++;
					result[cell_id.index] = value / wsum;
				} else {
					if( ! save_flags[cell_id.index] ) undone_count[tid] ++;
				}
			}
		});
		//
		if( ! std::accumulate(undone_count.begin(),undone_count.end(),0)) break;
	}
}
//
static void compute_laplacian( const grid3 &grid, const std::vector<Real> &scalar, std::vector<Real> &result, std::vector<unsigned char> &active_flags, double coef ) {
	//
	// Compute the gradient of the scalar
	std::vector<Real> grad_scalar(grid.face_count);
	//
	grid.iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
		double sum (0.0);
		grid.get_gradient(face_id,[&]( const cell_id3 &cell_id, double value, const grid3::gradient_info3 &info ) {
			sum += scalar[cell_id.index] * value;
		});
		grad_scalar[face_id.index] = sum;
	});
	//
	result.clear();
	result.resize(grid.cell_count);
	//
	std::vector<Real> scalar0(grid.cell_count);
	active_flags.resize(grid.cell_count);
	//
	grid.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
		if( grid.is_surface_cell(cell_id) ) {
			const double dx = grid.get_cell_dx(cell_id);
			const vec3d p = grid.get_cell_position(cell_id);
			const double value = 
				(grid.sample_face(p+0.5*vec3d(dx,0.0,0.0),0,grad_scalar)-grid.sample_face(p-0.5*vec3d(dx,0.0,0.0),0,grad_scalar)+
				 grid.sample_face(p+0.5*vec3d(0.0,dx,0.0),1,grad_scalar)-grid.sample_face(p-0.5*vec3d(0.0,dx,0.0),1,grad_scalar)+
				 grid.sample_face(p+0.5*vec3d(0.0,0.0,dx),2,grad_scalar)-grid.sample_face(p-0.5*vec3d(0.0,0.0,dx),2,grad_scalar)) / dx;
			scalar0[cell_id.index] = std::abs(value);
		} else {
			scalar0[cell_id.index] = 0.0;
		}
	});
	//
	// Interpolate exactly on the surface
	grid.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
		//
		const Real levelset0 = scalar[cell_id.index];
		double sum (0.0);
		double wsum (0.0);
		//
		grid.iterate_cell_neighbors(cell_id,[&]( char dim, const cell_id3 &cell_neigh_id ) {
			const double levelset1 = scalar[cell_neigh_id.index];
			if( levelset0 * levelset1 < 0.0 ) {
				const Real v0 = scalar0[cell_id.index];
				const Real v1 = scalar0[cell_neigh_id.index];
				const double d0 = std::abs(levelset0);
				const double d1 = std::abs(levelset1);
				const double value = (d1*v0+d0*v1) / (d0+d1);
				sum += value;
				wsum += 1.0;
			}
		});
		//
		if( wsum ) {
			result[cell_id.index] = coef * sum / wsum;
			active_flags[cell_id.index] = 1;
		} else {
			active_flags[cell_id.index] = 0;
		}
	});
}
//
void macoctreesizingfunc3::compute_sizing_function( const grid3 &grid0, const grid3 &grid1, double dt,
													std::function<double( const vec3d &p )> solid_func,
													std::function<vec3d( const vec3d &p )> additional_velocity_func ) {
	//
	std::vector<Real> save_laplacian (m_laplacian);
	//
	std::vector<unsigned char> active_flags;
	::compute_laplacian(grid1,grid1.levelset,m_laplacian,active_flags,m_param.curvature_base);
	//
	if( solid_func ) {
		std::vector<Real> solid_levelset(grid1.cell_count);
		std::vector<Real> solid_laplacian(grid1.cell_count);
		grid1.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
			solid_levelset[cell_id.index] = solid_func(grid1.get_cell_position(cell_id))-0.5*grid1.get_cell_dx(cell_id);
		});
		std::vector<unsigned char> active_flags1;
		::compute_laplacian(grid1,solid_levelset,solid_laplacian,active_flags1,m_param.curvature_base);
		for( size_t n=0; n<grid1.cell_count; ++n ) {
			if( active_flags[n] && active_flags1[n] ) m_laplacian[n] += solid_laplacian[n];
		}
	}
	//
	grid1.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
		if( active_flags[cell_id.index] ) {
			vec3d v;
			int count (0);
			const double dx = grid1.get_cell_dx(cell_id);
			grid1.get_unmofidied_divergence(cell_id,[&]( const face_id3 &face_id, double value ) {
				double u = grid1.velocity[face_id.index];
				if( additional_velocity_func ) {
					u += additional_velocity_func(grid1.get_face_position(face_id))[face_id.dim];
				}
				v[face_id.dim] += value * u / (dx*dx*dx);
				count ++;
			});
			if( count == 6 ) {
				m_laplacian[cell_id.index] += v.len() * m_param.velocity_base * m_param.velocity_scale;
			}
		}
	});
	//
	::diffuse(grid1,m_laplacian,active_flags,m_param.diffuse_count);
	//
	if( dt && m_param.decay_rate && grid0.cell_count && grid0.cell_count == save_laplacian.size()) {
		const double r = pow(m_param.decay_rate,dt/m_param.decay_time);
		grid1.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
			const Real v = m_laplacian[cell_id.index];
			const vec3d p = grid1.get_cell_position(cell_id);
			const vec3d u = grid0.sample_velocity(p);
			const Real w = grid0.sample_cell(p-dt*u,save_laplacian);
			m_laplacian[cell_id.index] = std::max(v,(Real)r*w);
		});
	}
}
//
void macoctreesizingfunc3::activate_cells( const grid3 &grid0, grid3 &grid1,
										   std::function<double( const vec3d &p )> solid_func,
										   std::function<double( const vec3d &p )> additional_fluid_func ) const {
	//
	assert(m_laplacian.size() == grid0.cell_count);
	//
	if( ! solid_func ) solid_func = []( const vec3d &p ){ return 1.0; };
	if( ! additional_fluid_func ) additional_fluid_func = []( const vec3d &p ){ return 1.0; };
	//
	const double dx0 = grid1.layers[0]->dx;
	grid1.activate_cells([&]( char depth, const vec3d &p ) {
		//
		bool result (false);
		const double dx1 = grid1.layers[depth]->dx;
		//
		if( solid_func(p) > -dx1 ) {
			//
			double dx;
			if( m_param.sizing_strength > 0.0 && m_param.sizing_strength < 1.0 ) {
				dx = fasterpow2(fastlog(dx1/dx0)/fastlog(1.0+m_param.sizing_strength)) * dx0;
			} else if( m_param.sizing_strength <= 0.0 ) {
				dx = std::numeric_limits<double>::max();
			} else {
				dx = dx1;
			}
			const double ext = additional_fluid_func(p);
			const double levelset = std::min(grid0.sample_levelset(p),ext);
			const double d = std::abs(levelset);
			const bool pass = levelset > 0.0 ? d <= m_param.halfband_width*dx1 : d <= m_param.halfband_width*dx;
			//
			if( pass ) {
				//
				if( ext < m_param.halfband_width*dx ) {
					result = depth >= m_param.injection_depth;
				} else {
					result = grid0.sample_cell(p,m_laplacian) >= 1.0/dx;
				}
			}
		}
		//
		return result;
	});
}
//