/*
**	macoctreeproject3.h
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
#ifndef SHKZ_OCTREEPROJECT3_H
#define SHKZ_OCTREEPROJECT3_H
//
#include "macoctreegrid3.h"
#include <shiokaze/math/RCMatrix_interface.h>
#include <shiokaze/linsolver/RCMatrix_solver.h>
//
SHKZ_BEGIN_NAMESPACE
//
namespace macotreeliquid3_namespace {
	//
	class macoctreeproject3 : public recursive_configurable, public credit {
	public:
		//
		LONG_NAME("MAC Octee Project 3D")
		ARGUMENT_NAME("macoctreeproject")
		//
		virtual void configure( configuration &config ) override;
		virtual void post_initialize( bool initialized_from_file ) override;
		//
		macoctreeproject3 ( recursive_configurable *parent ) {
			if( parent ) parent->add_child(this);
			else setup_now();
		}
		//
		struct matrix3 {
			RCMatrix_ptr<size_t,double> Lhs, G, D;
			bool allocated {false};
			bool assembled {false};
		};
		matrix3 m_matrix;
		//
		struct Parameters {
			bool volume_correction {true};
			bool volume_correct_skip_surfaces {false};
			bool check_symmetric {true};
			bool check_positive_diag {true};
			bool fix_divergence {false};
			bool debug_assemble {false};
			double volume_recover_ratio {0.9};
			double minimal_volume_correct_timestep {0.01};
			double np {25.0};
			bool remove_one_degrees_of_freedom {false};
			int special_boundary_condition {0};
		};
		Parameters m_param;
		//
		void assemble_matrix( grid3 &grid );
		void clear_matrix();
		void set_moving_solid( std::function<double(const vec3d &p)> moving_solid_func );
		void project( grid3 &grid, double dt, std::function<vec3d(const vec3d &p)> solid_velocity=nullptr, std::vector<Real> *pressure=nullptr );
		void project( grid3 &grid, double dt,
					  size_t region_count,
					  const std::vector<uint_type> &regions,
					  const std::vector<Real> &current_volumes,
					  const std::vector<Real> &target_volumes,
					  std::vector<Real> &y_list,
					  std::function<vec3d(const vec3d &p)> solid_velocity=nullptr,
					  std::vector<Real> *pressure=nullptr );
		//
		double m_initial_volume, m_y_prev;
		parallel_driver m_parallel{this};
		RCMatrix_factory_driver<size_t,double> m_factory{this,"RCMatrix"};
		RCMatrix_solver_driver<size_t,double> m_solver{this,"amg"};
		std::function<double(const vec3d &p)> m_moving_solid_func {nullptr};
		//
		virtual void initialize( const filestream &file ) override {
			file.r(m_initial_volume);
			file.r(m_y_prev);
		}
		virtual void serialize( const filestream &file ) const override {
			file.w(m_initial_volume);
			file.w(m_y_prev);
		}
	};
};
//
SHKZ_END_NAMESPACE
//
#endif
//