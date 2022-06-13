/*
**	macbackwardflip3.h
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on April 9, 2017.
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
#ifndef SHKZ_BACKWARDFLIP3_H
#define SHKZ_BACKWARDFLIP3_H
//
#include <shiokaze/backwardflip/macbackwardflip3_interface.h>
#include <shiokaze/parallel/parallel_driver.h>
#include <deque>
#include <functional>
#include <memory>
//
SHKZ_BEGIN_NAMESPACE
//
class macbackwardflip3 : public macbackwardflip3_interface {
protected:
	//
	virtual void configure( configuration &config ) override;
	virtual void initialize( const shape3 &shape, double dx ) override;
	virtual void post_initialize( bool initialized_from_file ) override;
	//
	virtual bool backtrace( const array3<Real> &solid, const array3<Real> &fluid ) override;
	virtual bool fetch( macarray3<Real> &u_reconstructed ) const override;
	virtual bool fetch( array3<Real> &density_reconstructed ) const override;
	//
	virtual void register_buffer(
						const macarray3<Real> &u1,				// Velocity at the end of the step
						const macarray3<Real> &u0,				// Velocity at the beggining of the step
						const macarray3<Real> *u_reconstructed,	// Reconstructed dirty velocity of the beggining of the step - can be nullptr
						const macarray3<Real> *g,					// Pressure gradient and the external forces (scaled by dt) - can be nullptr
						const array3<Real> *d1,					// Density field of the end of the step - can be nullptr
						const array3<Real> *d0,					// Density field of the beggining of the step - can be nullptr
						const array3<Real> *d_added,				// Density field source of the current step - can be nullptr
						double dt ) override;						// Time-step size of the current step
	//
	virtual void draw( graphics_engine &g ) const override;
	//
	struct Parameters {
		unsigned max_layers {4};
		unsigned max_velocity_layers {4};
		unsigned r_sample {2};
		double decay_rate {0.9};
		double decay_truncate {1e-2};
		bool use_hachisuka {false};
		bool use_temporal_adaptivity {false};
		bool use_accumulative_buffer {true};
		bool use_spatial_adaptivity {true};
		unsigned max_temporal_adaptivity_level {6};
		double temporal_adaptive_rate {0.75};
		double spatial_adaptive_rate {0.5};
		double spatial_density_threshold {0.01};
		double inject_diff {0.9};
	};
	Parameters m_param;
	//
	typedef struct {
		//
		std::shared_ptr<macarray3<Real> > u;
		std::shared_ptr<macarray3<Real> > u_reconstructed;
		std::shared_ptr<macarray3<Real> > g;
		std::shared_ptr<array3<Real> > d;
		std::shared_ptr<array3<Real> > d_added;
		//
		double dt;
		double time;
		bool allocated {false};
		//
		void allocate () {
			if( ! allocated ) {
				u = std::make_shared<macarray3<Real> >();
				u_reconstructed = std::make_shared<macarray3<Real> >();
				g = std::make_shared<macarray3<Real> >();
				d = std::make_shared<array3<Real> >();
				d_added = std::make_shared<array3<Real> >();
				allocated = true;
			}
		}
		//
	} layer3;
	//
	typedef struct {
		std::vector<vec3r> p;
		std::vector<vec3r> u;
		std::vector<Real> mass;
		std::vector<std::vector<Real> > adaptivity_rate;
		std::vector<Real> s;
	} tracers3;
	tracers3 m_tracer;
	//
	typedef struct {
		std::vector<Real> wsum;
		std::vector<vec3r> vel;
		std::vector<vec3r> g;
	} accumulator3;
	accumulator3 m_accumulator;
	//
	macarray3<Real> m_u_reconstructed{this};
	array3<Real> m_density_reconstructed{this};
	//
	bool m_exist_gradient;
	bool m_exist_density;
	//
	unsigned m_step_back_limit {0};			// Followed by the method of Hachisuka [H06]
	array3<vec3r> m_forward_tracers{this};	// Followed by the method of Hachisuka [H06]
	array3<vec3r> m_g_integrated{this};		// Followed by the method of Hachisuka [H06]
	//
	virtual void reset_forward_tracers();
	virtual void integrate_forward_tracers( const macarray3<Real> &velocity0, const macarray3<Real> &velocity1, const macarray3<Real> &g, double dt );
	//
	std::deque<layer3> m_buffers;
	layer3 m_back_buffer;
	//
	std::vector<std::deque<layer3> > m_coarse_buffers;
	std::vector<unsigned> m_level_stored;
	array3<char> m_spatial_adaptivity{this};
	//
	shape3 m_shape;
	double m_dx;
	unsigned m_step {0};
	macarray3<Real> m_velocity{this};
	array3<Real> m_density{this};
	macarray3<Real> m_u_diff{this};
	//
	std::vector<vec3r> m_original_seed_vector;
	std::vector<Real> m_original_seed_mass;
	array3<std::vector<unsigned> > m_seed_cell{this};
	macarray3<std::vector<unsigned> > m_seed_face{this};
	parallel_driver m_parallel{this};
	//
	void backtrace(	std::vector<vec3r> &p, std::vector<vec3r> &u, const std::vector<Real> &mass, std::vector<std::vector<Real> > &adaptivity_rate, std::vector<Real> *d );
};
//
SHKZ_END_NAMESPACE
//
#endif
