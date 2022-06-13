/*
**	macstats3.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on July 21, 2017.
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
#include <shiokaze/utility/macstats3_interface.h>
#include <shiokaze/utility/macutility3_interface.h>
#include <shiokaze/array/array_utility3.h>
#include <shiokaze/array/array_interpolator3.h>
#include <shiokaze/core/console.h>
#include <shiokaze/core/timer.h>
#include <algorithm>
//
SHKZ_USING_NAMESPACE
using namespace array_utility3;
using namespace array_interpolator3;
//
class macstats3 : public macstats3_interface {
protected:
	//
	virtual void dump_stats( const array3<Real> &solid, const array3<Real> &fluid, const macarray3<Real> &velocity, const timestepper_interface *tmstepper ) const override {
		//
		global_timer::pause();
		//
		// Count the number of active cells
		unsigned num_active_fluid (0);
		if( levelset_exist(solid)) {
			fluid.const_serial_actives([&](int i, int j, int k, const auto &it) {
				if( it() < 0.0 && interpolate<Real>(solid,vec3i(i,j,k).cell()) > 0.0 ) num_active_fluid ++;
			});
		} else {
			fluid.const_serial_inside([&](const auto &it) {
				num_active_fluid ++;
			});
		}
		console::write("macstats3_number_active_cells", num_active_fluid);
		//
		// Measure the kinetic energy
		if( m_param.report_kinetic_energy ) {
			double kinetic_energy = m_macutility->get_kinetic_energy(solid,fluid,velocity);
			if( m_param.report_console ) {
				if( has_different_values(fluid)) console::dump( "Report: active fluid cells = %d, kinetic energy = %.3e\n", num_active_fluid, kinetic_energy );
				else console::dump( "Report: kinetic energy = %.3e\n", kinetic_energy );
			}
			console::write("macstats3_kinetic_energy", kinetic_energy);
		}
		//
		global_timer::resume();
	}
	//
	virtual void configure( configuration &config ) override {
		config.get_bool("ReportConsole",m_param.report_console,"Whether to report in console");
		config.get_bool("ReportKineticEnergy",m_param.report_console,"Whether to report kinetic energy");
	}
	//
	virtual void initialize( const shape3 &shape, double dx ) override {
		m_shape = shape;
		m_dx = dx;
	}
	virtual void initialize( const filestream &file ) override {
		file.r(m_shape);
		file.r(m_dx);
	}
	virtual void serialize( const filestream &file ) const override {
		file.w(m_shape);
		file.w(m_dx);
	}
	//
	struct Parameters {
		bool report_console {true};
		bool report_kinetic_energy {true};
	};
	Parameters m_param;
	//
	macutility3_driver m_macutility{this,"macutility3"};
	//
	shape3 m_shape;
	double m_dx;
};
//
extern "C" module * create_instance() {
	return new macstats3();
}
//
extern "C" const char *license() {
	return "MIT";
}
//