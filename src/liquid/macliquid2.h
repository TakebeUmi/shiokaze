/*
**	macliquid2.h
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
#ifndef SHKZ_MACLIQUID2_H
#define SHKZ_MACLIQUID2_H
//
#include <shiokaze/array/array2.h>
#include <shiokaze/array/macarray2.h>
#include <shiokaze/ui/drawable.h>
#include <shiokaze/utility/gridutility2_interface.h>
#include <shiokaze/advection/macadvection2_interface.h>
#include <shiokaze/utility/macutility2_interface.h>
#include <shiokaze/utility/macstats2_interface.h>
#include <shiokaze/visualizer/macvisualizer2_interface.h>
#include <shiokaze/visualizer/gridvisualizer2_interface.h>
#include <shiokaze/projection/macproject2_interface.h>
#include <shiokaze/surfacetracker/maclevelsetsurfacetracker2_interface.h>
#include <shiokaze/timestepper/timestepper_interface.h>
#include <shiokaze/utility/graphplotter_interface.h>
#include <shiokaze/core/dylibloader.h>
#include <shiokaze/rigidbody/rigidworld2_interface.h>
#include <shiokaze/utility/graphplotter_interface.h>
#include <functional>
//
SHKZ_BEGIN_NAMESPACE
//
class macliquid2 : public drawable {
public:
	//
	macliquid2();
	LONG_NAME("MAC Liquid 2D")
	ARGUMENT_NAME("Liquid")
	//
protected:
	//
	virtual void drag( double x, double y, double z, double u, double v, double w ) override;
	virtual void idle() override;
	virtual void setup_window( std::string &name, int &width, int &height ) const override;
	virtual void draw( graphics_engine &g ) const override;
	virtual bool should_quit() const override { return m_timestepper->should_quit(); }
	virtual bool should_screenshot() const override { return m_timestepper->should_export_frame(); }
	//
	virtual void load( configuration &config ) override;
	virtual void configure( configuration &config ) override;
	virtual void post_initialize( bool initialized_from_file ) override;
	//
	macarray2<Real> m_velocity{this};
	macarray2<Real> m_solid_velocity{this};
	macarray2<Real> m_external_force{this};
	array2<Real> m_fluid{this};
	array2<Real> m_solid{this};
	//
	environment_setter arg_shape{this,"shape",&m_shape};
	environment_setter arg_dx{this,"dx",&m_dx};
	//
	macproject2_driver m_macproject{this,"macpressuresolver2"};
	macadvection2_driver m_macadvection{this,"macadvection2"};
	macsurfacetracker2_driver m_macsurfacetracker{this,"maclevelsetsurfacetracker2"};
	timestepper_driver m_timestepper{this,"timestepper"};
	gridutility2_driver m_gridutility{this,"gridutility2"};
	macutility2_driver m_macutility{this,"macutility2"};
	macstats2_driver m_macstats{this,"macstats2"};
	gridvisualizer2_driver m_gridvisualizer{this,"gridvisualizer2"};
	macvisualizer2_driver m_macvisualizer{this,"macvisualizer2"};
	graphplotter_driver m_graphplotter{this,"graphplotter"};
	dylibloader m_dylib;
	//
	shape2 m_shape;
	double m_dx;
	double m_target_volume;
	bool m_force_exist;
	unsigned m_graph_lists[4];
	//
	virtual void initialize( const filestream &file ) override {
		file.r(m_shape);
		file.r(m_dx);
		file.r(m_force_exist);
		file.r(m_target_volume);
	}
	virtual void serialize( const filestream &file ) const override {
		file.w(m_shape);
		file.w(m_dx);
		file.w(m_force_exist);
		file.w(m_target_volume);
	}
	//
	virtual bool keyboard( int key, int action, int mods ) override {
		if( action == UI_interface::PRESS ) {
			if ( key == UI_interface::KEY_S ) {
				filestream file("macliquid2.dat",filestream::WRITE);
				recursive_serialize(file);
				return true;
			}
		}
		return false;
	}
	//
	std::function<bool( double dx, double dt, double time, unsigned step )> m_check_inject_func;
	std::function<bool( const vec2d &p, double dx, double dt, double time, unsigned step, double &fluid, vec2d &velocity )> m_inject_func;
	std::function<void( double dx, double dt, double time, unsigned step, double &volume_change )> m_post_inject_func;
	std::function<vec2d(double)> m_gravity_func;
	//
	struct Parameters {
		vec2d gravity;
		double surftens_k;
		bool volume_correction;
		bool show_graph;
		double volume_change_tol_ratio;
	};
	//
	Parameters m_param;
	//
	virtual void inject_external_force( macarray2<Real> &velocity, double dt, bool clear=true );
	virtual void inject_external_fluid( array2<Real> &fluid, macarray2<Real> &velocity, double dt, double time );
	virtual size_t do_inject_external_fluid( array2<Real> &fluid, macarray2<Real> &velocity, double dt, double time, unsigned step );
	virtual void set_volume_correction( macproject2_interface *macproject );
	virtual void extend_both( int w=2 );
	virtual void add_to_graph();
};
//
SHKZ_END_NAMESPACE
//
#endif