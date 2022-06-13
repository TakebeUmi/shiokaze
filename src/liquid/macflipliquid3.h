/*
**	macflipliquid3.h
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on June 2, 2017.
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
#ifndef SHKZ_MACFLIPLIQUID3_H
#define SHKZ_MACFLIPLIQUID3_H
//
#include "macliquid3.h"
#include <shiokaze/flip/macflip3_interface.h>
#include <shiokaze/particlerasterizer/particlerasterizer3_interface.h>
#include <shiokaze/redistancer/redistancer3_interface.h>
//
SHKZ_BEGIN_NAMESPACE
//
class macflipliquid3 : public macliquid3 {
public:
	//
	macflipliquid3();
	LONG_NAME("MAC FLIP Liquid 3D")
	//
protected:
	//
	virtual size_t do_inject_external_fluid( array3<Real> &fluid, macarray3<Real> &velocity, double dt, double time, unsigned step ) override;
	virtual void idle() override;
	virtual void draw( graphics_engine &g ) const override;
	//
	virtual void configure( configuration &config ) override;
	virtual void post_initialize( bool initialized_from_file ) override;
	virtual void do_export_mesh( unsigned frame ) const override;
	virtual void render_mesh( unsigned frame ) const override;
	//
	struct Parameters {
		double PICFLIP;
		bool disable_resample {false};
		unsigned levelset_half_bandwidth_count {3};
		bool preupdate_FLIP {true};
	};
	Parameters m_param;
	//
	macflip3_driver m_flip{this,"macexnbflip3"};
	redistancer3_driver m_redistancer{this,"pderedistancer3"};
	//
	double interpolate_fluid( const vec3d &p ) const;
	double interpolate_solid( const vec3d &p ) const;
	vec3d interpolate_velocity( const vec3d &p ) const;
};
//
SHKZ_END_NAMESPACE
//
#endif