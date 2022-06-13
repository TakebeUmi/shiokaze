/*
**	macflipsmoke2.h
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on Aug 22, 2017.
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
#ifndef SHKZ_MACFLIPSMOKE2_H
#define SHKZ_MACFLIPSMOKE2_H
//
#include "macsmoke2.h"
#include <shiokaze/flip/macflip2_interface.h>
//
SHKZ_BEGIN_NAMESPACE
//
class macflipsmoke2 : public macsmoke2 {
public:
	//
	macflipsmoke2();
	LONG_NAME("MAC FLIP Smoke 2D")
	//
protected:
	//
	virtual void idle() override;
	virtual void draw( graphics_engine &g ) const override;
	//
	virtual void configure( configuration &config ) override;
	virtual void post_initialize( bool initialized_from_file ) override;
	//
	double interpolate_solid( const vec2d &p ) const;
	vec2d interpolate_velocity( const vec2d &p ) const;
	//
	macflip2_driver m_flip{this,"macnbflip2"};
	//
	struct Parameters {
		double PICFLIP;
		double gridmass;
	};
	Parameters m_param;
	//
};
//
SHKZ_END_NAMESPACE
//
#endif