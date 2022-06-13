/*
**	movingsolid2.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on December 12, 2019.
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
#include <shiokaze/math/vec.h>
#include <shiokaze/core/configuration.h>
#include <shiokaze/graphics/graphics_engine.h>
#include <shiokaze/graphics/graphics_utility.h>
#include <string>
#include <cmath>
#include <tuple>
//
SHKZ_USING_NAMESPACE
//
static double g_water_level (0.245);
static double g_speed (4.0);
static double g_amplitude (0.25);
static vec2d g_center (0.5,0.25);
static double g_r (0.1);
static double g_flux_flow (0.0);
//
extern "C" void configure( configuration &config ) {
	//
	configuration::auto_group group(config,"Moving Solid Scene 2D","MovingSolid");
	config.get_double("WaterLevel",g_water_level,"Water level");
	config.get_double("Speed",g_speed,"Speed of obstacle");
	config.get_double("Amplitude",g_amplitude,"Speed of obstacle");
	config.get_double("Radius",g_r,"Radius");
	config.get_vec2d("Center",g_center.v,"Center position");
	config.get_double("Flux",g_flux_flow,"Flux velocity on walls");
}
//
extern "C" std::pair<double,vec2d> moving_solid( double time, const vec2d &p ) {
	//
	const vec2d center = g_center+vec2d(-g_amplitude*cos(g_speed*time),0.0);
	const double d = (center-p).len()-g_r;
	//
	vec2d u;
	if( d <= 0.0 ) {
		u = vec2d(g_amplitude*g_speed*sin(g_speed*time),0.0);
	}
	return std::make_pair(d,u);
}
//
extern "C" vec2d velocity( const vec2d &p ) {
	return vec2d(g_flux_flow,0.0);
}
//
extern "C" void set_boundary_flux( double time, Real flux[DIM2][2] ) {
	//
	flux[0][0] = g_flux_flow;
	flux[0][1] = g_flux_flow;
}
//
extern "C" double fluid( const vec2d &p ) {
	return p[1]-g_water_level;
}
//
extern "C" void draw( graphics_engine &g, double time ) {
	//
	const vec2d center = g_center+vec2d(-g_amplitude*cos(g_speed*time),0.0);
	//
	g.color4(0.5,0.5,0.4,1.0);
	graphics_utility::draw_circle(g,center.v,g_r,graphics_engine::MODE::TRIANGLE_FAN);
	//
	g.color4(1.0,1.0,1.0,1.0);
	graphics_utility::draw_circle(g,center.v,g_r,graphics_engine::MODE::LINE_LOOP);
}
//
extern "C" const char *license() {
	return "MIT";
}
//