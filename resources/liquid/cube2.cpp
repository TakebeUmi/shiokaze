/*
**	cube2.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on October 29, 2019.
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
#include <shiokaze/utility/utility.h>
#include <shiokaze/core/configuration.h>
#include <string>
#include <cmath>
//
SHKZ_USING_NAMESPACE
//
static double g_width (0.2);
static unsigned g_default_gn (64);
static vec2d g_center(0.5,0.5);
static int g_version (0);
//
extern "C" void configure( configuration &config ) {
	configuration::auto_group group(config,"Cube Scene 2D","Cube");
	config.get_double("Width",g_width,"Width of cube");
	config.get_vec2d("Center",g_center.v,"Center of cube");
	config.get_integer("Version",g_version,"Version");
}
//
extern "C" std::map<std::string,std::string> get_default_parameters() {
	std::map<std::string,std::string> dictionary;
	dictionary["ResolutionX"] = std::to_string(g_default_gn);
	dictionary["ResolutionY"] = std::to_string(g_default_gn);
	dictionary["Gravity"] = "0.0,0.0";
	dictionary["SurfaceTension"] = "5e-3";
	dictionary["RegionalVolumeCorrection"] = "Yes";
	dictionary["CFL"] = "80";
	dictionary["UseSquaredCFL"] = "Yes";
	return dictionary;
}
//
extern "C" double fluid( const vec2d &p ) {
	if( g_version == 0 ) {
		return utility::box(p,g_center-vec2d(g_width,g_width),g_center+vec2d(g_width,g_width));
	} else if( g_version == 1 ) {
		vec2d center0 = g_center-vec2d(g_width,0.0);
		vec2d center1 = g_center+vec2d(g_width,0.0);
		double w = 0.5 * g_width;
		double value (1.0);
		value = std::min(value,utility::box(p,center0-vec2d(w,w),center0+vec2d(w,w)));
		value = std::min(value,utility::box(p,center1-vec2d(w,w),center1+vec2d(w,w)));
		return value;
	} else {
		return 1.0;
	}
}
//
extern "C" vec2d velocity( const vec2d &p ) {
	if( g_version == 0 ) {
		return vec2d();
	} else {
		if( p[0] < 0.5 ) return vec2d(0.4,0.04);
		else return vec2d(-0.4,-0.04);
	}
	return vec2d();
}
//
extern "C" const char *license() {
	return "MIT";
}