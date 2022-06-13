/*
**	river2.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on January 16, 2020.
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
#include <shiokaze/utility/utility.h>
#include <string>
#include <cmath>
//
SHKZ_USING_NAMESPACE
//
static double g_flux_flow (0.15);
//
extern "C" void configure( configuration &config ) {
	configuration::auto_group group(config,"River Scene 2D","River");
}
//
extern "C" std::map<std::string,std::string> get_default_parameters() {
	std::map<std::string,std::string> dictionary;
	dictionary["CollisionDomainBoundary"] = "No";
	dictionary["VolumeCorrection"] = "No";
	dictionary["SpecialBoundaryCondition"] = "1";
	dictionary["VolumeCorrection"] = "No";
	return dictionary;
}
//
extern "C" void set_boundary_flux( double time, Real flux[DIM2][2] ) {
	flux[0][0] = g_flux_flow;
}
//
extern "C" double fluid( const vec2d &p ) {
	double value (1.0);
	value = std::min(value,utility::box(p,vec2d(-1.0,-1.0),vec2d(0.425,0.25)));
	value = std::min(value,p[1]-1.5/16.0);
	return value;
}
//
extern "C" double solid( const vec2d &p ) {
	double value (1.0);
	value = std::min(value,utility::box(p,vec2d(-1.0,-1.0),vec2d(0.485,0.12)));
	value = std::min(value,utility::box(p,vec2d(0.4,-1.0),vec2d(0.485,0.18)));
	return value;
}
//
extern "C" const char *license() {
	return "MIT";
}
//