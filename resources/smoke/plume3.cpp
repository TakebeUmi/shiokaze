/*
**	plume3.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on August 15, 2017.
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
#include <map>
#include <string>
#include <cmath>
//
SHKZ_USING_NAMESPACE
// static double g_container_thickness (0.03);
// static double g_container_radius (0.5);
// static double g_container_height (0.3);
// static vec3d g_center (0.5,0.37,0.5);
// static double g_radius (0.075);
// static double g_level (0.245);
// static bool g_container (false);

// extern "C" void configure( configuration &config ) {
// 	configuration::auto_group group(config,"Waterdrop Scene 3D","Waterdrop");
// 	config.get_bool("Container",g_container,"Whether to place g_container");
// 	if( g_container ) {
// 		config.get_double("ContainerRadius",g_container_radius,"Radius of the solid hemisphere container");
// 		config.get_double("ContainerThickness",g_container_thickness,"Thickness of the solid hemisphere container");
// 		config.get_double("ContainerHeight",g_container_height,"Height of the solid hemisphere container");
// 	}
// 	config.get_double("Radius",g_radius,"Radius of water");
// 	config.get_vec3d("Center",g_center.v,"g_center of spherical water");
// 	config.get_double("Level",g_level,"Level of static water pool");
// }
//
static vec3d g_center (0.15,0.15,0.5);
static double g_radius (0.075);
static double g_level (0.245);

extern "C" std::map<std::string,std::string> get_default_parameters () {
	std::map<std::string,std::string> dictionary;
	return dictionary;
}
extern "C" double fluid( const vec3d &p) {
		vec3d center (0.15,0.15,0.5);
		double d = 0.0;
		if( (p-center).len() < 0.1 ) {
			// double dist = (p-center).len();
			// double r (0.075);
			// double s (10.0);
			// double v = std::min(10.0,std::max(0.0,s*(r-dist)/r));
			// u = 2.0 * vec3d(dt*v,0.0,0.0);
			d = 0.04; //* dt;
		}
		//d = 4.0;
	return d;
}
//
extern "C" void add ( const vec3d &p, vec3d &u, double &d, double time, double dt ) {
	if( time < 1.0 ) {
		vec3d center (0.15,0.15,0.5);
		if( (p-center).len() < 0.1 ) {
			double dist = (p-center).len();
			double r (0.075);
			double s (10.0);
			double v = std::min(10.0,std::max(0.0,s*(r-dist)/r));
			u = 2.0 * vec3d(1.5,0.0,0.0);
			d = 4.0 * dt;
		}
	}
}
//
extern "C" const char *license() {
	return "MIT";
}
//