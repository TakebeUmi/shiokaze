/*
**	rain3.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on November 26, 2019.
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
#include <string>
#include <cmath>
#include <random>
//
SHKZ_USING_NAMESPACE
//
static double g_water_radius (0.022);
static double g_water_level (0.18);
static double g_inject_height (0.4);
static std::mt19937 g_rand_src(3);
static double g_inject_speed (7.5);
static bool g_fix_volume (true);
static vec3d g_inject_center;
//
extern "C" std::map<std::string,std::string> get_default_parameters() {
	std::map<std::string,std::string> dictionary;
	dictionary["ResolutionX"] = std::to_string(384);
	dictionary["ResolutionY"] = std::to_string(192);
	dictionary["ResolutionZ"] = std::to_string(384);
	dictionary["FPS"] = std::to_string(450);
	dictionary["MaxFrame"] = std::to_string(450);
	dictionary["CFL"] = std::to_string(2);
	dictionary["MacCormack"] = "No";
	return dictionary;
}
//
extern "C" void configure( configuration &config ) {
	configuration::auto_group group(config,"Rain Scene 3D","Rain");
	config.get_double("Radius",g_water_radius,"Radius of water");
	config.get_double("WaterLevel",g_water_level,"Water level");
	config.get_double("InjectHeight",g_inject_height,"Injection height");
	config.get_double("InjectSpeed",g_inject_speed,"Injection speed");
	config.get_bool("FixVolume",g_fix_volume,"Fix total volume");
}
//
extern "C" double fluid( const vec3d &p ) {
	return std::min(p[1]-g_water_level,(p-vec3d(0.5,g_inject_height,0.5)).len()-g_water_radius);
}
//
extern "C" vec3d velocity( const vec3d &p ) {
	if( p[1] > 0.25 ) return vec3d(0.0,-g_inject_speed,0.0);
	else return vec3d();
}
//
extern "C" bool check_inject( double dx, double dt, double time, unsigned step ) {
	//
	static unsigned next_frame (1);
	if( 0.015 * next_frame < time ) {
		std::uniform_real_distribution<double> rand_dist(0.1, 0.9);
		g_inject_center = vec3d(rand_dist(g_rand_src),g_inject_height,rand_dist(g_rand_src));
		next_frame ++;
		return true;
	}
	return false;
}
//
extern "C" bool inject( const vec3d &p, double dx, double dt, double time, unsigned step, double &fluid, vec3d &velocity ) {
	//
	fluid = (p-g_inject_center).len()-g_water_radius;
	velocity = vec3d(0.0,-g_inject_speed,0.0);
	return true;
}
//
extern "C" void post_inject( double dx, double dt, double time, unsigned step, double &volume_change ) {
	if( g_fix_volume ) volume_change = 0.0;
}
//
extern "C" const char *license() {
	return "MIT";
}
//