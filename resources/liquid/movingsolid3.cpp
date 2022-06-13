/*
**	movingsolid3.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on January 5, 2020.
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
#define _USE_MATH_DEFINES
#include <shiokaze/math/vec.h>
#include <shiokaze/core/configuration.h>
#include <shiokaze/graphics/graphics_engine.h>
#include <shiokaze/graphics/graphics_utility.h>
#include <string>
#include <cmath>
#include <tuple>
#include <utility>
#include <vector>
//
SHKZ_USING_NAMESPACE
//
static double g_water_level (0.245);
static double g_speed (4.0);
static double g_amplitude (0.25);
static vec3d g_center (0.5,0.0,0.5);
static double g_r (0.05);
static double g_height (0.45);
static double g_flux_flow (0.0);
static bool g_debug_first_frame {false};
static unsigned g_subdiv_num {20};
const static double g_bottom_margin {0.01};
//
extern "C" void configure( configuration &config ) {
	//
	configuration::auto_group group(config,"Moving Solid Scene 3D","MovingSolid");
	config.get_double("WaterLevel",g_water_level,"Water level");
	config.get_double("Speed",g_speed,"Speed of obstacle");
	config.get_double("Amplitude",g_amplitude,"Speed of obstacle");
	config.get_double("Radius",g_r,"Radius");
	config.get_double("Height",g_height,"Height");
	config.get_vec3d("Center",g_center.v,"Center position");
	config.get_double("Flux",g_flux_flow,"Flux velocity on walls");
	config.get_unsigned("SubdivNum",g_subdiv_num,"Subdivision number for cylinder");
	config.get_bool("DebugFirstFrame",g_debug_first_frame,"Debug first frame");
}
//
static double cylinder_levelset( const vec3d &p, const vec3d &center, double height, double r ) {
	return std::max(std::hypot(center[0]-p[0],center[2]-p[2])-r,p[1]-height);
}
//
static void cylinder_mesh( const vec3d &center, double height, double r,
						   std::vector<vec3d> &vertices, std::vector<std::vector<size_t> > &faces ) {
	//
	vertices.clear();
	faces.clear();
	//
	const double theta = 2.0 * M_PI / g_subdiv_num;
	for( unsigned n=0; n<g_subdiv_num; ++n ) {
		//
		const double t = theta * n;
		vertices.push_back(center+vec3d(r*cos(t),0.0,r*sin(t)));
		vertices.push_back(center+vec3d(r*cos(t),height,r*sin(t)));
	}
	//
	for( unsigned n=0; n<g_subdiv_num; ++n ) {
		//
		std::vector<size_t> face;
		face.push_back(2*n+1);
		face.push_back(2*n);
		const unsigned m = (n+1) % g_subdiv_num;
		face.push_back(2*m);
		face.push_back(2*m+1);
		faces.push_back(face);
	}
	//
	const unsigned idx0 = vertices.size();
	vertices.push_back(center);
	const unsigned idx1 = vertices.size();
	vertices.push_back(center+vec3d(0.0,height,0.0));
	//
	for( unsigned n=0; n<g_subdiv_num; ++n ) {
		//
		const unsigned m = (n+1) % g_subdiv_num;
		std::vector<size_t> face0;
		face0.push_back(2*m);
		face0.push_back(2*n);
		face0.push_back(idx0);
		faces.push_back(face0);
		//
		std::vector<size_t> face1;
		face1.push_back(2*n+1);
		face1.push_back(2*m+1);
		face1.push_back(idx1);
		faces.push_back(face1);
	}
}
//
extern "C" std::pair<double,vec3d> moving_solid( double time, const vec3d &p ) {
	//
	const vec3d center = g_center+vec3d(-g_amplitude*cos(g_speed*time),-g_bottom_margin,0.0);
	const double d = cylinder_levelset(p,center,g_height,g_r);
	//
	vec3d u;
	if( d <= 0.0 ) {
		u = vec3d(g_amplitude*g_speed*sin(g_speed*time),0.0,0.0);
	}
	return std::make_pair(d,u);
}
//
extern "C" vec3d velocity( const vec3d &p ) {
	return vec3d(g_flux_flow,0.0,0.0);
}
//
extern "C" void set_boundary_flux( double time, Real flux[DIM3][2] ) {
	//
	flux[0][0] = g_flux_flow;
	flux[0][1] = g_flux_flow;
}
//
extern "C" double fluid( const vec3d &p ) {
	return p[1]-g_water_level;
}
//
extern "C" double solid( const vec3d &p ) {
	if( g_debug_first_frame ) {
		return moving_solid(0.0,p).first;
	} else {
		return 1.0;
	}
}
//
using polygon_list3 = std::vector<std::pair<std::vector<vec3d>,std::vector<std::vector<size_t> > > >;
extern "C" void export_moving_poygon( polygon_list3 &polygons ) {
	//
	std::vector<vec3d> vertices;
	std::vector<std::vector<size_t> > faces;
	cylinder_mesh(vec3d(),g_height,g_r,vertices,faces);
	polygons.push_back(std::make_pair(vertices,faces));
}
//
static vec3d get_center( double time ) {
	return g_center+vec3d(-g_amplitude*cos(g_speed*time),-g_bottom_margin,0.0);
}
//
extern "C" void get_moving_polygon_transforms( double time, std::vector<vec3d> &translations, std::vector<vec3d> &rotations ) {
	//
	translations.push_back(get_center(time));
	rotations.push_back(vec3d());
}
//
extern "C" void draw( graphics_engine &g, double time ) {
	//
	const vec3d center = get_center(time);
	std::vector<vec3d> vertices;
	std::vector<std::vector<size_t> > faces;
	cylinder_mesh(center,g_height,g_r,vertices,faces);
	//
	for( size_t i=0; i<faces.size(); i++ ) {
		g.begin(graphics_engine::MODE::LINE_LOOP);
		for( unsigned j=0; j<faces[i].size(); j++ ) g.vertex3v(vertices[faces[i][j]].v);
		g.end();
	}
}
//
extern "C" const char *license() {
	return "MIT";
}
//