/*
**	river3.cpp
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
#include <shiokaze/polygon/polygon3_interface.h>
#include <shiokaze/polygon/polygon3_utility.h>
#include <shiokaze/graphics/graphics_engine.h>
#include <shiokaze/graphics/graphics_utility.h>
#include <shiokaze/meshlevelset/meshlevelset_interface.h>
#include <string>
#include <utility>
#include <cmath>
//
SHKZ_USING_NAMESPACE
//
static double g_flux_flow (0.15);
static polygon3_ptr g_polygon;
static std::vector<vec3d> g_rock_vertices, g_fish_vertices;
static std::vector<std::vector<size_t> > g_rock_faces, g_fish_faces;
static meshlevelset_ptr g_rock_levelset, g_fish_levelset;
static bool g_setup_done (false);
static std::string name("River Scene 3D"), argname("River");
//
static const unsigned g_fish_number (5);
static vec3d g_position_list[] = { vec3d(0.5,0.07,0.25), vec3d(0.5,0.07,0.6), vec3d(0.5,0.07,1.0), vec3d(0.48,0.07,1.25), vec3d(0.48,0.07,1.6) };
static double g_global_time_shift (0.4);
static double g_timeshift_list[] = { 0.4, 0.2, 0.0, 0.3, 0.5 };
//
extern "C" void load( configuration &config ) {
	//
	configuration::auto_group group(config,name,argname);
	g_polygon = polygon3_interface::quick_load_module(config,"polygon3");
	g_rock_levelset = meshlevelset_interface::quick_load_module(config,"SDFGen");
	g_fish_levelset = meshlevelset_interface::quick_load_module(config,"SDFGen");
}
//
extern "C" void unload () {
	g_polygon.reset();
	g_rock_levelset.reset();
	g_fish_levelset.reset();
}
//
extern "C" void configure( configuration &config ) {
	//
	g_polygon->recursive_configure(config);
	g_rock_levelset->recursive_configure(config);
	g_fish_levelset->recursive_configure(config);
	//
	configuration::auto_group group(config,name,argname);
	config.get_double("GlobalTimeShift",g_global_time_shift,"Global time shift");
}
//
static vec3d get_fish_transform( unsigned n, double time ) {
	//
	const double start_time (0.79);
	const double arc (0.5);
	const double jump (1.0);
	const vec3d &base_position = g_position_list[n];
	//
	time = time+g_global_time_shift-g_timeshift_list[n];
	vec3d position (base_position+vec3d(1.0,0.0,0.0));
	position += time * vec3d(-1.0,0.0,0.0);
	//
	auto parabolic_func = []( const double x ) {
		return x*(1.0-x);
	};
	if( time > start_time ) position += vec3d(0.0,jump*parabolic_func((time-start_time)/arc),0.0);
	if( position[0] < 0.5 ) position[1] = std::max(0.215,position[1]);
	//
	return position;
}
//
extern "C" void initialize( const shape3 &shape, double dx ) {
	if( ! g_setup_done ) {
		//
		// Rock part
		g_rock_levelset->recursive_initialize({{"dx",&dx}});
		g_polygon->load_mesh(filesystem::find_resource_path("objects","rocks.ply"));
		g_polygon->get_mesh(g_rock_vertices,g_rock_faces);
		polygon3_utility::transform(g_rock_vertices,vec3d(0.45,0.05,0.17),0.28);
		g_rock_levelset->set_mesh(g_rock_vertices,g_rock_faces);
		g_rock_levelset->generate_levelset();
		//
		// Fish part
		g_fish_levelset->recursive_initialize({{"dx",&dx}});
		g_polygon->load_mesh(filesystem::find_resource_path("objects","fish.ply"));
		g_polygon->get_mesh(g_fish_vertices,g_fish_faces);
		//
		const vec3d center = polygon3_utility::get_center_of_gravity(g_fish_vertices,g_fish_faces);
		for( auto &v : g_fish_vertices ) {
			v -= center;
			v *= 0.05;
			std::swap(v[0],v[2]);
		}
		//
		g_fish_levelset->set_mesh(g_fish_vertices,g_fish_faces);
		g_fish_levelset->generate_levelset();
	}
	g_setup_done = true;
}
//
extern "C" std::map<std::string,std::string> get_default_parameters() {
	std::map<std::string,std::string> dictionary;
	dictionary["CollisionDomainBoundary"] = "No";
	dictionary["VolumeCorrection"] = "No";
	dictionary["SpecialBoundaryCondition"] = "1";
	dictionary["VolumeCorrection"] = "No";
	dictionary["ResolutionX"] = "64";
	dictionary["ResolutionY"] = "32";
	dictionary["ResolutionZ"] = "128";
	dictionary["ViewScale"] = "2.0";
	dictionary["TargetPos"] = "0.5,0.05,0.8";
	dictionary["OriginPos"] = "1.7,0.75,-1.25";
	return dictionary;
}
//
extern "C" void set_boundary_flux( double time, Real flux[DIM3][2] ) {
	flux[0][0] = g_flux_flow;
}
//
extern "C" double fluid( const vec3d &p ) {
	double value (1.0);
	value = std::min(value,utility::box(p,vec3d(-1.0,-1.0,-1.0),vec3d(0.425,0.22,4.0)));
	value = std::min(value,p[1]-0.065);
	return value;
}
//
extern "C" double solid( const vec3d &p ) {
	double value (1.0);
	value = std::min(value,utility::box(p,vec3d(-1.0,-1.0,-1.0),vec3d(0.48,0.12,2.0)));
	value = std::min(value,std::min(g_rock_levelset->get_levelset(p),g_rock_levelset->get_levelset(p-vec3d(0.0,0.0,1.0))));
	return value;
}
//
extern "C" std::pair<double,vec3d> moving_solid( double time, const vec3d &p ) {
	//
	const double dt (0.01);
	double value (1.0);
	vec3d u;
	//
	for( unsigned n=0; n<g_fish_number; ++n ) {
		const vec3d pos = get_fish_transform(n,time);
		if((pos-p).len() < 0.075 ) {
			value = g_fish_levelset->get_levelset(p-pos);
			u = (get_fish_transform(n,time+dt)-get_fish_transform(n,time)) / dt;
		}
	}
	if( p[0] > 0.1 ) {
		return std::make_pair(value,u);
	} else {
		return std::make_pair(1.0,vec3d());
	}
}
//
extern "C" double solid_visualize( const vec3d &p ) {
	const double eps (0.01);
	return std::max(solid(p),utility::box(p,vec3d(eps,-1.0,eps),vec3d(1.0,1.0,2.0-eps)));
}
//
using polygon_list3 = std::vector<std::pair<std::vector<vec3d>,std::vector<std::vector<size_t> > > >;
extern "C" void export_moving_poygon( polygon_list3 &polygons ) {
	for( unsigned n=0; n<g_fish_number; ++n ) {
		polygons.push_back(std::make_pair(g_fish_vertices,g_fish_faces));
	}
}
//
extern "C" void get_moving_polygon_transforms( double time, std::vector<vec3d> &translations, std::vector<vec3d> &rotations ) {
	for( unsigned n=0; n<g_fish_number; ++n ) {
		translations.push_back(get_fish_transform(n,time));
		rotations.push_back(vec3d());
	}
}
//
extern "C" void draw( graphics_engine &g, double time ) {
	//
	auto draw_mesh = [&]( const vec3d &center, const std::vector<vec3d> &vertices, const std::vector<std::vector<size_t> > &faces ) {
		g.color4(1.0,1.0,1.0,0.5);
		for( size_t i=0; i<faces.size(); i++ ) {
			g.begin(graphics_engine::MODE::LINE_LOOP);
			for( unsigned j=0; j<faces[i].size(); j++ ) {
				g.vertex3v((vertices[faces[i][j]]+center).v);
			}
			g.end();
		}
	};
	//
	for( unsigned n=0; n<g_fish_number; ++n ) {
		const vec3d &p = get_fish_transform(n,time);
		draw_mesh(p,g_fish_vertices,g_fish_faces);
	}
}
//
extern "C" const char *license() {
	return "MIT";
}
//