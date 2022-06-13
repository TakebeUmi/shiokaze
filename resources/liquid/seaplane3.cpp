/*
**	seaplane3.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on January 8, 2020.
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
#include <shiokaze/utility/utility.h>
#include <shiokaze/core/configuration.h>
#include <shiokaze/core/filesystem.h>
#include <shiokaze/polygon/polygon3_interface.h>
#include <shiokaze/polygon/polygon3_utility.h>
#include <shiokaze/meshlevelset/meshlevelset_interface.h>
#include <shiokaze/graphics/graphics_engine.h>
#include <shiokaze/graphics/graphics_utility.h>
#include <string>
#include <cmath>
//
SHKZ_USING_NAMESPACE
//
static bool g_debug {false};
static double g_time_shift {3.7};
static double g_level {13.0};
static double g_jitter_eps {2.0};
static unsigned g_subdiv_num (5);
static double g_leg_r (0.1);
static vec3d g_center(45.0,g_level+2.0,30.0);
static double g_flux_flow = 30.0;
static double g_flux_flow_start = 5.0;
static bool g_setup_done (false);
static double g_start (1000.0);
static double g_take_off_time(11.3);
static bool g_lock_rotation {false};
static polygon3_ptr g_polygon;
//
struct object3 {
	std::vector<vec3d> vertices;
	std::vector<std::vector<size_t> > faces;
	meshlevelset_ptr levelset;
};
//
object3 g_body;
object3 g_wings;
//
struct cylinder_obj3 {
	cylinder_obj3( const vec3d &p0, const vec3d &p1, double r ) {
		this->p0 = p0;
		this->p1 = p1;
		this->r = r;
	}
	vec3d p0, p1;
	double r;
};
//
static std::vector<cylinder_obj3> g_cylinders;
static std::string name ("Seaplane Scene 3D","Seaplane"), argname ("Seaplane3");
//
extern "C" void load( configuration &config ) {
	//
	configuration::auto_group group(config,name,argname);
	g_polygon = polygon3_interface::quick_load_module(config,"polygon3");
	g_body.levelset = meshlevelset_interface::quick_load_module(config,"SDFGen");
	g_wings.levelset = meshlevelset_interface::quick_load_module(config,"SDFGen");
}
//
extern "C" void unload () {
	g_polygon.reset();
	g_body.levelset.reset();
	g_wings.levelset.reset();
}
//
extern "C" void configure( configuration &config ) {
	//
	g_polygon->recursive_configure(config);
	g_body.levelset->recursive_configure(config);
	g_wings.levelset->recursive_configure(config);
	//
	double level (g_level);
	//
	configuration::auto_group group(config,name,argname);
	config.get_double("Level",level,"Water level");
	config.get_double("JitterEps",g_jitter_eps,"Jitetr eps for water level");
	config.get_double("TimeShift",g_time_shift,"Time shift");
	config.get_double("FluxFlow",g_flux_flow,"Flux flow");
	config.get_double("StopStart",g_start,"Stop start time");
	config.get_double("TakeOffTime",g_take_off_time,"Take off time");
	config.get_double("FluxFlowStart",g_flux_flow_start,"Flux flow at start");
	config.get_bool("LockRotation",g_lock_rotation,"Lock rotation");
	config.get_bool("DebugMode",g_debug,"Debug mode");
	g_center[1] += g_jitter_eps;
	//
	const double diff = level-g_level;
	g_level += diff;
	g_center[1] += diff;
}
//
// **************************************************************
// Specification from: https://en.wikipedia.org/wiki/Maule_M-5
// **************************************************************
// Length: 23 ft 6 in (7.16 m)
// Wingspan: 30 ft 10 in (9.40 m)
// Landing length: (approx) 200m
// **************************************************************
//
extern "C" std::map<std::string,std::string> get_default_parameters() {
	//
	std::map<std::string,std::string> dictionary;
	dictionary["ResolutionX"] = "128";
	dictionary["ResolutionY"] = "64";
	dictionary["ResolutionZ"] = "128";
	dictionary["ViewScale"] = "60.0";
	dictionary["MaxFrame"] = "800";
	dictionary["CollisionDomainBoundary"] = "No";
	dictionary["FPS"] = "60.0";
	dictionary["TargetPos"] = "47,11,30";
	dictionary["OriginPos"] = "-15,22,30";
	dictionary["ZPosition"] = "28.4";
	dictionary["MeshPadding"] = "32";
	dictionary["VelocityScale"] = "0.008";
	return dictionary;
}
//
static double cylinder_levelset( const vec3d &p, const cylinder_obj3 &obj ) {
	//
	const vec3d nvec = (obj.p1-obj.p0).normal();
	const double h0 = obj.p0 * nvec;
	const double h1 = obj.p1 * nvec;
	const double d0 = p * nvec;
	const double d1 = ((obj.p0-p) ^ nvec).len();
	return std::max(h0-d0,std::max(d0-h1,d1-obj.r));
}
//
static void cylinder_mesh( const cylinder_obj3 &obj,
						   std::vector<vec3d> &vertices, std::vector<std::vector<size_t> > &faces ) {
	//
	vertices.clear();
	faces.clear();
	//
	const vec3d nvec = (obj.p1-obj.p0).normal();
	const vec3d xvec = nvec ^ vec3d(1.0,0.0,0.0);
	const vec3d yvec = nvec ^ xvec;
	//
	const double theta = 2.0 * M_PI / g_subdiv_num;
	for( unsigned n=0; n<g_subdiv_num; ++n ) {
		//
		const double t = theta * n;
		vertices.push_back(obj.p0+obj.r*(xvec*cos(t)+yvec*sin(t)));
		vertices.push_back(obj.p1+obj.r*(xvec*cos(t)+yvec*sin(t)));
	}
	//
	for( unsigned n=0; n<g_subdiv_num; ++n ) {
		//
		std::vector<size_t> face;
		face.push_back(2*n);
		face.push_back(2*n+1);
		const unsigned m = (n+1) % g_subdiv_num;
		face.push_back(2*m+1);
		face.push_back(2*m);
		faces.push_back(face);
	}
	//
	const unsigned idx0 = vertices.size();
	vertices.push_back(obj.p0);
	const unsigned idx1 = vertices.size();
	vertices.push_back(obj.p1);
	//
	for( unsigned n=0; n<g_subdiv_num; ++n ) {
		//
		const unsigned m = (n+1) % g_subdiv_num;
		std::vector<size_t> face0;
		face0.push_back(2*n);
		face0.push_back(2*m);
		face0.push_back(idx0);
		faces.push_back(face0);
		//
		std::vector<size_t> face1;
		face1.push_back(2*m+1);
		face1.push_back(2*n+1);
		face1.push_back(idx1);
		faces.push_back(face1);
	}
}
//
static void setup_seaplane( double scale ) {
	//
	double local_dx = 9.4 / (g_debug ? 8.0 : 128.0);
	g_body.levelset->recursive_initialize({{"dx",&local_dx}});
	g_wings.levelset->recursive_initialize({{"dx",&local_dx}});
	//
	g_polygon->load_mesh(filesystem::find_resource_path("objects","seaplane_body.ply"));
	g_polygon->get_mesh(g_body.vertices,g_body.faces);
	//
	const vec3d center_of_gravity(1.45,2.85,0.0);
	const double real_scale = scale * 0.05;
	for( auto &v : g_body.vertices ) v = (v-center_of_gravity)*real_scale;
	//
	g_polygon->load_mesh(filesystem::find_resource_path("objects","seaplane_wings.ply"));
	g_polygon->get_mesh(g_wings.vertices,g_wings.faces);
	for( auto &v : g_wings.vertices ) v = (v-center_of_gravity)*real_scale;
	//
	g_body.levelset->set_mesh(g_body.vertices,g_body.faces);
	g_body.levelset->generate_levelset();
	//
	g_wings.levelset->set_mesh(g_wings.vertices,g_wings.faces);
	g_wings.levelset->generate_levelset();
	//
	// Measure the scale
	std::vector<vec3d> merged_vertices;
	merged_vertices.insert(merged_vertices.end(),g_body.vertices.begin(),g_body.vertices.end());
	merged_vertices.insert(merged_vertices.end(),g_wings.vertices.begin(),g_wings.vertices.end());
	vec3d corner0, corner1;
	polygon3_utility::compute_AABB(merged_vertices,corner0,corner1);
	//
	printf( "(Length: %.2fm, Wingspan: %.2fm, Height: %.2fm)...", (corner1-corner0)[0], (corner1-corner0)[2], (corner1-corner0)[1] );
	//
	g_cylinders.push_back(cylinder_obj3(real_scale*vec3d(0.0,-2.4,-1.6),real_scale*vec3d(1.0,-0.8,-0.5),real_scale*g_leg_r));
	g_cylinders.push_back(cylinder_obj3(real_scale*vec3d(3.0,-2.4,-1.6),real_scale*vec3d(2.0,-0.8,-0.5),real_scale*g_leg_r));
	g_cylinders.push_back(cylinder_obj3(real_scale*vec3d(3.0,-2.4,-1.6),real_scale*vec3d(1.0,-0.8,-0.5),real_scale*g_leg_r));
	g_cylinders.push_back(cylinder_obj3(real_scale*vec3d(0.0,-2.4,1.6),real_scale*vec3d(1.0,-0.8,0.5),real_scale*g_leg_r));
	g_cylinders.push_back(cylinder_obj3(real_scale*vec3d(3.0,-2.4,1.6),real_scale*vec3d(2.0,-0.8,0.5),real_scale*g_leg_r));
	g_cylinders.push_back(cylinder_obj3(real_scale*vec3d(3.0,-2.4,1.6),real_scale*vec3d(1.0,-0.8,0.5),real_scale*g_leg_r));
	g_cylinders.push_back(cylinder_obj3(real_scale*vec3d(2.0,0.0,-0.7),real_scale*vec3d(2.0,1.3,-2.4),real_scale*g_leg_r));
	g_cylinders.push_back(cylinder_obj3(real_scale*vec3d(2.0,0.0,0.7),real_scale*vec3d(2.0,1.3,2.4),real_scale*g_leg_r));
}
//
extern "C" void initialize( const shape3 &shape, double dx ) {
	if( ! g_setup_done ) {
		setup_seaplane(14.0);
		g_setup_done = true;
	}
}
//
extern "C" vec3d velocity( const vec3d &p ) {
	return vec3d(-g_flux_flow,0.0,0.0);
}
//
static double get_flux( double time ) {
	return std::min(0.0,(g_flux_flow/g_start)*std::max(0.0,time-g_start)-g_flux_flow);
}
//
extern "C" void set_boundary_flux( double time, Real flux[DIM3][2] ) {
	//
	const double f = get_flux(time);
	flux[0][0] = f;
	flux[0][1] = f;
}
//
extern "C" vec3d gravity( double time ) {
	//
	static const double dt = 0.01;
	const double a = (get_flux(time+dt)-get_flux(time)) / dt;
	return vec3d(a,-9.8,0.0);
}
//
static vec3d convert_position( const vec3d &p, const vec3d &center, double theta, bool forward ) {
	vec3d p1;
	if( forward ) {
		p1[0] = p[0]*cos(theta)-p[1]*sin(theta);
		p1[1] = p[0]*sin(theta)+p[1]*cos(theta);
		p1[2] = p[2];
		return center+p1;
	} else {
		vec3d p0 = p-center;
		p1[0] = p0[0]*cos(theta)+p0[1]*sin(theta);
		p1[1] = -p0[0]*sin(theta)+p0[1]*cos(theta);
		p1[2] = p0[2];
		return p1;
	}
}
//
static double get_rotation( double time ) {
	if( g_lock_rotation ) time = std::min(2.18,time);
	time += g_time_shift;
	return 0.13+exp(-std::max(0.0,time-7.0))*0.06*exp(sin(time))*sin(1.5*time);
}
//
static vec3d get_center( double time ) {
	time += g_time_shift;
	const double s (0.5);
	return g_center+vec3d(0.0,s*std::max(0.0,5.0-time),0.0) + vec3d(0.0,std::max(0.0,time-g_take_off_time),0.0);
}
//
extern "C" std::pair<double,vec3d> moving_solid( double time, const vec3d &p ) {
	//
	double value (1.0);
	vec3d symm_p (p);
	// Force symmetrize
	if( symm_p[2] > g_center[2] ) symm_p[2] = g_center[2]-(symm_p[2]-g_center[2]);
	vec3d p1 = convert_position(symm_p,get_center(time),get_rotation(time),false);
	value = std::min(value,g_body.levelset->get_levelset(p1));
	value = std::min(value,g_wings.levelset->get_levelset(p1));
	//if( ! time ) value = std::min(value,(p-g_center).len()-2.0);
	return std::make_pair(value,vec3d());
}
//
extern "C" double fluid( const vec3d &p ) {
	return p[1]-g_level-g_jitter_eps;
}
//
using polygon_list3 = std::vector<std::pair<std::vector<vec3d>,std::vector<std::vector<size_t> > > >;
extern "C" void export_moving_poygon( polygon_list3 &polygons ) {
	//
	std::vector<vec3d> vertices;
	std::vector<std::vector<size_t> > faces;
	//
	auto add_mesh = [&]( const std::vector<vec3d> &_vertices, const std::vector<std::vector<size_t> > &_faces ) {
		size_t head = vertices.size();
		std::vector<std::vector<size_t> > new_faces(_faces);
		for( auto &f : new_faces ) for( auto &e : f ) e += head;
		vertices.insert(vertices.end(),_vertices.begin(),_vertices.end());
		faces.insert(faces.end(),new_faces.begin(),new_faces.end());
	};
	//
	add_mesh(g_body.vertices,g_body.faces);
	add_mesh(g_wings.vertices,g_wings.faces);
	//
	for( const auto &cylinder : g_cylinders ) {
		std::vector<vec3d> tmp_vertices;
		std::vector<std::vector<size_t> > tmp_faces;
		cylinder_mesh(cylinder,tmp_vertices,tmp_faces);
		add_mesh(tmp_vertices,tmp_faces);
	}
	//
	polygons.push_back(std::make_pair(vertices,faces));
}
//
extern "C" void get_moving_polygon_transforms( double time, std::vector<vec3d> &translations, std::vector<vec3d> &rotations ) {
	//
	translations.push_back(get_center(time));
	rotations.push_back(vec3d(0.0,0.0,get_rotation(time)/M_PI*180.0));
}
//
extern "C" void draw( graphics_engine &g, double time ) {
	//
	auto draw_mesh = [&]( const vec3d &center, double theta, const std::vector<vec3d> &vertices, const std::vector<std::vector<size_t> > &faces ) {
		g.color4(1.0,1.0,1.0,0.5);
		for( size_t i=0; i<faces.size(); i++ ) {
			g.begin(graphics_engine::MODE::LINE_LOOP);
			for( unsigned j=0; j<faces[i].size(); j++ ) {
				g.vertex3v(convert_position(vertices[faces[i][j]],center,theta,true).v);
			}
			g.end();
		}
	};
	//
	// XY rotation (z-axis)
	const double z_theta = get_rotation(time);
	const vec3d center = get_center(time);
	//
	draw_mesh(center,z_theta,g_body.vertices,g_body.faces);
	draw_mesh(center,z_theta,g_wings.vertices,g_wings.faces);
	//
	for( const auto &cylinder : g_cylinders ) {
		//
		std::vector<vec3d> vertices;
		std::vector<std::vector<size_t> > faces;
		cylinder_mesh(cylinder,vertices,faces);
		draw_mesh(center,z_theta,vertices,faces);
	}
}
//
extern "C" const char *license() {
	return "MIT";
}