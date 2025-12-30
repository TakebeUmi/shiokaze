/*
**	cloud_oc.h
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on April 10, 2017.
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
#ifndef SHKZ_cloud_oc_H
#define SHKZ_cloud_oc_H
//
//#include <shiokaze/array/array3.h>
//#include <shiokaze/array/macarray3.h>
#include <shiokaze/parallel/parallel_driver.h>
#include <shiokaze/ui/drawable.h>
//#include <shiokaze/advection/macadvection3_interface.h>
#include <shiokaze/utility/gridutility3_interface.h>
#include <shiokaze/utility/macutility3_interface.h>
#include <shiokaze/utility/macstats3_interface.h>
#include <shiokaze/visualizer/gridvisualizer3_interface.h>
#include <shiokaze/visualizer/macvisualizer3_interface.h>
#include <shiokaze/projection/macproject3_interface.h>
#include <shiokaze/timestepper/timestepper_interface.h>
#include <shiokaze/utility/graphplotter_interface.h>
#include <shiokaze/core/dylibloader.h>
//
#include <shiokaze/visualizer/gridvisualizer3_cloud_interface.h>
//

#include "macoctreegrid3.h"
#include "macoctreeproject3.h"
//#include "macoctreehelper3.h"
#include "macoctreesegregator3.h"
#include "macoctreesizingfunc3.h"
#include <shiokaze/meshexporter/meshexporter3_interface.h>
#include <shiokaze/cellmesher/cellmesher3_interface.h>
#include <shiokaze/graphics/graphics_interface.h>
#include "macoctreemesher3.h"
//
SHKZ_BEGIN_NAMESPACE
//
using namespace macotreeliquid3_namespace;

class cloud_oc : public drawable {
public:
	//
	cloud_oc();
	LONG_NAME("MAC Smoke 3D")
	ARGUMENT_NAME("Smoke")
	//
protected:
	//
	void write_to_txt(std::string type) const;
	virtual void setup_window( std::string &name, int &width, int &height ) const override;
	//virtual void drag( double x, double y, double z, double u, double v, double w ) override;
	virtual void idle() override;
	virtual void draw( graphics_engine &g ) const override;
	virtual bool should_quit() const override; /*{ return m_timestepper->should_quit(); }*/
	virtual bool should_screenshot() const override { return m_timestepper->should_export_frame(); }
	virtual void load( configuration &config ) override;
	virtual void configure( configuration &config ) override;
	virtual void post_initialize( bool initialized_from_file ) override;
	//
	// macarray3<Real> m_velocity{this};
	// macarray3<Real> m_solid_velocity{this};
	// macarray3<Real> m_external_force{this};
	// array3<Real> m_density{this};
	// array3<Real> m_accumulation{this};
	//
	// array3<Real> m_fluid{this};
	// array3<Real> m_solid{this};
	//
	std::vector<vec3d> m_dust_particles;
	//
	shape3 m_shape;
	double m_dx;
	//bool m_force_exist;
	unsigned m_graph_id;
	shape3 m_solid_shape;
	double m_solid_dx;
	bool m_do_inject;
	double m_injected_volume;

	//gridvisualizer3_driver m_solid_gridvisualizer{this,"gridvisualizer3"};
	//array3<Real> m_solid_visualize{this};

	//
	//追加
    grid3 m_grid_0{this};
    grid3 m_grid_1{this};
    grid3 *m_grid {&m_grid_0};
    grid3 *m_grid_prev {&m_grid_1};
    double m_accumulated_CFL {0.0};
	UnionFind uf;
    macoctreeproject3 m_macoctreeproject{this};

	struct Parameters {
		bool use_FLIP {false};
		bool mouse_interaction {false};
		bool use_dust {false};
		double minimal_density {0.01};
		unsigned r_sample {4};
		bool show_graph {false};
		unsigned extrapolated_width {3};
		double buoyancy_factor {2.0};
		bool render_density {false};
		unsigned render_sample_count {128};
		double volume_scale {40.0};
		//
        unsigned min_resolution {64};//ここをいじると最小解像度が変化する
        bool use_sizing_func {true};
        unsigned initial_refinement {3};
        double maximal_CFL_accumulation {1.0};

	double surftens_k {0.0};
	unsigned erode_width {0};
	bool render_mesh {false};
	bool render_wireframe {false};
	bool render_grid {false};
	bool remove_quater {false};
	double z {0.5};
	bool export_svg {true};
	//unsigned render_sample_count {8};
	unsigned save_interval {100};
	double scale {12000.0};
	vec3d target {0.5*scale,0.5*scale,0.5*scale};
	vec3d origin {1.3*scale,0.3*scale,1.7*scale};
	vec3d gravity {0.0,-9.8,0.0};
	double PICFLIP {0.98};
	bool render_transparent {false};
	bool volume_correction {true};
	bool regional_volume_correction {false};
	bool maccormack {false};
	unsigned render_transparent_sample_count {32};
	bool transfer_file {true};
	int debug_mode {0};

	double radius {0.2};
	vec2d center {0.5,0.5};
	double source_coeff {1.0};
	double Time_coeff {30.0};
	};
	//
	Parameters m_param;
	//
	environment_setter arg_shape{this,"shape",&m_shape};
	environment_setter arg_dx{this,"dx",&m_dx};
	//
	macproject3_driver m_macproject{this,"macpressuresolver3"};
	//macadvection3_driver m_macadvection{this,"macadvection3"};
	gridvisualizer3_driver m_gridvisualizer{this,"gridvisualizer3"};
	gridvisualizer3_cloud_driver m_gridvisualizer_oc{this,"gridvisualizer3_cloud"};
	gridutility3_driver m_gridutility{this,"gridutility3"};
	graphplotter_driver m_graphplotter{this,"graphplotter"};
	macstats3_driver m_macstats{this,"macstats3"};
	macvisualizer3_driver m_macvisualizer{this,"macvisualizer3"};
	timestepper_driver m_timestepper{this,"timestepper"};
	//macutility3_driver m_macutility{this,"macutility3"};
	//
	parallel_driver m_parallel{this};
	dylibloader m_dylib;
	//added
	bool m_should_quit_on_save {false}; 
	//
	//
	//added
		std::function<double(const vec3d &)> m_solid_func, m_combined_solid_func;
		std::function<void(double,Real [DIM3][2])> m_set_boundary_flux;
		std::function<vec3d(double)> m_gravity_func;
		std::function<void(graphics_engine &,double)> m_draw_func;
		std::function<std::pair<double,vec3d>( double, const vec3d &)> m_moving_solid_func;
		//
		using polygon_list3 = std::vector<std::pair<std::vector<vec3d>,std::vector<std::vector<size_t> > > >;
		std::function<void( polygon_list3 &polygons )> m_export_moving_poygon_func;
		std::function<void( double time, std::vector<vec3d> &, std::vector<vec3d> &)> m_get_moving_polygon_transforms_func;
		//
		std::string m_export_path;
		//

		macoctreesegregator3 m_macoctreesegregator{this};
		macoctreemesher3 m_macoctreemesher{this};
		meshutility3_driver m_meshutility{this,"meshutility3"};
		macoctreesizingfunc3 m_macoctreesizingfunc{this};
		meshexporter3_driver m_mesh_exporter{this,"meshexporter3"};
		cellmesher3_driver m_solid_mesher{this,"marchingcubes"};
		graphics_interface_driver m_svg_writer{this,"graphics_svg"};
		//
	//added
	//
	//virtual void inject_external_force( macarray3<Real> &velocity );
	//virtual void add_buoyancy_force( macarray3<Real> &velocity, const array3<Real> &density, double dt );
	//virtual void advect_dust_particles( const macarray3<Real> &velocity, double dt );
	//virtual void add_source ( macarray3<Real> &velocity, array3<Real> &density, double time, double dt );
	void add_source_cloud ( double time, double dt, grid3 &grid );
	// virtual void rasterize_dust_particles( array3<Real> &rasterized_density );
	// virtual void draw_dust_particles( graphics_engine &g ) const;
	virtual void export_density () const;
	virtual void do_export_density( int frame ) const;
	//virtual void add_to_graph();
	virtual void render_density( int frame ) const;
	//virtual void add_source_oc ( double time, double dt );
	//
	//added
		void export_moving_polygon();
		void do_export_solid_mesh( const array3<Real> &solid );
		void do_export_empty_solid_mesh( bool force=false );
		void export_mesh( int frame );
		void render_mesh( unsigned frame ) const;
		void save_state();

		//virtual void do_inject_external_density( double dt, double time, unsigned step );
	void microphysics_cloud(grid3 &grid, double dt);
	void print_position(grid3 &grid) const;
	double perlin_noise_source(const vec3d &p) const;
	double circle_source(const vec3d &p) const;
	//added
};
//
SHKZ_END_NAMESPACE
//
#endif