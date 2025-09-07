/*
**	macoctreeliquid3.h
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on Februrary 8, 2019. All rights reserved.
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
#ifndef SHKZ_ADAPTIVELIQUID3_H
#define SHKZ_ADAPTIVELIQUID3_H
//
#include <shiokaze/ui/drawable.h>
#include <shiokaze/timestepper/timestepper_interface.h>
#include <shiokaze/visualizer/gridvisualizer3_interface.h>
#include <shiokaze/cellmesher/cellmesher3_interface.h>
#include <shiokaze/utility/macutility3_interface.h>
#include <shiokaze/utility/gridutility3_interface.h>
#include <shiokaze/utility/meshutility3_interface.h>
#include <shiokaze/graphics/graphics_interface.h>
#include <shiokaze/flip/macflip3_interface.h>
#include <shiokaze/meshexporter/meshexporter3_interface.h>
#include <shiokaze/array/array3.h>
#include <shiokaze/array/macarray3.h>
#include <shiokaze/core/dylibloader.h>
#include <vector>
#include <memory>
#include <string>
#include <utility>
//
#include "macoctreegrid3.h"
#include "macoctreemesher3.h"
#include "macoctreeproject3.h"
#include "macoctreesegregator3.h"
#include "macoctreesizingfunc3.h"
//
SHKZ_BEGIN_NAMESPACE
//
namespace macotreeliquid3_namespace {
	//
	class macoctreeliquid3 : public drawable {
	public:
		//
		LONG_NAME("MAC Octree Liquid 3D")
		ARGUMENT_NAME("macoctreeliquid")
		//
	protected:
		//
		shape3 m_shape;
		double m_dx;
		double m_narrowband_depth;
		double m_initial_volume, m_current_volume, m_y_prev;
		//
		shape3 m_solid_shape;
		double m_solid_dx;
		//
		gridvisualizer3_driver m_solid_gridvisualizer{this,"gridvisualizer3"};
		array3<Real> m_solid_visualize{this};
		//
		grid3 m_grid_0{this};
		grid3 m_grid_1{this};
		//
		grid3 *m_grid {&m_grid_0};
		grid3 *m_grid_prev {&m_grid_1};
		//
		// -------- Reginal Volume Correction Members --------
		std::vector<Real> m_y_list;
		uint_type m_region_count;
		std::vector<uint_type> m_regions;
		std::vector<Real> m_current_volumes;
		std::vector<Real> m_volumes;
		// ---------------------------------------------------
		//
		macflip3_driver m_flip{this,"macexnbflip3"};
		//
		virtual void setup_window( std::string &name, int &width, int &height ) const override;
		virtual void load( configuration &config ) override;
		virtual void configure( configuration &config ) override;
		virtual void post_initialize( bool initialized_from_file ) override;
		//
		virtual void drag( double x, double y, double z, double u, double v, double w ) override;
		virtual void idle() override;
		virtual void draw( graphics_engine &g ) const override;
		virtual bool should_quit() const override;
		virtual bool should_screenshot() const override { return m_timestepper->should_export_frame(); }
		//
		void export_moving_polygon();
		void do_export_solid_mesh( const array3<Real> &solid );
		void do_export_empty_solid_mesh( bool force=false );
		void export_mesh( int frame );
		void render_mesh( unsigned frame ) const;
		void save_state();
		//
		timestepper_driver m_timestepper{this,"timestepper"};
		parallel_driver m_parallel{this};
		dylibloader m_dylib;
		bool m_should_quit_on_save {false}; // should not be saved
		//
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
		macoctreeproject3 m_macoctreeproject{this};
		macoctreesegregator3 m_macoctreesegregator{this};
		macoctreemesher3 m_macoctreemesher{this};
		meshutility3_driver m_meshutility{this,"meshutility3"};
		macoctreesizingfunc3 m_macoctreesizingfunc{this};
		meshexporter3_driver m_mesh_exporter{this,"meshexporter3"};
		cellmesher3_driver m_solid_mesher{this,"marchingcubes"};
		graphics_interface_driver m_svg_writer{this,"graphics_svg"};
		//
		environment_setter arg_shape{this,"shape",&m_shape};
		environment_setter arg_dx{this,"dx",&m_dx};
		//
		struct Parameters {
			bool use_FLIP {true};
			unsigned min_resolution {16};
			unsigned erode_width {0};
			double surftens_k {0.0};
			bool use_sizing_func {true};
			unsigned initial_refinement {3};
			double maximal_CFL_accumulation {1.0};
			bool render_mesh {false};
			bool render_wireframe {false};
			bool render_grid {false};
			bool remove_quater {false};
			double z {0.5};
			bool export_svg {true};
			unsigned render_sample_count {8};
			unsigned save_interval {100};
			vec3d target {0.5,0.15,0.5};
			vec3d origin {0.5,1.5,3.0};
			vec3d gravity {0.0,-9.8,0.0};
			double PICFLIP {0.98};
			bool render_transparent {false};
			bool volume_correction {true};
			bool regional_volume_correction {false};
			bool maccormack {false};
			unsigned render_transparent_sample_count {32};
			bool transfer_file {true};
			int debug_mode {0};
		};
		//
		Parameters m_param;
		//
		std::function<bool( double dx, double dt, double time, unsigned step )> m_check_inject_func;
		std::function<bool( const vec3d &p, double dx, double dt, double time, unsigned step, double &fluid, vec3d &velocity )> m_inject_func;
		std::function<void( double dx, double dt, double time, unsigned step, double &volume_change )> m_post_inject_func;
		//
		bool m_do_inject;
		double m_injected_volume;
		double m_accumulated_CFL;
		//
		virtual void begin_inject_external_fluid( double dt, double time, unsigned step );
		virtual void do_inject_external_fluid( double dt, double time, unsigned step );
		virtual void end_inject_external_fluid( double dt, double time, unsigned step );
		//
		virtual void initialize( const filestream &file ) override {
			file.r(m_shape);
			file.r(m_dx);
			file.r(m_narrowband_depth);
			file.r(m_initial_volume);
			file.r(m_current_volume);
			file.r(m_y_prev);
			file.r(m_solid_shape);
			file.r(m_solid_dx);
			file.read(m_y_list);
			file.r(m_region_count);
			file.read(m_regions);
			file.read(m_current_volumes);
			file.read(m_volumes);
			file.r(m_do_inject);
			file.r(m_injected_volume);
			file.r(m_accumulated_CFL);
		}
		virtual void serialize( const filestream &file ) const override {
			file.w(m_shape);
			file.w(m_dx);
			file.w(m_narrowband_depth);
			file.w(m_initial_volume);
			file.w(m_current_volume);
			file.w(m_y_prev);
			file.w(m_solid_shape);
			file.w(m_solid_dx);
			file.write(m_y_list);
			file.w(m_region_count);
			file.write(m_regions);
			file.write(m_current_volumes);
			file.write(m_volumes);
			file.w(m_do_inject);
			file.w(m_injected_volume);
			file.w(m_accumulated_CFL);
		}
	};
};
//
SHKZ_END_NAMESPACE
//
#endif
//