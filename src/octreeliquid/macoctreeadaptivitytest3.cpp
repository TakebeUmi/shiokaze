/*
**	macoctreeadaptivitytest3.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on January 1, 2020. All rights reserved.
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
#ifndef SHKZ_OCTREEADAPTIVITYTEST3_H
#define SHKZ_OCTREEADAPTIVITYTEST3_H
//
#include <shiokaze/ui/drawable.h>
#include <shiokaze/core/filesystem.h>
#include <shiokaze/core/console.h>
#include <shiokaze/meshexporter/meshexporter3_interface.h>
//
#include "macoctreesizingfunc3.h"
#include "macoctreegrid3.h"
#include "macoctreemesher3.h"
#include "mitsuba_xml.h"
//
SHKZ_BEGIN_NAMESPACE
//
namespace macotreeliquid3_namespace {
	//
	class macoctreeadaptivitytest3 : public drawable {
	protected:
		//
		LONG_NAME("MAC Octree Adaptivity Test 3D")
		ARGUMENT_NAME("macadaptivitytest")
		//
		virtual void configure( configuration &config ) override {
			//
			double scale (1.0);
			config.get_double("ResolutionScale",scale,"Resolution doubling scale");
			config.get_integer("Scene",m_param.scene,"Scene");
			config.get_unsigned("MinResolution",m_param.min_resolution,"Minimal resolution");
			config.get_bool("ReuseGrid",m_param.reuse_grid,"Reuse grid");
			config.get_bool("RenderMesh",m_param.render_mesh,"Render mesh");
			config.get_unsigned("SampleCount",m_param.sample_count,"Sample count");
			//
			m_shape *= scale;
			m_dx = m_shape.dx();
			//
			set_environment("shape",&m_shape);
			set_environment("dx",&m_dx);
			//
			config.set_default_bool("SteepAdaptivity",false);
		}
		//
		virtual void post_initialize( bool initialized_from_file ) override {
			//
			m_step = 0;
			//
			refine([&]() {
				m_grid->activate_cells([&](char depth, const vec3d &p) {
					return depth > 3;
				});
			});
			update(false);
			m_camera->set_bounding_box(vec2d().v,m_shape.box(m_dx).v);
		}
		//
		void export_mesh( unsigned frame ) {
			//
			std::string export_path = console::get_root_path() + "/mesh";
			if( ! filesystem::is_exist(export_path)) {
				filesystem::create_directory(export_path);
				//
				const std::string sensor {
				"<sensor type=\"perspective\">\n"
				"	<float name=\"focusDistance\" value=\"6\"/>\n"
				"	<float name=\"fov\" value=\"30\"/>\n"
				"	<string name=\"fovAxis\" value=\"x\"/>\n"
				"	<transform name=\"toWorld\">\n"
				"		<lookat target=\"0.5,0.4,0.5\" origin=\"$x,$y,$z\" up=\"0,1,0\"/>\n"
				"	</transform>\n"
				"	<sampler type=\"halton\">\n"
				"		<integer name=\"sampleCount\" value=\"%d\"/>\n"
				"	</sampler>\n"
				"	<film type=\"hdrfilm\">\n"
				"		<boolean name=\"banner\" value=\"false\"/>\n"
				"		<integer name=\"width\" value=\"$w\"/>\n"
				"		<integer name=\"height\" value=\"$h\"/>\n"
				"		<string name=\"pixelFormat\" value=\"rgb\"/>\n"
				"		<rfilter type=\"gaussian\"/>\n"
				"	</film>\n"
				"</sensor>\n"
				};
				//
				const std::string wireframe_mesh { 
				"<shape type=\"serialized\">\n"
				"	<string name=\"filename\" value=\"$filename\"/>\n"
				"	<bsdf type=\"twosided\"> <bsdf type=\"plastic\">\n"
				"		<texture type=\"wireframe\" name=\"diffuseReflectance\">\n"
				"			<float name=\"lineWidth\" value=\"0.0005\"/>\n"
				"			<srgb name=\"interiorColor\" value=\"%g, %g, %g\"/>\n"
				"			<srgb name=\"edgeColor\" value=\"0.0\"/>\n"
				"		</texture>\n"
				"		<float name=\"intIOR\" value=\"1.2\"/>\n"
				"	</bsdf></bsdf>\n"
				"</shape>\n"
				};
				//
				const std::string sunlight {
				"<emitter type=\"sunsky\">\n"
				"	<spectrum name=\"albedo\" value=\"0\"/>\n"
				"	<vector name=\"sunDirection\" x=\"-0.2\" y=\"0.4\" z=\"0.3\"/>\n"
				"	<float name=\"sunScale\" value=\"2.0\"/><float name=\"skyScale\" value=\"6.0\"/>\n"
				"</emitter>\n"
				};
				//
				// Export XML file
				FILE * file = fopen((export_path+"/mesh.xml").c_str(),"w");
				fprintf(file,mitsuba_xml::header.c_str());
				fprintf(file,sensor.c_str(),m_param.sample_count);
				fprintf(file,sunlight.c_str());
				fprintf(file,mitsuba_xml::floor.c_str());
				fprintf(file,wireframe_mesh.c_str(),0.7,0.8,1.0);
				fprintf(file,mitsuba_xml::footer.c_str());
				fclose(file);
			}
			//
			// Export fluid mesh
			std::vector<vec3d> vertices;
			std::vector<std::vector<size_t> > faces;
			m_macoctreemesher.generate_mesh(*m_grid,m_grid->levelset,0.0,nullptr,vertices,faces,false);
			//
			m_mesh_exporter->set_mesh(vertices,faces);
			m_mesh_exporter->export_ply(console::format_str("%s/%d_mesh.ply",export_path.c_str(),frame));
			m_mesh_exporter->export_mitsuba(console::format_str("%s/%d_mesh.serialized",export_path.c_str(),frame));
		}
		//
		void render_mesh( int step ) {
			//
			if( m_param.render_mesh && console::system("mitsuba > /dev/null 2>&1") == 0 ) {
				//
				const std::string export_path = console::get_root_path() + "/mesh";
				assert(filesystem::is_exist(export_path));
				//
				const std::string img_path = console::get_root_path() + "/img";
				if( ! filesystem::is_exist(img_path)) {
					filesystem::create_directory(img_path);
				}
				//
				for( int step1=0; step1<=step; ++step1 ) {
					//
					if( ! filesystem::is_exist(console::format_str("%s/%u_mesh.exr",img_path.c_str(),step1))) {
						std::string command = console::format_str("mitsuba -Dx=0.5 -Dy=1.65 -Dz=3.5 -Dw=1280 -Dh=720 -Dfilename=%s/%u_mesh.serialized -o %s/%u_mesh.exr %s/mesh.xml",
							export_path.c_str(),step1,img_path.c_str(),step1,export_path.c_str());
						//
						console::dump("Running command: %s\n", command.c_str());
						console::system(command);
					}
					if( ! filesystem::is_exist(console::format_str("%s/%u_mesh.png",img_path.c_str(),step1))) {
						console::system("mtsutil tonemap -g 3 %s/%u_mesh.exr",img_path.c_str(),step1);
					}
				}
			}
		}
		//
		void refine( std::function<void()> activate_func ) {
			//
			std::swap(m_grid,m_grid_prev);
			m_grid->clear();
			shape3 shape(m_shape);
			double dx (m_dx);
			while( shape.max() >= m_param.min_resolution ) {
				m_grid->add_layer(shape,dx);
				shape = shape/2;
				dx *= 2.0;
			}
			activate_func();
			m_grid->balance_layers();
			m_grid->assign_indices();
		}
		//
		void update( bool use_prev_grid ) {
			//
			auto solid_func = [&](const vec3d &p) { return 1.0; };
			//
			if( use_prev_grid ) {
				m_grid->assign_levelset([&](const vec3d &p) {
					return m_grid_prev->sample_levelset(p);
				},solid_func);
			} else {
				if( m_param.scene == 0 ) {
					m_grid->assign_levelset([&](const vec3d &p) {
						double value = (p-vec3d(0.5,0.5,0.5)).len()-0.2;
						value = std::min(value,p[1]-0.51);
						return value;
					},solid_func);
				} else if( m_param.scene == 1 ) {
					m_grid->assign_levelset([&](const vec3d &p) {
						double value = (p-vec3d(0.5,0.7,0.5)).len()-0.2;
						value = std::min(value,p[1]-0.51);
						return value;
					},solid_func);
				}
			}
			//
			if( m_param.scene == 0 ) {
				m_grid->set_velocity([&]( const vec3d &p, char dim ) { return 0.0; });
			} else if( m_param.scene == 1 ) {
				m_grid->set_velocity([&]( const vec3d &p, char dim ) { return p[1]<0.5 ? 0.0 : vec2d(0.0,-1.0)[dim]; });
			}
			//
			m_macoctreesizingfunc.compute_sizing_function(*m_grid_prev,*m_grid,0.0);
			if( console::get_root_path().size()) {
				export_mesh(m_step);
				render_mesh(m_step);
			}
			m_step ++;
		}
		//
		virtual bool keyboard( int key, int action, int mods ) override {
			//
			if( action == UI_interface::PRESS && key == UI_interface::KEY_Q ) {
				refine([&]() {
					m_macoctreesizingfunc.activate_cells(*m_grid_prev,*m_grid);
				});
				update(m_param.reuse_grid);
				return true;
			}
			return false;
		}
		//
		virtual void draw( graphics_engine &g ) const override {
			//
			// Draw surface
			std::vector<vec3d> vertices;
			std::vector<std::vector<size_t> > faces;
			m_macoctreemesher.generate_mesh(*m_grid,m_grid->levelset,0.0,nullptr,vertices,faces,false);
			//
			g.color4(1.0,1.0,1.0,0.5);
			for( unsigned i=0; i<faces.size(); i++ ) {
				g.begin(graphics_engine::MODE::LINE_LOOP);
				for( unsigned j=0; j<faces[i].size(); j++ ) g.vertex3v(vertices[faces[i][j]].v);
				g.end();
			}
		}
		//
		virtual void setup_window( std::string &name, int &width, int &height ) const override {
			height = width;
		}
		//
		shape3 m_shape {64,64,64};
		double m_dx;
		unsigned m_step;
		//
		std::function<double(const vec3d &p)> m_fluid_func;
		//
		struct Parameters {
			int scene {0};
			bool reuse_grid {false};
			unsigned min_resolution {8};
			bool render_mesh {false};
			unsigned sample_count {64};
		};
		Parameters m_param;
		//
		grid3 m_grid0{this};
		grid3 m_grid1{this};
		//
		grid3 *m_grid {&m_grid0};
		grid3 *m_grid_prev {&m_grid1};
		//
		macoctreesizingfunc3 m_macoctreesizingfunc{this};
		macoctreemesher3 m_macoctreemesher{this};
		meshexporter3_driver m_mesh_exporter{this,"meshexporter3"};
		std::vector<Real> lap_levelset;
	};
	//
	extern "C" module * create_instance() {
		return new macoctreeadaptivitytest3;
	}
};
//
SHKZ_END_NAMESPACE
//
#endif