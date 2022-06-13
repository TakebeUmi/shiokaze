/*
**	macoctreeaccuracytest3.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on December 30, 2019. All rights reserved.
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
#ifndef SHKZ_OCTREESTANDINGWATERTEST3_H
#define SHKZ_OCTREESTANDINGWATERTEST3_H
//
#define _USE_MATH_DEFINES
#include <shiokaze/ui/drawable.h>
#include <shiokaze/graphics/graphics_interface.h>
#include <shiokaze/core/filesystem.h>
#include <shiokaze/core/console.h>
#include <shiokaze/meshexporter/meshexporter3_interface.h>
#include <utility>
#include <cmath>
//
#include "macoctreemesher3.h"
#include "macoctreegrid3.h"
#include "macoctreeproject3.h"
#include "mitsuba_xml.h"
//
SHKZ_BEGIN_NAMESPACE
//
namespace macotreeliquid3_namespace {
	//
	class macoctreeaccuracytest3 : public drawable {
	protected:
		//
		LONG_NAME("MAC Octree Standing Water Test 3D")
		ARGUMENT_NAME("macoctreeaccuracytest")
		//
		virtual void configure( configuration &config ) override {
			//
			m_dx = m_shape.dx();
			//
			set_environment("shape",&m_shape);
			set_environment("dx",&m_dx);
			//
			config.set_default_unsigned("CanvasX",1280);
			config.set_default_unsigned("CanvasY",640);
			config.set_default_integer("Threads",1);
			config.set_default_double("Residual",1e-18);
			config.set_default_double("EpsFluid",0.0);
			config.set_default_double("AccuracyEps",0.0);
			config.set_default_bool("SteepAdaptivity",false);
			config.set_default_bool("DegradeWarning",true);
			config.set_default_bool("VolumeCorrection",false);
			config.set_default_unsigned("DilateCount",0);
			config.set_default_unsigned("MaxIterations",300000);
			//
			config.get_bool("Debug",m_param.debug,"Debug mode");
			config.get_bool("OnlyVertical",m_param.only_vertical,"Allow only vertical movement");
			m_param.debug = m_param.debug || config.is_parameter_set("LinSolver");
			config.get_double("Level",m_param.water_level,"Water level");
			config.get_double("Scale",m_param.scale,"Force scale");
			config.get_bool("RenderMesh",m_param.render_mesh,"Render mesh");
			config.get_unsigned("SampleCount",m_param.sample_count,"Sample count");
			config.get_unsigned("MaxFrame",m_param.max_frame,"Maxiam frame");
		}
		//
		virtual void post_initialize( bool initialized_from_file ) override {
			//
			m_step = 0;
			m_time = 0.0;
			m_grid.clear();
			unsigned n (1);
			for( unsigned k=0; k<4; ++k ) {
				shape3 shape = m_shape/n;
				m_grid.add_layer(shape,n*m_dx);
				n *= 2;
			}
			auto adaptivity_func = [&](char depth, const vec3d &p) {
				return 4*p[0] < depth+1;
			};
			m_grid.activate_cells(adaptivity_func);
			m_grid.balance_layers();
			m_grid.assign_indices();
			//
			m_camera->set_bounding_box(vec3d().v,m_shape.box(m_dx).v);
			idle();
		}
		//
		virtual void idle() override {
			//
			m_time += 0.01;
			vec3d tilt_vec;
			if( m_param.only_vertical ) {
				tilt_vec = vec3d(0.0,1.0,0.0);
			} else {
				tilt_vec = vec3d(0.1*cos(3.0*m_time),1.0,0.1*sin(3.0*m_time)).normal();
			}
			//
			auto fluid_func = [&](const vec3d &p) { return (tilt_vec*p)-m_param.water_level+0.05*sin(m_time); };
			auto solid_func = [&](const vec3d &p) { return 1.0; };
			//
			m_grid.assign_levelset(fluid_func,solid_func);
			m_grid.set_velocity([&]( const vec3d &p, char dim ) { return -m_param.scale*tilt_vec[dim]; });
			if( m_param.debug ) {
				m_macoctreeproject.assemble_matrix(m_grid);
			}
			m_macoctreeproject.project(m_grid,1.0);
			//
			m_velocity.clear();
			m_grid.serial_iterate_active_cells([&]( const cell_id3 &cell_id ) {
				//
				if( m_grid.levelset[cell_id.index] < 0.0 ) {
					const vec3d &p = m_grid.get_cell_position(cell_id);
					const vec3d &u = m_grid.sample_velocity(p);
					if( u.len() > 1e-2 ) {
						m_velocity.push_back(std::make_pair(p,u));
					}
				}
			});
			//
			m_has_velocity.push_back(m_velocity.empty() ? false : true);
			if( console::get_root_path().size()) {
				export_mesh(m_step);
				if( m_param.render_mesh ) render_mesh(m_step);
			}
			m_step ++;
			if( m_step > m_param.max_frame+10 ) exit(0);
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
				"		<lookat target=\"0.5,0.15,0.5\" origin=\"$x,$y,$z\" up=\"0,1,0\"/>\n"
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
				"			<srgb name=\"interiorColor\" value=\"%g, %g, %g\"/>\n"
				"			<srgb name=\"edgeColor\" value=\"0.0\"/>\n"
				"		</texture>\n"
				"		<float name=\"intIOR\" value=\"1.2\"/>\n"
				"	</bsdf></bsdf>\n"
				"</shape>\n"
				};
				//
				const std::string arrow_mesh { 
				"<shape type=\"serialized\">\n"
				"	<string name=\"filename\" value=\"$arrowname\"/>\n"
				"	<bsdf type=\"diffuse\">\n"
				"		<srgb name=\"reflectance\" value=\"#D1F000\"/>\n"
				"	</bsdf>\n"
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
				fprintf(file,wireframe_mesh.c_str(),0.5,0.6,1.0);
				fprintf(file,arrow_mesh.c_str());
				fprintf(file,mitsuba_xml::footer.c_str());
				fclose(file);
				//
				file = fopen((export_path+"/arrow.xml").c_str(),"w");
				fprintf(file,mitsuba_xml::header.c_str());
				fprintf(file,sensor.c_str(),m_param.sample_count);
				fprintf(file,sunlight.c_str());
				fprintf(file,mitsuba_xml::floor.c_str());
				fprintf(file,arrow_mesh.c_str());
				fprintf(file,mitsuba_xml::footer.c_str());
				fclose(file);
			}
			//
			// Export fluid mesh
			std::vector<vec3d> vertices;
			std::vector<std::vector<size_t> > faces;
			m_macoctreemesher.generate_mesh(m_grid,m_grid.levelset,0.0,nullptr,vertices,faces,false);
			//
			m_mesh_exporter->set_mesh(vertices,faces);
			m_mesh_exporter->export_ply(console::format_str("%s/%d_mesh.ply",export_path.c_str(),frame));
			m_mesh_exporter->export_mitsuba(console::format_str("%s/%d_mesh.serialized",export_path.c_str(),frame));
			//
			// Export arrow mesh
			std::vector<vec3d> arrow_vertices;
			std::vector<std::vector<size_t> > arrow_faces;
			//
			const double r (0.0025);
			for( const auto &e : m_velocity ) {
				//
				const vec3d &dir_vec = e.second;
				const double len = dir_vec.len();
				vec3d normal = dir_vec ^ vec3d(1.0,0.0,0.0);
				if( normal.empty()) normal = dir_vec ^ vec3d(0.0,1.0,0.0);
				normal.normalize();
				//
				const vec3d &x_vec = normal;
				const vec3d &y_vec = (normal ^ dir_vec).normal();
				//
				const unsigned count (10);
				const double theta = 2.0 * M_PI / count;
				const double tip0 (0.8);
				const double tip1 (3.0);
				//
				size_t prev_indices0[3];
				size_t prev_indices[3];
				size_t tip_indices[2];
				//
				const vec3d &p0 = e.first;
				const vec3d &p4 = p0+dir_vec;
				//
				tip_indices[0] = arrow_vertices.size();
				tip_indices[1] = tip_indices[0]+1;
				arrow_vertices.push_back(p0);
				arrow_vertices.push_back(p4);
				//
				for( unsigned n=0; n<count; ++n ) {
					//
					const double t = n*theta;
					const vec3d &nvec = cos(t)*x_vec+sin(t)*y_vec;
					const vec3d &p1 = p0+r*nvec;
					const vec3d &p2 = p1+tip0*dir_vec;
					const vec3d &p3 = p2+tip1*r*nvec;
					//
					size_t indices[3];
					indices[0] = arrow_vertices.size();
					indices[1] = indices[0]+1;
					indices[2] = indices[0]+2;
					//
					if( n >= 1 ) {
						arrow_faces.push_back({indices[0],prev_indices[0],tip_indices[0]});
						arrow_faces.push_back({prev_indices[0],indices[0],indices[1],prev_indices[1]});
						arrow_faces.push_back({prev_indices[1],indices[1],indices[2],prev_indices[2]});
						arrow_faces.push_back({prev_indices[2],indices[2],tip_indices[1]});
					}
					//
					if( n == count-1 ) {
						arrow_faces.push_back({prev_indices0[0],indices[0],tip_indices[0]});
						arrow_faces.push_back({prev_indices0[1],indices[1],indices[0],prev_indices0[0]});
						arrow_faces.push_back({prev_indices0[2],indices[2],indices[1],prev_indices0[1]});
						arrow_faces.push_back({indices[2],prev_indices0[2],tip_indices[1]});
					}
					//
					prev_indices[0] = indices[0];
					prev_indices[1] = indices[1];
					prev_indices[2] = indices[2];
					//
					if( n == 0 ) {
						prev_indices0[0] = prev_indices[0];
						prev_indices0[1] = prev_indices[1];
						prev_indices0[2] = prev_indices[2];
					}
					//
					arrow_vertices.push_back(p1);
					arrow_vertices.push_back(p2);
					arrow_vertices.push_back(p3);
				}
			}
			//
			if( arrow_vertices.empty()) {
				arrow_vertices.push_back(vec3d(-1e4+1e-4,-1e4,-1e4));
				arrow_vertices.push_back(vec3d(-1e4,-1e4+1e-4,-1e4));
				arrow_vertices.push_back(vec3d(-1e4,-1e4,-1e4));
				arrow_faces.push_back({0,1,2});
			}
			//
			for( const auto &f : arrow_faces ) {
				for( const auto &index : f ) assert( index < arrow_vertices.size());
			}
			m_mesh_exporter->set_mesh(arrow_vertices,arrow_faces);
			m_mesh_exporter->export_ply(console::format_str("%s/%d_arrow.ply",export_path.c_str(),frame));
			m_mesh_exporter->export_mitsuba(console::format_str("%s/%d_arrow.serialized",export_path.c_str(),frame));
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
				if( step >= 10 ) {
					for( int step1=0; step1<=step-10; ++step1 ) {
						//
						double distance_exist_velocity_frame (1.0);
						for( int i=-10; i<=10; ++i ) {
							int frame_examine = step1+i;
							if( frame_examine >= 10 && m_has_velocity[frame_examine]) {
								distance_exist_velocity_frame = std::min(distance_exist_velocity_frame,std::abs(i/10.0));
							}
						}
						//
						if( ! filesystem::is_exist(console::format_str("%s/%u_mesh.exr",img_path.c_str(),step1))) {
							std::string command = console::format_str("mitsuba -Dx=0.5 -Dy=0.85 -Dz=3.5 -Dw=1280 -Dh=720 -Darrowname=%s/%u_arrow.serialized -Dfilename=%s/%u_mesh.serialized -o %s/%u_mesh.exr %s/mesh.xml",
								export_path.c_str(),step1,export_path.c_str(),step1,img_path.c_str(),step1,export_path.c_str());
							//
							console::dump("Running command: %s\n", command.c_str());
							console::system(command);
						}
						if( ! filesystem::is_exist(console::format_str("%s/%u_mesh.png",img_path.c_str(),step1))) {
							console::system("mtsutil tonemap -g 3 %s/%u_mesh.exr",img_path.c_str(),step1);
						}
						//
						if( ! filesystem::is_exist(console::format_str("%s/%u_arrow.exr",img_path.c_str(),step1))) {
							std::string command = console::format_str("mitsuba -Dx=0.5 -Dy=0.85 -Dz=3.5 -Dw=1280 -Dh=720 -Darrowname=%s/%u_arrow.serialized -o %s/%u_arrow.exr %s/arrow.xml",
								export_path.c_str(),step1,img_path.c_str(),step1,export_path.c_str());
							//
							console::dump("Running command: %s\n", command.c_str());
							console::system(command);
						}
						if( ! filesystem::is_exist(console::format_str("%s/%u_arrow.png",img_path.c_str(),step1))) {
							console::system("mtsutil tonemap -g 3 %s/%u_arrow.exr",img_path.c_str(),step1);
						}
						//
						std::string exp = img_path;
						std::string final_path = console::format_str("%s/%d_composite.jpg",exp.c_str(),step1);
						if( ! filesystem::is_exist(final_path)) {
							const double blending = 100.0 * distance_exist_velocity_frame + (1.0-distance_exist_velocity_frame) * 70.0;
							console::system("composite -blend %d %s/%d_mesh.png %s/%d_arrow.png %s",
									(int)blending, exp.c_str(),step1,exp.c_str(),step1,final_path.c_str());
						}
					}
					//
					if( step && step % 10 == 0 ) {
						console::system("rm -f %s/composite.mp4; avconv -r 60 -i %s/%%d_composite.jpg -b:v 12000k -pix_fmt yuv420p %s/composite.mp4",
									img_path.c_str(), img_path.c_str(), img_path.c_str() );
					}
				}
			}
		}
		//
		virtual void draw( graphics_engine &g ) const override {
			//
			// Draw surface
			std::vector<vec3d> vertices;
			std::vector<std::vector<size_t> > faces;
			m_macoctreemesher.generate_mesh(m_grid,m_grid.levelset,0.0,nullptr,vertices,faces,false);
			//
			g.color4(1.0,1.0,1.0,0.5);
			for( unsigned i=0; i<faces.size(); i++ ) {
				g.begin(graphics_engine::MODE::LINE_LOOP);
				for( unsigned j=0; j<faces[i].size(); j++ ) g.vertex3v(vertices[faces[i][j]].v);
				g.end();
			}
			//
			// Draw velocity
			g.color4(1.0,1.0,0.25,0.5);
			g.begin(graphics_engine::MODE::LINES);
			for( const auto &e : m_velocity ) {
				const vec3d &p = e.first;
				const vec3d &u = e.second;
				g.vertex3v(p.v);
				g.vertex3v((p+u).v);
			}
			g.end();
		}
		//
		virtual void setup_window( std::string &name, int &width, int &height ) const override {
			height = width;
		}
		//
		shape3 m_shape {64,32,64};
		double m_dx;
		double m_time;
		unsigned m_step;
		std::vector<bool> m_has_velocity;
		//
		struct Parameters {
			bool debug {true};
			bool only_vertical {false};
			double water_level {0.25};
			double scale {5.0};
			bool render_mesh {false};
			unsigned sample_count {8};
			unsigned max_frame {420};
		};
		Parameters m_param;
		//
		grid3 m_grid{this};
		macoctreeproject3 m_macoctreeproject{this};
		macoctreemesher3 m_macoctreemesher{this};
		meshexporter3_driver m_mesh_exporter{this,"meshexporter3"};
		//
		std::vector<std::pair<vec3d,vec3d> > m_velocity;
	};
	//
	extern "C" module * create_instance() {
		return new macoctreeaccuracytest3;
	}
};
//
SHKZ_END_NAMESPACE
//
#endif