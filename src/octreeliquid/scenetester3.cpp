/*
**	scenetester3.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on January 10, 2020.
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
#include <shiokaze/ui/drawable.h>
#include <shiokaze/math/shape.h>
#include <shiokaze/core/dylibloader.h>
#include <shiokaze/utility/utility.h>
#include <shiokaze/core/console.h>
//
SHKZ_USING_NAMESPACE
//
class scenetester3 : public drawable {
public:
	//
	LONG_NAME("Scene Tester 3D")
	ARGUMENT_NAME("SceneTester")
	//
	virtual void load( configuration &config ) override {
		//
		std::string name("waterdrop3"); config.get_string("Name",name,"Scene file name");
		m_dylib.open_library(filesystem::resolve_libname(name));
		m_dylib.load(config);
		m_dylib.overwrite(config);
	}
	//
	virtual void configure( configuration &config ) override {
		//
		m_dylib.configure(config);
		//
		config.get_unsigned("ResolutionX",m_shape[0],"Resolution towards X axis");
		config.get_unsigned("ResolutionY",m_shape[1],"Resolution towards Y axis");
		config.get_unsigned("ResolutionZ",m_shape[2],"Resolution towards Y axis");
		//
		double resolution_scale (1.0);
		config.get_double("ResolutionScale",resolution_scale,"Resolution doubling scale");
		//
		double view_scale (1.0);
		config.get_double("ViewScale",view_scale,"View scale");
		//
		m_shape *= resolution_scale;
		m_dx = view_scale * m_shape.dx();
		//
		set_environment("shape",&m_shape);
		set_environment("dx",&m_dx);
	}
	//
	virtual void post_initialize( bool initialized_from_file ) override {
		//
		auto initialize_func = reinterpret_cast<void(*)(const shape3 &m_shape, double m_dx)>(m_dylib.load_symbol("initialize"));
		if( initialize_func ) {
			initialize_func(m_shape,m_dx);
		}
		m_draw_func = reinterpret_cast<void(*)(graphics_engine &,double)>(m_dylib.load_symbol("draw"));
		m_camera->set_bounding_box(vec3d().v,m_shape.box(m_dx).v);
		//
		m_time0 = utility::get_seconds();
	}
	//
	virtual bool keyboard( int key, int action, int mods ) override {
		return false;
	}
	//
	virtual void draw( graphics_engine &g ) const override {
		//
		const float seconds = utility::get_seconds()-m_time0;
		//
		if( m_draw_func ) {
			m_draw_func(g,seconds);
		}
		g.set_2D_coordinate(0.0,1.0,0,1.0);
		g.color4(1.0,1.0,1.0,1.0);
		static double test (0.1);
		std::string str = console::format_str("Time: %.2f seconds",seconds);
		g.draw_string(vec2d(0.02,0.95).v,str,18);
	}
	//
	virtual void setup_window( std::string &name, int &width, int &height ) const override {
		double ratio = m_shape[1] / (double) m_shape[0];
		height = ratio * width;
	}
	//
	shape3 m_shape {32,32,32};
	double m_dx;
	double m_time0;
	//
	std::function<void(graphics_engine &,double)> m_draw_func;
	dylibloader m_dylib;
};
//
extern "C" module * create_instance() {
	return new scenetester3;
}
//
extern "C" const char *license() {
	return "MIT";
}
//