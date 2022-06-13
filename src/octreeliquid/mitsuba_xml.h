/*
**	mitsuba_xml.h
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on November 17, 2019. All rights reserved.
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
#ifndef SHKZ_MITSUBA_XML_H
#define SHKZ_MITSUBA_XML_H
//
#include <shiokaze/core/common.h>
//
SHKZ_BEGIN_NAMESPACE
//
namespace mitsuba_xml {
	//
	std::string header {
	"<?xml version='1.0' encoding='utf-8'?>\n"
	"<scene version=\"0.5.0\">\n"
	"<integrator type=\"direct\">\n"
	"</integrator>\n"
	};
	//
	std::string sensor {
	"<sensor type=\"perspective\">\n"
	"	<float name=\"focusDistance\" value=\"6\"/>\n"
	"	<float name=\"fov\" value=\"30\"/>\n"
	"	<string name=\"fovAxis\" value=\"x\"/>\n"
	"	<transform name=\"toWorld\">\n"
	"		<lookat target=\"0.5,0.45,0.5\" origin=\"$x,$y,$z\" up=\"0,1,0\"/>\n"
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
	std::string sunlight {
	"<emitter type=\"sunsky\">\n"
	"	<spectrum name=\"albedo\" value=\"0\"/>\n"
	"	<vector name=\"sunDirection\" x=\"-0.2\" y=\"0.4\" z=\"0.3\"/>\n"
	"	<float name=\"sunScale\" value=\"2.0\"/><float name=\"skyScale\" value=\"6.0\"/>\n"
	"</emitter>\n"
	"<shape type=\"sphere\">\n"
	"	<point name=\"center\" x=\"0.8\" y=\"0.5\" z=\"3.0\"/>\n"
	"	<float name=\"radius\" value=\"0.5\"/>\n"
	"	<emitter type=\"area\">\n"
	"		<spectrum name=\"radiance\" value=\"20\"/>\n"
	"	</emitter>\n"
	"</shape>\n"
	};
	//
	std::string inlight {
	"<shape type=\"sphere\">\n"
	"	<point name=\"center\" x=\"0.4\" y=\"0.5\" z=\"0.8\"/>\n"
	"	<float name=\"radius\" value=\"0.03\"/>\n"
	"	<emitter type=\"area\">\n"
	"		<spectrum name=\"radiance\" value=\"50\"/>\n"
	"	</emitter>\n"
	"</shape>\n"
	};
	//
	std::string floor {
	"<shape type=\"rectangle\">\n"
	"	<transform name=\"toWorld\">\n"
	"		<rotate x=\"1\" angle=\"-90\"/>\n"
	"		<scale value=\"100\"/>\n"
	"	</transform>\n"
	"	<bsdf type=\"diffuse\">\n"
	"		<srgb name=\"reflectance\" value=\"#FFFFFF\"/>\n"
	"	</bsdf>\n"
	"</shape>\n"
	};
	//
	std::string mesh { 
	"<shape type=\"serialized\">\n"
	"	<string name=\"filename\" value=\"%s\"/>\n"
	"	<bsdf type=\"twosided\"> <bsdf type=\"plastic\">\n"
	"		<spectrum name=\"diffuseReflectance\" value=\"%g, %g, %g\"/>\n"
	"		<float name=\"intIOR\" value=\"1.2\"/>\n"
	"	</bsdf></bsdf>\n"
	"</shape>\n"
	};
	//
	std::string wireframe_mesh { 
	"<shape type=\"serialized\">\n"
	"	<string name=\"filename\" value=\"%s\"/>\n"
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
	std::string sphere {
	"<shape type=\"sphere\">\n"
	"	<point name=\"center\" x=\"%g\" y=\"%g\" z=\"%g\"/>\n"
	"	<float name=\"radius\" value=\"%g\"/>\n"
	"	<bsdf type=\"diffuse\">\n"
	"		<srgb name=\"reflectance\" value=\"#333333\"/>\n"
	"	</bsdf>\n"
	"</shape>\n"
	};
	//
	std::string footer {"</scene>\n"};
};
//
SHKZ_END_NAMESPACE
//
#endif
//