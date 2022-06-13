/*
**	octree2.h
**
**	This is part of Shiokaze fluid solver, a research-oriented fluid solver designed for collaborative projects.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on April 10, 2017. All rights reserved.
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
#ifndef SHKZ_OCTREE2_H
#define SHKZ_OCTREE2_H
//
#include <shiokaze/math/vec.h>
#include <shiokaze/graphics/graphics_engine.h>
#include <functional>
#include <vector>
//
SHKZ_BEGIN_NAMESPACE
//
class octree2 {
public:
	// Octree cell structure
	typedef struct _leaf2 {
		_leaf2 *children[2][2];
		vec2i position;
		bool subdivided;
		vec2i center;
		unsigned depth;
		unsigned dx;
		unsigned corners[2][2];	// index reference to the corner indices
		unsigned index;
	} leaf2;
	//
	virtual ~octree2();
	//
	void copy( const octree2 &octree );
	void clear();
	//
	void build_octree( std::function<double(const vec2d &p)> hint, unsigned maxdepth );
	unsigned hit_test( vec2d p, leaf2 *leaf=nullptr ) const;
	void draw_octree( graphics_engine &g ) const;
	//
	leaf2 *m_root {nullptr};
	std::vector<leaf2 *> m_terminals;	// Octree terminal list
	std::vector<vec2d> m_nodes;			// Octree corner list
	//
	unsigned m_maxdepth;
	unsigned m_resolution;
	void subdivide( leaf2 *leaf, std::function<double(vec2d p)> hint, unsigned maxdepth );
	//
private:
	//
	void count_num_terminal( leaf2 *leaf, unsigned &index );
	void build_array( leaf2 *leaf, unsigned &index );
	void build_nodes();
	bool release_children( leaf2 *leaf );
	leaf2* alloc_leaf( vec2i center, unsigned depth, vec2i position );
	void draw_octree( graphics_engine &g, const leaf2 *leaf ) const;
	void copy( leaf2 *src, leaf2 *dest );
	uint64_t compute_corner_index( vec2i p ) const;
	bool check_subdivision( vec2i pos, unsigned dx, std::function<double(vec2d p)> hint, int threshold, int depth, unsigned max_nest, unsigned maxdepth ) const;
	//
};
//
SHKZ_END_NAMESPACE
//
#endif
