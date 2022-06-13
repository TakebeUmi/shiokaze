/*
**	unstructured_extrapolator2.h
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on Februrary 6, 2018.
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
#ifndef SHKZ_UNSTRUCTURED_EXTRAPOLATOR2_H
#define SHKZ_UNSTRUCTURED_EXTRAPOLATOR2_H
//
#include <vector>
#include <algorithm>
#include <limits>
#include <numeric>
#include <shiokaze/math/vec.h>
#include <functional>
//
SHKZ_BEGIN_NAMESPACE
//
template<class T> class unstructured_extrapolator2 {
public:
	//
	static void extrapolate(
					  std::function<vec2r( size_t i )> position_func,
					  std::function<void( size_t i, std::function<void( size_t j )> func )> iterate_connections,
					  const std::vector<Real> &levelset,
					  std::vector<char> &fixed,
					  std::vector<T> &values,
					  double distance ) {
		//
		// Compute approximate volume
		std::vector<Real> volumes(levelset.size());
		for( size_t n=0; n<levelset.size(); ++n ) {
			double v (0.0);
			iterate_connections(n,[&]( size_t m ) {
				v += (position_func(n)-position_func(m)).len();
			});
			volumes[n] = v;
		}
		//
		// Sort by distance
		std::vector<size_t> order_map(levelset.size());
		std::iota(order_map.begin(), order_map.end(), 0);
		std::sort(order_map.begin(), order_map.end(), [&](size_t a, size_t b){ return std::abs(levelset[a]) < std::abs(levelset[b]); });
		//
		// Propagate
		for( size_t _n=0; _n<order_map.size(); ++_n ) {
			size_t n = order_map[_n];
			if( ! fixed[n] ) {
				T sum = T();
				double wsum (0.0);
				iterate_connections(n,[&]( size_t index ) {
					if( fixed[index] ) {
						double w = volumes[index] / (std::abs(levelset[index])+1e-8);
						sum += w*values[index];
						wsum += w;
					}
				});
				if( wsum ) {
					values[n] = sum / wsum;
					fixed[n] = true;
				}
				if( std::abs(levelset[n]) > distance ) break;
			}
		}
	}
};
//
SHKZ_END_NAMESPACE
//
#endif