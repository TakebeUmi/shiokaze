/*
**	ann3.h
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on Feb 15, 2017.
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
#include <cstdint>
#include <ann/ANN.h>
#include <shiokaze/math/vec.h>
#include <shiokaze/core/console.h>
//
#ifndef SHKZ_ANN3_H
#define SHKZ_ANN3_H
//
/** @file */
/// \~english @brief Class that searches k-nearest neighbors of points.
/// \~japanese @brief ポイント群の近傍探索を行うクラス。
SHKZ_BEGIN_NAMESPACE
//
class ann3 {
public:
	/**
	 \~english @brief Constructor for ann2.
	 \~japanese @brief ann2 のコンストラクタ。
	 */
	ann3() {
		numbers = 0;
		kdtree = nullptr;
	}
	/**
	 \~english @brief Destructor for ann2.
	 \~japanese @brief ann2 のデストラクタ。
	 */
	virtual ~ann3() {
		clear();
	}
	/**
	 \~english @brief Clear out internal data.
	 \~japanese @brief 内部情報をクリアする。
	 */
	void clear() {
		if( numbers ) {
			delete kdtree;
			delete [] dataPts[0];
			delete [] dataPts;
		}
		numbers = 0;
	}
	/**
	 \~english @brief Sort and construct a kdtree for this set of points.
	 @param[in] array Input array of points.
	 \~japanese @brief ポイント群をソートして kdtree を作成する。
	 @param[in] array ポイント群の配列の入力。
	 */
	void sort( const std::vector<vec3d> &array ) {
		clear();
		numbers = (size_t) array.size();
		if( ! numbers ) {
			console::dump( "Trying to sort empty array !\n");
			exit(0);
			return;
		}
		// Create data point set
		dataPts = annAllocPts(numbers,3);
		// Assign point information
		for( size_t n=0; n<numbers; n++ ) {
			for( size_t dim=0; dim<3; dim++) dataPts[n][dim] = array[n][dim];
		}
		// Build KD Tree
		kdtree = new ANNkd_tree(dataPts,numbers,3);
	}
	/**
	 \~english @brief Get a set of nearst points from an input position.
	 @param[in] p Input position.
	 @param[in] n Number of query points.
	 \~japanese @brief 入力の位置の近傍の粒子の一覧を取得する。
	 @param[in] p 入力位置。
	 @param[in] n 取得するポイントの数。
	 */
	std::vector<size_t> get_neighbors( const vec3d &p, size_t n ) const {
		std::vector<size_t> res(n);
		if( ! numbers ) {
			console::dump( "ANN not initialized !\n");
			exit(0);
			return std::vector<size_t>();
		}
		ANNpoint queryPt = annAllocPt(3);
		for( size_t dim=0; dim<3; dim++ ) queryPt[dim] = p[dim];
		ANNidx *nnIdx = new ANNidx[n];
		ANNdist *dists = new ANNdist[n];
		kdtree->annkSearch(queryPt,n,nnIdx,dists,0.0);
		for( size_t k=0; k<n; k++ ) {
			res[k] = nnIdx[k];
		}
		delete [] queryPt;
		delete [] nnIdx;
		delete [] dists;
		return std::move(res);
	}
private:
	size_t numbers;
	ANNpointArray dataPts;
	ANNkd_tree * kdtree;
};
//
SHKZ_END_NAMESPACE
//
#endif
//