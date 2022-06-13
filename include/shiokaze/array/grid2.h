/*
**	grid2.h
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on November 17, 2019.
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
#ifndef SHKZ_GRID2_H
#define SHKZ_GRID2_H
//
#include <shiokaze/array/array2.h>
#include <shiokaze/math/coordsys.h>
//
SHKZ_BEGIN_NAMESPACE
//
/** @file */
/// \~english @brief Two dimensional grid class designed to be defined as instance member in recursive_configurable class.
/// \~japanese @brief recursive_configurable インスタンスのメンバーインスタンスとして定義可能な2次元グリッドクラス。
template<class T> class grid2 : public array2<T> {
public:
	//
	using array2<T>::array2;
	using array2<T>::initialize;
	/**
	 \~english @brief Allocate grid memory with value.
	 @param[in] shape Shape of the grid.
	 @param[in] coord Coordinate system the grid.
	 @param[in] value Initial value
	 \~japanese @brief グリッドを値でメモリに展開する。
	 @param[in] shape グリッドの形。
	 @param[in] coord 座標系。
	 @param[in] value 初期値。
	 */
	void initialize( const shape2 &shape, const coordsys2 &coord, const  T value=T()) {
		array2<T>::initialize(shape,value);
		m_coordsys.copy_coord(coord);
	}
	/**
	 \~english @brief Deep copy operation for grid2.
	 @param[in] array Reference to an instance of array to copy from.
	 \~japanese @brief grid2 のディープコピー演算子。
	 @param[in] コピーする grid2 のインスタンスへの参照。
	 */
	grid2& operator=(const grid2 &grid) {
		if( this != &grid ) {
			copy(grid);
		}
		return *this;
	}
	/**
	 \~english @brief Deep copy function for grid2.
	 @param[in] grid Reference to an instance of grid to copy from.
	 \~japanese @brief grid2 のディープコピー関数。
	 @param[in] コピーする grid2 のインスタンスへの参照。
	 */
	void copy( const grid2 &grid ) {
		array2<T>::copy(grid);
		m_coordsys.copy_coord(grid.coordsys());
	}
	/**
	 \~english @brief Return if the grid is different from an input grid.
	 @param[in] grid Target grid to compare.
	 @return \c true if the grid is different from the input grid and \c false otherwise.
	 \~japanese @brief グリッドが入力されたグリッドと違うかどうか返す。
	 @param[in] grid 目標とする比べるグリッド。
	 @return もしグリッドが入力と違うグリッドなら \c true そうでなければ \c false を返す。
	 */
	bool operator!=( const grid2<T> &grid ) const {
		return ! (*this == grid);
	}
	/**
	 \~english @brief Return if the grid is same to an input grid.
	 @param[in] v Target grid to compare.
	 @return \c true if the grid is the same to the input and \c false otherwise.
	 \~japanese @brief グリッドが入力されたグリッドと同じかどうか返す。
	 @param[in] v 目標とする比べるグリッド。
	 @return もしグリッドが入力と同じグリッドなら \c true そうでなければ \c false を返す。
	 */
	bool operator==(const grid2<T> &grid) const {
		return array2<T>::operator==(grid) && m_coordsys.equal(grid.coordsys());
	}
	/**
	 \~english @brief Get coordinate system.
	 @return Coordinate system.
	 \~japanese @brief 座標系を得る。
	 @return 座標系。
	 */
	const coordsys2& coordsys() const { return m_coordsys; }
protected:
	//
	coordsys2 m_coordsys;
};
//
SHKZ_END_NAMESPACE
//
#endif
