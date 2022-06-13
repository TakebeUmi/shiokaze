/*
**	coordsys.h
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
#ifndef SHKZ_COORDSYS_H
#define SHKZ_COORDSYS_H
//
#include <shiokaze/math/vec.h>
//
SHKZ_BEGIN_NAMESPACE
//
/** @file */
/// \~english @brief 2D Cartesian coordinate system structure.
/// \~japanese @brief 2次元カーテシアン座標系を定義する構造体。
struct coordsys2 {
	//
	/// \~english @brief Origin position.
	/// \~japanese @brief 原点。
	vec2d origin;
	//
	/// \~english @brief Reference position.
	/// \~japanese @brief 概念的な原点。
	vec2d reference;
	//
	/// \~english @brief Grid spacing.
	/// \~japanese @brief 格子幅。
	double dx;
	//
	/// \~english @brief Grid coord_type.
	/// \~japanese @brief 格子の種類。
	enum coord_type2 { CELL, NODAL, FACE };
	//
	/// \~english @brief Grid coord_type.
	/// \~japanese @brief 格子の種類。
	coord_type2 coord_type {NODAL};
	//
	/// \~english @brief Default constructor.
	/// \~japanese @brief デフォルトコンストラクタ。
	coordsys2() = default;
	/**
	 \~english @brief Constructor for coordsys2.
	 @param[in] coord Coordinate system instance.
	 \~japanese @brief coordsys2 のコンストラクタ。
	 @param[in] coord 座標系のインスタンス。
	 */
	coordsys2( const coordsys2 &coord ) {
		copy_coord(coord);
	}
	/**
	 \~english @brief Constructor for coordsys2.
	 @param[in] origin Origin position.
	 @param[in] dx Grid spacing.
	 @param[in] coord_type Coordinate coord_type.
	 \~japanese @brief coordsys2 のコンストラクタ。
	 @param[in] origin 原点。
	 @param[in] dx 格子幅。
	 @param[in] coord_type 座標系の種類。
	 */
	coordsys2( const vec2d &origin, double dx, coord_type2 coord_type=NODAL ) {
		set_coord(origin,dx,coord_type);
	}
	/**
	 \~english @brief Constructor for coordsys2.
	 @param[in] reference Reference origin position.
	 @param[in] origin Origin position.
	 @param[in] dx Grid spacing.
	 @param[in] coord_type Coordinate coord_type.
	 \~japanese @brief coordsys2 のコンストラクタ。
	 @param[in] reference リファレンス原点。
	 @param[in] origin 原点。
	 @param[in] dx 格子幅。
	 @param[in] coord_type 座標系の種類。
	 */
	coordsys2( const vec2d &reference, const vec2d &origin, double dx, coord_type2 coord_type=NODAL ) {
		set_coord(reference,origin,dx,coord_type);
	}
	/**
	 \~english @brief Copy coordinate system.
	 @param[in] coord Coordinate system instance.
	 \~japanese @brief 座標系をコピーする。
	 @param[in] coord 座標系のインスタンス。
	 */
	void copy_coord( const coordsys2 &coord ) {
		this->origin = coord.origin;
		this->reference = coord.reference;
		this->dx = coord.dx;
		this->coord_type = coord.coord_type;
	}
	/**
	 \~english @brief Set coordinate system.
	 @param[in] origin Origin position.
	 @param[in] dx Grid spacing.
	 @param[in] coord_type Coordinate coord_type.
	 \~japanese @brief 座標系を設定する。
	 @param[in] origin 原点。
	 @param[in] dx 格子幅。
	 @param[in] coord_type 座標系の種類。
	 */
	void set_coord( const vec2d &origin, double dx, coord_type2 coord_type=NODAL ) {
		this->reference = this->origin = origin;
		this->dx = dx;
		this->coord_type = coord_type;
	}
	/**
	 \~english @brief Set coordinate system.
	 @param[in] reference Reference origin position.
	 @param[in] origin Origin position.
	 @param[in] dx Grid spacing.
	 @param[in] coord_type Coordinate coord_type.
	 \~japanese @brief 座標系を設定する。
	 @param[in] reference リファレンス原点。
	 @param[in] origin 原点。
	 @param[in] dx 格子幅。
	 @param[in] coord_type 座標系の種類。
	 */
	void set_coord( const vec2d &reference, const vec2d &origin, double dx, coord_type2 coord_type=NODAL ) {
		this->origin = origin;
		this->reference = reference;
		this->dx = dx;
		this->coord_type = coord_type;
	}
	/**
	 \~english @brief Compare coordinate system.
	 @param[in] coord Coordinate system instance to compare.
	 @return Return \c true if identical and \c false otherwise.
	 \~japanese @brief 座標系を比較する。
	 @param[in] coord 比較する座標系のインスタンス。
	 @return もし同じ座標系なら \c true 違えば \c false を返す。
	 */
	const bool equal( const coordsys2 &coord ) const {
		return dx == coord.dx && origin == coord.origin && reference == coord.reference && coord_type == coord.coord_type;
	}
	/**
	 \~english @brief Check if an input coordinate systems lie on the same consistent world space.
	 @param[in] coord Coordinate system instance to compare.
	 @return Return \c true if consistent and \c false otherwise.
	 \~japanese @brief 入力の座標系が同じ一貫性のあるワールド座標系にあるか確認する。
	 @param[in] coord 比較する座標系のインスタンス。
	 @return もし同じ一貫性のある座標系なら \c true 違えば \c false を返す。
	 */
	const bool consistent( const coordsys2 &coord ) const {
		return dx == coord.dx && reference == coord.reference;
	}
	/**
	 \~english @brief Get the world coordinate from local index coordinate.
	 @param[in] i X index coordinate.
	 @param[in] j Y index coordinate.
	 @return World coordinate.
	 \~japanese @brief ローカルインデックス座標からワールド座標を得る。
	 @param[in] i X 座標系のインデックス。
	 @param[in] j Y 座標系のインデックス。
	 @return ワールド座標。
	 */
	vec2d world_coord( int i, int j ) const {
		return dx*vec2d(i,j)+origin;
	}
	/**
	 \~english @brief Get the world coordinate from local index coordinate.
	 @param[in] pi Local index coordinate.
	 @return World coordinate.
	 \~japanese @brief ローカルインデックス座標からワールド座標を得る。
	 @param[in] pi ローカル座標系のインデックス。
	 @return ワールド座標。
	 */
	vec2d world_coord( const vec2i &pi ) const {
		return dx*vec2d(pi)+origin;
	}
	/**
	 \~english @brief Get the local index coordinate from the world coordinate.
	 @param[in] x World coordinate in X dimension.
	 @param[in] y World coordinate in Y dimension.
	 @return Local index coordinate.
	 \~japanese @brief ワールド座標からローカルインデックス座標を得る。
	 @param[in] x X次元でのワールド座標。
	 @param[in] y Y次元でのワールド座標。
	 @return ローカルインデックス座標。
	 */
	vec2i local_coord( double x, double y ) const {
		return vec2i((x-origin[0])/dx-0.5,(y-origin[1])/dx-0.5);
	}
	/**
	 \~english @brief Get the local index coordinate from the world coordinate.
	 @param[in] p World coordinate.
	 @return Local index coordinate.
	 \~japanese @brief ワールド座標からローカルインデックス座標を得る。
	 @param[in] p ワールド座標。
	 @return ローカルインデックス座標。
	 */
	vec2i local_coord( const vec2d &p ) const {
		return local_coord(p[0],p[1]);
	}
	/**
	 \~english @brief Get the cell coordinate.
	 @param[in] dx Grid spacing.
	 @return Cell coordinate system.
	 \~japanese @brief セル座標系を得る。
	 @param[in] dx 格子幅。
	 @return セル格子座標系。
	 */
	static coordsys2 cell( double dx ) { return cell(vec2d(),dx); }
	/**
	 \~english @brief Get the cell coordinate.
	 @param[in] reference Reference position.
	 @param[in] dx Grid spacing.
	 @return Cell coordinate system.
	 \~japanese @brief セル座標系を得る。
	 @param[in] reference リファレンス原点。
	 @param[in] dx 格子幅。
	 @return セル格子座標系。
	 */
	static coordsys2 cell( const vec2d &reference, double dx ) {
		return coordsys2(reference,reference+dx*vec2d().cell(),dx,CELL);
	}
	/**
	 \~english @brief Get the nodal coordinate.
	 @param[in] dx Grid spacing.
	 @return Nodal coordinate system.
	 \~japanese @brief 節点系の座標系を得る。
	 @param[in] dx 格子幅。
	 @return 節点系の格子座標系。
	 */
	static coordsys2 nodal( double dx ) { return nodal(vec2d(),dx); }
	/**
	 \~english @brief Get the nodal coordinate.
	 @param[in] reference Reference position.
	 @param[in] dx Grid spacing.
	 @return Nodal coordinate system.
	 \~japanese @brief 節点系の座標系を得る。
	 @param[in] reference リファレンス原点。
	 @param[in] dx 格子幅。
	 @return 節点系の格子座標系。
	 */
	static coordsys2 nodal( const vec2d &reference, double dx ) {
		return coordsys2(reference,reference+dx*vec2d().nodal(),dx,NODAL);
	}
	/**
	 \~english @brief Get facial coordinate.
	 @param[in] dim Dimension number.
	 @param[in] dx Grid spacing.
	 @return Facial coordinate system.
	 \~japanese @brief 格子面の座標系を得る。
	 @param[in] dim 次元数。
	 @param[in] dx 格子幅。
	 @return 格子面の格子座標系。
	 */
	static coordsys2 face( int dim, double dx ) { return face(dim,vec2d(),dx); }
	/**
	 \~english @brief Get facial coordinate.
	 @param[in] dim Dimension number.
	 @param[in] reference Reference position.
	 @param[in] dx Grid spacing.
	 @return Facial coordinate system.
	 \~japanese @brief 格子面の座標系を得る。
	 @param[in] dim 次元数。
	 @param[in] reference リファレンス原点。
	 @param[in] dx 格子幅。
	 @return 格子面の格子座標系。
	 */
	static coordsys2 face( int dim, const vec2d &reference, double dx ) {
		return coordsys2(reference,reference+dx*vec2d().face(dim),dx,FACE);
	}
	/**
	 \~english @brief Check if multiples coordinate systems lie on the same consistent world space.
	 @param[in] coords A set of coordiante systems to examine.
	 @return \c true if the multiple coordinate systems lie on the same consistent world space.
	 \~japanese @brief 複数の座標系が同じ一貫性のあるワールド座標にあるか確認する。
	 @param[in] coords 確認する複数の座標系。
	 @return もし複数の座標系が一貫性のあるワールド座標系にあれば \c true を返し、そうでなければ \c false を返す。
	 */
	static bool consistent_check( std::initializer_list<coordsys2> coords ) {
		bool consistent (true);
		for( auto it0=coords.begin(); it0!=coords.end(); ++it0 ) for( auto it1=it0+1; it1!=coords.end(); ++it1 ) {
			if( ! it0->consistent(*it1) ) consistent = false;
		}
		return consistent;
	}
};
/** @file */
/// \~english @brief 3D Cartesian coordinate system structure.
/// \~japanese @brief 3次元カーテシアン座標系を定義する構造体。
struct coordsys3 {
	//
	/// \~english @brief Origin position.
	/// \~japanese @brief 原点。
	vec3d origin;
	//
	/// \~english @brief Reference position.
	/// \~japanese @brief 概念的な原点。
	vec3d reference;
	//
	/// \~english @brief Grid spacing.
	/// \~japanese @brief 格子幅。
	double dx;
	//
	/// \~english @brief Grid coord_type.
	/// \~japanese @brief 格子の種類。
	enum coord_type3 { CELL, NODAL, FACE, EDGE };
	//
	/// \~english @brief Grid coord_type.
	/// \~japanese @brief 格子の種類。
	coord_type3 coord_type {NODAL};
	//
	/// \~english @brief Default constructor.
	/// \~japanese @brief デフォルトコンストラクタ。
	coordsys3() = default;
	/**
	 \~english @brief Constructor for coordsys3.
	 @param[in] coord Coordinate system instance.
	 \~japanese @brief coordsys3 のコンストラクタ。
	 @param[in] coord 座標系のインスタンス。
	 */
	coordsys3( const coordsys3 &coord ) {
		copy_coord(coord);
	}
	/**
	 \~english @brief Constructor for coordsys3.
	 @param[in] origin Origin position.
	 @param[in] dx Grid spacing.
	 @param[in] coord_type Coordinate coord_type.
	 \~japanese @brief coordsys3 のコンストラクタ。
	 @param[in] origin 原点。
	 @param[in] dx 格子幅。
	 @param[in] coord_type 座標系の種類。
	 */
	coordsys3( const vec3d &origin, double dx, coord_type3 coord_type=NODAL ) {
		set_coord(origin,dx,coord_type);
	}
	/**
	 \~english @brief Constructor for coordsys3.
	 @param[in] reference Reference origin position.
	 @param[in] origin Origin position.
	 @param[in] dx Grid spacing.
	 @param[in] coord_type Coordinate coord_type.
	 \~japanese @brief coordsys3 のコンストラクタ。
	 @param[in] reference リファレンス原点。
	 @param[in] origin 原点。
	 @param[in] dx 格子幅。
	 @param[in] coord_type 座標系の種類。
	 */
	coordsys3( const vec3d &reference, const vec3d &origin, double dx, coord_type3 coord_type=NODAL ) {
		set_coord(reference,origin,dx,coord_type);
	}
	/**
	 \~english @brief Copy coordinate system.
	 @param[in] coord Coordinate system instance.
	 \~japanese @brief 座標系をコピーする。
	 @param[in] coord 座標系のインスタンス。
	 */
	void copy_coord( const coordsys3 &coord ) {
		this->origin = coord.origin;
		this->reference = coord.reference;
		this->dx = coord.dx;
	}
	/**
	 \~english @brief Set coordinate system.
	 @param[in] origin Origin position.
	 @param[in] dx Grid spacing.
	 @param[in] coord_type Coordinate coord_type.
	 \~japanese @brief 座標系を設定する。
	 @param[in] origin 原点。
	 @param[in] dx 格子幅。
	 @param[in] coord_type 座標系の種類。
	 */
	void set_coord( const vec3d &origin, double dx, coord_type3 coord_type=NODAL ) {
		this->reference = this->origin = origin;
		this->dx = dx;
		this->coord_type = coord_type;
	}
	/**
	 \~english @brief Set coordinate system.
	 @param[in] reference Reference origin position.
	 @param[in] origin Origin position.
	 @param[in] dx Grid spacing.
	 @param[in] coord_type Coordinate coord_type.
	 \~japanese @brief 座標系を設定する。
	 @param[in] reference リファレンス原点。
	 @param[in] origin 原点。
	 @param[in] dx 格子幅。
	 @param[in] coord_type 座標系の種類。
	 */
	void set_coord( const vec3d &reference, const vec3d &origin, double dx, coord_type3 coord_type=NODAL ) {
		this->origin = origin;
		this->reference = reference;
		this->dx = dx;
		this->coord_type = coord_type;
	}
	/**
	 \~english @brief Compare coordinate system.
	 @param[in] coord Coordinate system instance to compare.
	 @return Return \c true if identical and \c false otherwise.
	 \~japanese @brief 座標系を比較する。
	 @param[in] coord 比較する座標系のインスタンス。
	 @return もし同じ座標系なら \c true 違えば \c false を返す。
	 */
	const bool equal( const coordsys3 &coord ) const {
		return dx == coord.dx && origin == coord.origin && reference == coord.reference && coord_type == coord.coord_type;
	}
	/**
	 \~english @brief Check if an input coordinate systems lie on the same consistent world space.
	 @param[in] coord Coordinate system instance to compare.
	 @return Return \c true if consistent and \c false otherwise.
	 \~japanese @brief 入力の座標系が同じ一貫性のあるワールド座標系にあるか確認する。
	 @param[in] coord 比較する座標系のインスタンス。
	 @return もし同じ一貫性のある座標系なら \c true 違えば \c false を返す。
	 */
	const bool consistent( const coordsys3 &coord ) const {
		return dx == coord.dx && reference == coord.reference;
	}
	/**
	 \~english @brief Get the world coordinate from local index coordinate.
	 @param[in] i X index coordinate.
	 @param[in] j Y index coordinate.
	 @param[in] k Z index coordinate.
	 @return World coordinate.
	 \~japanese @brief ローカルインデックス座標からワールド座標を得る。
	 @param[in] i X 座標系のインデックス。
	 @param[in] j Y 座標系のインデックス。
	 @param[in] k Z 座標系のインデックス。
	 @return ワールド座標。
	 */
	vec3d world_coord( int i, int j, int k ) const {
		return dx*vec3d(i,j,k)+origin;
	}
	/**
	 \~english @brief Get the world coordinate from local index coordinate.
	 @param[in] pi Local index coordinate.
	 @return World coordinate.
	 \~japanese @brief ローカルインデックス座標からワールド座標を得る。
	 @param[in] pi ローカル座標系のインデックス。
	 @return ワールド座標。
	 */
	vec3d world_coord( const vec3i &pi ) const {
		return dx*vec3d(pi)+origin;
	}
	/**
	 \~english @brief Get the local index coordinate from the world coordinate.
	 @param[in] x World coordinate in X dimension.
	 @param[in] y World coordinate in Y dimension.
	 @param[in] z World coordinate in Z dimension.
	 @return Local index coordinate.
	 \~japanese @brief ワールド座標からローカルインデックス座標を得る。
	 @param[in] x X次元でのワールド座標。
	 @param[in] y Y次元でのワールド座標。
	 @param[in] z Z次元でのワールド座標。
	 @return ローカルインデックス座標。
	 */
	vec3i local_coord( double x, double y, double z ) const {
		return vec3i((x-origin[0])/dx-0.5,(y-origin[1])/dx-0.5,(z-origin[2])/dx-0.5);
	}
	/**
	 \~english @brief Get the local index coordinate from the world coordinate.
	 @param[in] p World coordinate.
	 @return Local index coordinate.
	 \~japanese @brief ワールド座標からローカルインデックス座標を得る。
	 @param[in] p ワールド座標。
	 @return ローカルインデックス座標。
	 */
	vec3i local_coord( const vec3d &p ) const {
		return local_coord(p[0],p[1],p[2]);
	}
	/**
	 \~english @brief Get the cell coordinate.
	 @param[in] dx Grid spacing.
	 @return Cell coordinate system.
	 \~japanese @brief セル座標系を得る。
	 @param[in] dx 格子幅。
	 @return セル格子座標系。
	 */
	static coordsys3 cell( double dx ) { return cell(vec3d(),dx); }
	/**
	 \~english @brief Get the cell coordinate.
	 @param[in] reference Reference position.
	 @param[in] dx Grid spacing.
	 @return Cell coordinate system.
	 \~japanese @brief セル座標系を得る。
	 @param[in] reference リファレンス原点。
	 @param[in] dx 格子幅。
	 @return セル格子座標系。
	 */
	static coordsys3 cell( const vec3d &reference, double dx ) {
		return coordsys3(reference,reference+dx*vec3d().cell(),dx,CELL);
	}
	/**
	 \~english @brief Get the nodal coordinate.
	 @param[in] dx Grid spacing.
	 @return Nodal coordinate system.
	 \~japanese @brief 節点系の座標系を得る。
	 @param[in] dx 格子幅。
	 @return 節点系の格子座標系。
	 */
	static coordsys3 nodal( double dx ) { return nodal(vec3d(),dx); }
	/**
	 \~english @brief Get the nodal coordinate.
	 @param[in] reference Reference position.
	 @param[in] dx Grid spacing.
	 @return Nodal coordinate system.
	 \~japanese @brief 節点系の座標系を得る。
	 @param[in] reference リファレンス原点。
	 @param[in] dx 格子幅。
	 @return 節点系の格子座標系。
	 */
	static coordsys3 nodal( const vec3d &reference, double dx ) {
		return coordsys3(reference,reference+dx*vec3d().nodal(),dx,NODAL);
	}
	/**
	 \~english @brief Get facial coordinate.
	 @param[in] dim Dimension number.
	 @param[in] dx Grid spacing.
	 @return Facial coordinate system.
	 \~japanese @brief 格子面の座標系を得る。
	 @param[in] dim 次元数。
	 @param[in] dx 格子幅。
	 @return 格子面の格子座標系。
	 */
	static coordsys3 face( int dim, double dx ) { return face(dim,vec3d(),dx); }
	/**
	 \~english @brief Get facial coordinate.
	 @param[in] dim Dimension number.
	 @param[in] reference Reference position.
	 @param[in] dx Grid spacing.
	 @return Facial coordinate system.
	 \~japanese @brief 格子面の座標系を得る。
	 @param[in] dim 次元数。
	 @param[in] reference リファレンス原点。
	 @param[in] dx 格子幅。
	 @return 格子面の格子座標系。
	 */
	static coordsys3 face( int dim, const vec3d &reference, double dx ) {
		return coordsys3(reference,reference+dx*vec3d().face(dim),dx,FACE);
	}
	/**
	 \~english @brief Get edge coordinate.
	 @param[in] dim Dimension number.
	 @param[in] dx Grid spacing.
	 @return Facial coordinate system.
	 \~japanese @brief エッジの座標系を得る。
	 @param[in] dim 次元数。
	 @param[in] dx 格子幅。
	 @return 格子面の格子座標系。
	 */
	static coordsys3 edge( int dim, double dx ) { return edge(dim,vec3d(),dx); }
	/**
	 \~english @brief Get edge coordinate.
	 @param[in] dim Dimension number.
	 @param[in] reference Reference position.
	 @param[in] dx Grid spacing.
	 @return Facial coordinate system.
	 \~japanese @brief エッジの座標系を得る。
	 @param[in] dim 次元数。
	 @param[in] reference リファレンス原点。
	 @param[in] dx 格子幅。
	 @return 格子面の格子座標系。
	 */
	static coordsys3 edge( int dim, const vec3d &reference, double dx ) {
		return coordsys3(reference,reference+dx*vec3d().edge(dim),dx,EDGE);
	}
	/**
	 \~english @brief Check if multiples coordinate systems lie on the same consistent world space.
	 @param[in] coords A set of coordiante systems to examine.
	 @return \c true if the multiple coordinate systems lie on the same consistent world space.
	 \~japanese @brief 複数の座標系が同じ一貫性のあるワールド座標にあるか確認する。
	 @param[in] coords 確認する複数の座標系。
	 @return もし複数の座標系が一貫性のあるワールド座標系にあれば \c true を返し、そうでなければ \c false を返す。
	 */
	static bool consistent_check( std::initializer_list<coordsys3> coords ) {
		bool consistent (true);
		for( auto it0=coords.begin(); it0!=coords.end(); ++it0 ) for( auto it1=it0+1; it1!=coords.end(); ++it1 ) {
			if( ! it0->consistent(*it1) ) consistent = false;
		}
		return consistent;
	}
};
//
SHKZ_END_NAMESPACE
//
#endif
//