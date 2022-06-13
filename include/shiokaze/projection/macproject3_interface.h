/*
**	macproject3_interface.h
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on April 26, 2017.
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
#ifndef SHKZ_MACPROJECT3_INTERFACE_H
#define SHKZ_MACPROJECT3_INTERFACE_H
//
#include <shiokaze/graphics/graphics_engine.h>
#include <shiokaze/core/recursive_configurable_module.h>
#include <shiokaze/array/macarray3.h>
#include "signed_rigidbody3_interface.h"
#include <functional>
//
SHKZ_BEGIN_NAMESPACE
//
/** @file */
/// \~english @brief Interface that makes compressible vector field incompressible.
/// \~japanese @brief 圧縮性のあるベクトル場を非圧縮場なベクトル場に変換するインターフェース。
class macproject3_interface : public recursive_configurable_module {
public:
	//
	DEFINE_MODULE(macproject3_interface,"MAC Project 3D","Projection","Projection module")
	/**
	 \~english @brief Set a target volume.
	 @param[in] current_volume Current volume.
	 @param[in] target_volume Target volume.
	 \~japanese @brief 目標の体積をセットする。
	 @param[in] current_volume 現在の体積。
	 @param[in] target_volume 目標とする体積。
	 */
	virtual void set_target_volume( double current_volume, double target_volume ) {}
	/**
	 \~english @brief Set solid density, allowing some fluid to pass through solid faces.
	 @param[in] solid_density Solid face density array.
	 \~japanese @brief 壁面の密度を設定する。流体の一部が壁を突き抜けられるようにする。
	 @param[in] solid_density 密度の格子。
	 */
	virtual void set_solid_density( const macarray3<Real> *solid_density ) {}
	/**
	 \~english @brief Project a vector field onto an imcompressible vector field.
	 @param[in] dt Time step size.
	 @param[in-out] velocity Velocity field.
	 @param[in] solid Solid level set.
	 @param[in] solid_velocity Solid face velocity.
	 @param[in] fluid Fluid level set.
	 @param[in] surface_tension Surface tension coefficient.
	 @param[in] moving_solid_func Moving solid function.
	 @param[in] rigidbodies Rigidbody list.
	 \~japanese @brief ベクトル場を非圧縮に投影する。
	 @param[in] dt タイムステップサイズ。
	 @param[in-out] velocity ベクトル場。
	 @param[in] solid 壁のレベルセット。
	 @param[in] solid_velocity 壁の面の速度。
	 @param[in] fluid 水のレベルセット。
	 @param[in] surface_tension 表面張力の係数。
	 @param[in] moving_solid_func 動的に動く壁の関数。
	 @param[in] rigidbodies 剛体リスト。
	 */
	virtual void project (
			double dt, macarray3<Real> &velocity,
			const array3<Real> &solid,
			const macarray3<Real> &solid_velocity,
			const array3<Real> &fluid,
			double surface_tension=0.0,
			const std::vector<signed_rigidbody3_interface *> *rigidbodies=nullptr
	) = 0;
	/**
	 \~english @brief Get pressure field if available.
	 @return Pressure field.
	 \~japanese @brief 圧力場を得る。
	 @return 圧力場。
	*/
	virtual const array3<Real>* get_pressure() const = 0;
	/**
	 \~english @brief Draw internal information.
	 @param[in] g Graphics engine.
	 \~japanese @brief 内部情報を描画する。
	 @param[in] g gグラフィックスエンジン。
	 */
	virtual void draw( graphics_engine &g ) const {}
	//
private:
	//
	virtual void initialize( const shape3 &shape, double dx ) = 0;
	virtual void initialize( const configurable::environment_map &environment ) override {
		//
		assert(check_set(environment,{"shape","dx"}));
		initialize(
			get_env<shape3>(environment,"shape"),
			get_env<double>(environment,"dx")
		);
	}
};
//
using macproject3_ptr = std::unique_ptr<macproject3_interface>;
using macproject3_driver = recursive_configurable_driver<macproject3_interface>;
//
SHKZ_END_NAMESPACE
//
#endif