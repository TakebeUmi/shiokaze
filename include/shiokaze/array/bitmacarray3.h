/*
**	bitmacarray3.h
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on Ausugst 8, 2018.
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
#ifndef SHKZ_BITMACARRAY3_H
#define SHKZ_BITMACARRAY3_H
//
#include <shiokaze/array/bitarray3.h>
#include <algorithm>
#include <array>
//
SHKZ_BEGIN_NAMESPACE
//
/** @file */
/// \~english @brief Three dimensional staggered bit grid class designed to be defined as instance member in recursive_configurable class.
/// \~japanese @brief recursive_configurable インスタンスのメンバーインスタンスとして定義可能なスタッガード3次元ビット配列クラス。
class bitmacarray3 : public recursive_configurable, public messageable {
public:
	/**
	 \~english @brief Constructor for bitmacarray3.
	 @param[in] parent Pointer to a parent recursive_configurable instance. Can be nullptr.
	 @param[in] shape Shape of the grid
	 @param[in] core_name Core module name. Default value is "lineararray_core3". Can be also "tiledarray_core3".
	 \~japanese @brief bitmacarray3 のコンストラクタ。
	 @param[in] parent 親 recursive_configurable のインスタンスへのポインタ。nullptr も可。
	 @param[in] shape グリッドの形
	 @param[in] core_name コア子ジュールの名前。デフォルトは "lineararray_core3"。"tiledarray_core3" も可能。
	 */
	bitmacarray3 ( recursive_configurable *parent, const shape3 &shape, std::string core_name="") : m_shape(shape), m_array_0(this,shape.face(0),core_name), m_array_1(this,shape.face(1),core_name), m_array_2(this,shape.face(2),core_name) {
		if( parent ) parent->add_child(this);
		else setup_now();
	}
	/**
	 \~english @brief Constructor for bitmacarray3.
	 @param[in] parent Pointer to a parent recursive_configurable instance. Can be nullptr.
	 @param[in] core_name Core module name. Default value is "lineararray_core3". Can be also "tiledarray_core3".
	 \~japanese @brief bitmacarray3 のコンストラクタ。
	 @param[in] parent 親 recursive_configurable のインスタンスへのポインタ。nullptr も可。
	 @param[in] core_name コア子ジュールの名前。デフォルトは "lineararray_core3"。"tiledarray_core3" も可能。
	 */
	bitmacarray3( recursive_configurable *parent, std::string core_name="" ) : bitmacarray3(parent,shape3(0,0,0),core_name) {}
	/**
	 \~english @brief Constructor for bitmacarray3.
	 @param[in] core_name Core module name. Default value is "lineararray_core3". Can be also "tiledarray_core3".
	 \~japanese @brief bitmacarray3 のコンストラクタ。
	 @param[in] core_name コア子ジュールの名前。デフォルトは "lineararray_core3"。"tiledarray_core3" も可能。
	 */
	bitmacarray3 ( std::string core_name="") : bitmacarray3(nullptr,shape3(0,0,0),core_name) {}
	/**
	 \~english @brief Constructor for bitmacarray3.
	 @param[in] shape Shape of the grid
	 @param[in] core_name Core module name. Default value is "lineararray_core3". Can be also "tiledarray_core3".
	 \~japanese @brief bitmacarray3 のコンストラクタ。
	 @param[in] shape グリッドの形
	 @param[in] core_name コア子ジュールの名前。デフォルトは "lineararray_core3"。"tiledarray_core3" も可能。
	 */
	bitmacarray3( const shape3 &shape, std::string core_name="") : bitmacarray3(nullptr,shape,core_name) {}
	/**
	 \~english @brief Copy constructor for bitmacarray3.
	 @param[in] array Reference to an instance of array to copy from.
	 \~japanese @brief bitmacarray3 のコピーコンストラクタ。
	 @param[in] コピーする bitmacarray3 のインスタンスへの参照。
	 */
	bitmacarray3 ( const bitmacarray3& v ) : m_array_0(this), m_array_1(this), m_array_2(this){
		copy(v);
	}
	/**
	 \~english @brief Send a message to the core module.
	 @param[in] message Message
	 @param[in] ptr Pointer to some value.
	 @return \c true if handled \c false otherwise.
	 \~japanese @brief コアモジュールにメッセージを送る
	 @param[in] message メッセージ
	 @param[in] ptr あるポインターの値
	 @return もし処理されたら \c true を、処理されなかったら \c false
	 */
	virtual bool send_message( std::string message, void *ptr ) override {
		bool handled (false);
		if( m_array_0.send_message(message,ptr)) handled = true;
		if( m_array_1.send_message(message,ptr)) handled = true;
		if( m_array_2.send_message(message,ptr)) handled = true;
		return handled;
	}
	/**
	 \~english @brief Send a message to the core module.
	 @param[in] message Message
	 @param[in] ptr Pointer to some value.
	 @return \c true if handled \c false otherwise.
	 \~japanese @brief コアモジュールにメッセージを送る
	 @param[in] message メッセージ
	 @param[in] ptr あるポインターの値
	 @return もし処理されたら \c true を、処理されなかったら \c false
	 */
	virtual bool const_send_message( std::string message, void *ptr ) const override {
		bool handled (false);
		if( m_array_0.const_send_message(message,ptr)) handled = true;
		if( m_array_1.const_send_message(message,ptr)) handled = true;
		if( m_array_2.const_send_message(message,ptr)) handled = true;
		return handled;
	}
	/**
	 \~english @brief Deep copy operation for bitmacarray3.
	 @param[in] array Reference to an instance of macarray to copy from.
	 \~japanese @brief bitmacarray3 のディープコピー演算子。
	 @param[in] コピーする bitmacarray3 のインスタンスへの参照。
	 */
	bitmacarray3& operator=(const bitmacarray3 &array) {
		copy(array);
		return *this;
	}
	/**
	 \~english @brief Deep copy function for bitmacarray3.
	 @param[in] array Reference to an instance of array to copy from.
	 \~japanese @brief bitmacarray3 のディープコピー関数。
	 @param[in] コピーする bitmacarray3 のインスタンスへの参照。
	 */
	void copy( const bitmacarray3 &array ) {
		if( this != &array ) {
			set_type(array.type());
			for( int dim : DIMS3 ) (*this)[dim].copy(array[dim]);
		}
	}
	/**
	 \~english @brief Allocate grid memory with value.
	 @param[in] shape Shape of the grid.
	 @param[in] value Initial value
	 \~japanese @brief グリッドを値でメモリに展開する。
	 @param[in] shape グリッドの形。
	 @param[in] value 初期値。
	 */
	void initialize ( const shape3 &shape ) {
		m_shape = shape;
		for( int dim : DIMS3 ) (*this)[dim].initialize(shape.face(dim));
	}
	/**
	 \~english @brief Function to count the number of active cells.
	 @return Active cell count.
	 \~japanese @brief アクティブセルの数を数える関数。
	 @return アクティブセルの数。
	 */
	size_t count () const {
		size_t sum (0);
		for( int dim : DIMS3 ) sum += (*this)[dim].count();
		return sum;
	}
	/**
	 \~english @brief Function to return the list of active cells positions.
	 @return The list of active cells positions.
	 \~japanese @brief アクティブセルの位置のリストを返す関数。
	 @return アクティブセルの位置のリスト。
	 */
	std::array<std::vector<vec3i>,DIM3> actives() const {
		std::array<std::vector<vec3i>,DIM3> result;
		m_parallel.for_each( DIM3, [&]( size_t dim ) {
			result[dim] = (*this)[dim].actives();
		});
		return result;
	}
	/**
	 \~english @brief Activate cells at the positons of active_entries.
	 @param[in] active_entries The list of target positions to activate.
	 @param[in] offset Offset applied to the active_entries.
	 \~japanese @brief active_entries と同じ場所のセルを offset だけずらして、アクティブにする。
	 @param[in] active_entries アクティブにするセルの場所のリスト。
	 @param[in] offset active_entries に適用されるオフセット。
	 */
	void activate( const std::array<std::vector<vec3i>,DIM3> &active_entries, const std::array<vec3i,DIM3> &offsets={vec3i(),vec3i(),vec3i()}  ) {
		m_parallel.for_each( DIM3, [&]( size_t dim ) {
			(*this)[dim].activate(active_entries[dim],offsets[dim]);
		});
	}
	/**
	 \~english @brief Activate cells at the same positons where an input array is active with an offset.
	 @param[in] array Target array.
	 @param[in] offset Offset applied to the target array.
	 \~japanese @brief 入力のグリッドのアクティブセルと同じ場所のセルを offset だけずらして、アクティブにする。
	 @param[in] array 目標となるグリッド。
	 @param[in] offset 目標となるグリッドに適用されるオフセット。
	 */
	void activate_as_bit( const bitmacarray3 &array, const std::array<vec3i,DIM3> &offsets={vec3i(),vec3i(),vec3i()} ) {
		m_parallel.for_each( DIM3, [&]( size_t dim ) {
			(*this)[dim].activate_as_bit(array[dim],offsets[dim]);
		});
	}
	/**
	 \~english @brief Activate cells at the same positons where an input array is active with an offset.
	 @param[in] array Target array.
	 @param[in] offset Offset applied to the target array.
	 \~japanese @brief 入力のグリッドのアクティブセルと同じ場所のセルを offset だけずらして、アクティブにする。
	 @param[in] array 目標となるグリッド。
	 @param[in] offset 目標となるグリッドに適用されるオフセット。
	 */
	template <class Y> void activate_as( const Y &array, const std::array<vec3i,DIM3> &offsets={vec3i(),vec3i(),vec3i()} ) {
		m_parallel.for_each( DIM3, [&]( size_t dim ) {
			(*this)[dim].activate_as(array[dim],offsets[dim]);
		});
	}
	/**
	 \~english @brief Activate all the cells.
	 \~japanese @brief 全てのセルをアクティブにする。
	 */
	void activate_all() {
		m_parallel.for_each( DIM3, [&]( size_t dim ) {
			(*this)[dim].activate_all();
		});
	}
	/**
	 \~english @brief Copy the states of active and inactive cells as same as input array with an offset.
	 @param[in] array Target input array from which the states to be copied.
	 @param[in] offset Offset
	 \~japanese @brief セルのアクティブと非アクティブステートの状態を入力のグリッドと同じようにセットする。
	 @param[in] array 目標となる状態をコピーする元となるグリッド。
	 @param[in] offset オフセット
	 */
	void copy_active_as( const bitmacarray3 &array, const vec3i &offset=vec3i() ) {
		m_parallel.for_each( DIM3, [&]( size_t dim ) {
			(*this)[dim].copy_active_as(array[dim],offset);
		});
	}
	/**
	 \~english @brief Get the shape of the array.
	 @return Shape of the array.
	 \~japanese @brief グリッドの形を返す。
	 @return グリッドの形。
	 */
	shape3 shape() const { return m_shape; }
	/**
	 \~english @brief Get the shape of the staggered grid of a specified dimension.
	 @return Shape of the grid of an input dimension.
	 \~japanese @brief 指定された次元のスタッガードグリッドの形を返す。
	 @return 指定された次元でのグリッドの形。
	 */
	shape3 shape(int dim) const { return (*this)[dim].shape(); }
	/**
	 \~english @brief Clear out the grid.
	 *
	 Note that size, the memory allocation, background values and the information regarding level set or fillable left intact.
	 \~japanese @brief グリッドを初期化する。
	 *
	 グリッドの大きさやメモリ確保、レベルセット情報、バックグラウンド値は変更されない。
	 */
	void clear() {
		for( int dim : DIMS3 ) (*this)[dim].clear();
	}
	/**
	 \~english @brief Return if the grid is different from an input array.
	 @param[in] array Target array to compare.
	 @return \c true if the array is different from the input array and \c false otherwise.
	 \~japanese @brief グリッドが入力されたグリッドと違うかどうか返す。
	 @param[in] array 目標とする比べるグリッド。
	 @return もしグリッドが入力と違うグリッドなら \c true そうでなければ \c false を返す。
	 */
	bool operator!=( const bitmacarray3 &v ) const {
		return ! (*this == v);
	}
	/**
	 \~english @brief Return if the grid is same to an input array.
	 @param[in] array Target array to compare.
	 @return \c true if the array is the same to the input and \c false otherwise.
	 \~japanese @brief グリッドが入力されたグリッドと同じかどうか返す。
	 @param[in] array 目標とする比べるグリッド。
	 @return もしグリッドが入力と同じグリッドなら \c true そうでなければ \c false を返す。
	 */
	bool operator==(const bitmacarray3 &v) const {
		for( int dim : DIMS3 ) {
			if ((*this)[dim] != v[dim]) return false;
		}
		return true;
	}
	/**
	 \~english @brief Get the read-only reference to the staggered array of a specified dimension.
	 @param[in] dim Dimensiton of the grid.
	 \~japanese @brief 指定した次元の読み込み処理だけ可能なスタッガード格子の参照を得る。
	 @param[in] dim グリッドの次元。
	 */
	const bitarray3& operator[](int dim) const {
		return dim==0 ? m_array_0 : (dim == 1 ? m_array_1 : m_array_2 );
	}
	/**
	 \~english @brief Get the reference to the staggered array of a specified dimension.
	 @param[in] dim Dimensiton of the grid.
	 \~japanese @brief 指定した次元のスタッガード格子の参照を得る。
	 @param[in] dim グリッドの次元。
	 */
	bitarray3& operator[](int dim) {
		return dim==0 ? m_array_0 : (dim == 1 ? m_array_1 : m_array_2 );
	}
	/**
	 \~english @brief Set the number of threads for parallel processing on this grid.
	 @param[in] number Number of threads.
	 \~japanese @brief 並列処理をするためのスレッドの数を設定する。
	 @param[in] number スレッドの数。
	 */
	void set_thread_num( int number ) {
		for( int dim : DIMS3 ) (*this)[dim].set_thread_num(number);
	}
	/**
	 \~english @brief Get the current number of threads for parallel processing on this grid.
	 @return number Number of threads.
	 \~japanese @brief 現在設定されている並列処理をするためのスレッドの数を得る。
	 @return number スレッドの数。
	 */
	int get_thread_num() const {
		return m_array_0.get_thread_num();
	}
	//
	enum { ACTIVES = true, ALL = false };
	/**
	 \~english @brief Loop over all the active cells in parallel.
	 @param[in] func Function that defines how each grid cell is processed.
	 \~japanese @brief アクティブセルを並列に処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 */
	inline void parallel_actives( std::function<void(typename bitarray3::iterator& it)> func ) { parallel_op(func,ACTIVES); }
	/**
	 \~english @brief Loop over all the cells in parallel.
	 @param[in] func Function that defines how each grid cell is processed.
	 \~japanese @brief 全てのセルを並列に処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 */
	inline void parallel_all( std::function<void(typename bitarray3::iterator& it)> func ) { parallel_op(func,ALL); }
	/**
	 \~english @brief Loop over cells in parallel.
	 @param[in] func Function that defines how each grid cell is processed.
	 @param[in] type Type of target cells. ACTIVE or ALL.
	 \~japanese @brief セルを並列に処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 @param[in] type ターゲットセルのタイプ. ACTIVE か ALL。
	 */
	inline void parallel_op( std::function<void(typename bitarray3::iterator& it)> func, bool type=ALL ) {
		parallel_op([func](int dim, int i, int j, int k, typename bitarray3::iterator& it, int thread_index){
			func(it);
		},type);
	}
	/**
	 \~english @brief Loop over all the active cells in parallel.
	 @param[in] func Function that defines how each grid cell is processed.
	 \~japanese @brief アクティブセルを並列に処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 */
	inline void parallel_actives( std::function<void(int dim, int i, int j, int k, typename bitarray3::iterator& it)> func ) { parallel_op(func,ACTIVES); }
	/**
	 \~english @brief Loop over all the cells in parallel.
	 @param[in] func Function that defines how each grid cell is processed.
	 \~japanese @brief 全てのセルを並列に処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 */
	inline void parallel_all( std::function<void(int dim, int i, int j, int k, typename bitarray3::iterator& it)> func ) { parallel_op(func,ALL); }
	/**
	 \~english @brief Loop over cells in parallel.
	 @param[in] func Function that defines how each grid cell is processed.
	 @param[in] type Type of target cells. ACTIVE or ALL.
	 \~japanese @brief セルを並列に処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 @param[in] type ターゲットセルのタイプ. ACTIVE か ALL。
	 */
	inline void parallel_op( std::function<void(int dim, int i, int j, int k, typename bitarray3::iterator& it)> func, bool type=ALL ) {
		parallel_op([func](int dim, int i, int j, int k, typename bitarray3::iterator& it, int thread_index){
			func(dim,i,j,k,it);
		},type);
	}
	/**
	 \~english @brief Loop over all the active cells in parallel.
	 @param[in] func Function that defines how each grid cell is processed.
	 \~japanese @brief アクティブセルを並列に処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 */
	inline void parallel_actives( std::function<void(int dim, int i, int j, int k, typename bitarray3::iterator& it, int thread_index)> func ) { parallel_op(func,ACTIVES); }
	/**
	 \~english @brief Loop over all the cells in parallel.
	 @param[in] func Function that defines how each grid cell is processed.
	 \~japanese @brief 全てのセルを並列に処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 */
	inline void parallel_all( std::function<void(int dim, int i, int j, int k, typename bitarray3::iterator& it, int thread_index)> func ) { parallel_op(func,ALL); }
	/**
	 \~english @brief Loop over cells in parallel.
	 @param[in] func Function that defines how each grid cell is processed.
	 @param[in] type Type of target cells. ACTIVE or ALL.
	 \~japanese @brief セルを並列に処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 @param[in] type ターゲットセルのタイプ. ACTIVE か ALL。
	 */
	inline void parallel_op( std::function<void(int dim, int i, int j, int k, typename bitarray3::iterator& it, int thread_index)> func, bool type=ALL ) {
		for( int dim : DIMS3 ) {
			(*this)[dim].parallel_op([&](int i, int j, int k, typename bitarray3::iterator& it, int thread_index) {
				func(dim,i,j,k,it,thread_index);
			},type);
		};
	}
	/**
	 \~english @brief Loop over all the cells in parallel by read-only fashion.
	 @param[in] func Function that defines how each grid cell is processed.
	 \~japanese @brief 全てのセルを並列に読み込みのみで処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 */
	inline void const_parallel_all( std::function<void(const typename bitarray3::const_iterator& it)> func ) const { const_parallel_op(func,ALL); }
	/**
	 \~english @brief Loop over cells in parallel by read-only fashion.
	 @param[in] func Function that defines how each grid cell is processed.
	 @param[in] type Type of target cells. ACTIVE or ALL.
	 \~japanese @brief セルを並列に読み込みのみで処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 @param[in] type ターゲットセルのタイプ. ACTIVE か ALL。
	 */
	inline void const_parallel_op( std::function<void(const typename bitarray3::const_iterator& it)> func, bool type=ALL ) const {
		const_parallel_op([func](int dim, int i, int j, int k, const typename bitarray3::const_iterator& it, int thread_index){
			func(it);
		},type);
	}
	/**
	 \~english @brief Loop over all the active cells in parallel by read-only fashion.
	 @param[in] func Function that defines how each grid cell is processed.
	 \~japanese @brief アクティブセルを並列に読み込みのみで処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 */
	inline void const_parallel_actives( std::function<void(int dim, int i, int j, int k)> func ) const {
		const_parallel_op([&](int dim, int i, int j, int k, const typename bitarray3::const_iterator& it){ func(dim,i,j,k); }, ACTIVES); }
	/**
	 \~english @brief Loop over all the cells in parallel by read-only fashion.
	 @param[in] func Function that defines how each grid cell is processed.
	 \~japanese @brief 全てのセルを並列に読み込みのみで処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 */
	inline void const_parallel_all( std::function<void(int dim, int i, int j, int k, const typename bitarray3::const_iterator& it)> func ) const { const_parallel_op(func,ALL); }
	/**
	 \~english @brief Loop over cells in parallel by read-only fashion.
	 @param[in] func Function that defines how each grid cell is processed.
	 @param[in] type Type of target cells. ACTIVE or ALL.
	 \~japanese @brief セルを並列に読み込みのみで処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 @param[in] type ターゲットセルのタイプ. ACTIVE か ALL。
	 */
	inline void const_parallel_op( std::function<void(int dim, int i, int j, int k, const typename bitarray3::const_iterator& it)> func, bool type=ALL ) const {
		const_parallel_op([func](int dim, int i, int j, int k, const typename bitarray3::const_iterator& it, int thread_index){
			func(dim,i,j,k,it);
		},type);
	}
	/**
	 \~english @brief Loop over all the active cells in parallel by read-only fashion.
	 @param[in] func Function that defines how each grid cell is processed.
	 \~japanese @brief アクティブセルを並列に読み込みのみで処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 */
	inline void const_parallel_actives( std::function<void(int dim, int i, int j, int k, int thread_index)> func ) const {
		const_parallel_op([&](int dim, int i, int j, int k, const typename bitarray3::const_iterator& it, int thread_index){ func(dim,i,j,k,thread_index); }, ACTIVES); }
	/**
	 \~english @brief Loop over all the cells in parallel by read-only fashion.
	 @param[in] func Function that defines how each grid cell is processed.
	 \~japanese @brief 全てのセルを並列に読み込みのみで処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 */
	inline void const_parallel_all( std::function<void(int dim, int i, int j, int k, const typename bitarray3::const_iterator& it, int thread_index)> func ) const { const_parallel_op(func,ALL); }
	/**
	 \~english @brief Loop over cells in parallel by read-only fashion.
	 @param[in] func Function that defines how each grid cell is processed.
	 @param[in] type Type of target cells. ACTIVE or ALL.
	 \~japanese @brief セルを並列に読み込みのみで処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 @param[in] type ターゲットセルのタイプ. ACTIVE か ALL。
	 */
	inline void const_parallel_op( std::function<void(int dim, int i, int j, int k, const typename bitarray3::const_iterator& it, int thread_index)> func, bool type=ALL ) const {
		for( int dim : DIMS3 ) {
			(*this)[dim].const_parallel_op([&](int i, int j, int k, const typename bitarray3::const_iterator& it, int thread_index) {
				func(dim,i,j,k,it,thread_index);
			},type);
		};
	}
	/**
	 \~english @brief Loop over all the active cells in serial order.
	 @param[in] func Function that defines how each grid cell is processed.
	 \~japanese @brief アクティブなセルをシリアルに処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 */
	inline void serial_actives( std::function<void(typename bitarray3::iterator& it)> func ) { serial_op(func,ACTIVES); }
	/**
	 \~english @brief Loop over all the cells in serial order.
	 @param[in] func Function that defines how each grid cell is processed.
	 \~japanese @brief 全てのセルをシリアルに処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 */
	inline void serial_all( std::function<void(typename bitarray3::iterator& it)> func ) { serial_op(func,ALL); }
	/**
	 \~english @brief Loop over cells in serial order.
	 @param[in] func Function that defines how each grid cell is processed.
	 @param[in] type Type of target cells. ACTIVE or ALL.
	 \~japanese @brief セルをシリアルに処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 @param[in] type ターゲットセルのタイプ. ACTIVE か ALL。
	 */
	inline void serial_op( std::function<void(typename bitarray3::iterator& it)> func, bool type=ALL ) {
		serial_op([func](int dim, int i, int j, int k, typename bitarray3::iterator& it) {
			func(it);
		},type);
	}
	/**
	 \~english @brief Loop over all the active cells in serial order.
	 @param[in] func Function that defines how each grid cell is processed.
	 \~japanese @brief アクティブなセルをシリアルに処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 */
	inline void serial_actives( std::function<void(int dim, int i, int j, int k, typename bitarray3::iterator& it)> func ) { serial_op(func,ACTIVES); }
	/**
	 \~english @brief Loop over all the cells in serial order.
	 @param[in] func Function that defines how each grid cell is processed.
	 \~japanese @brief 全てのセルをシリアルに処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 */
	inline void serial_all( std::function<void(int dim, int i, int j, int k, typename bitarray3::iterator& it)> func ) { serial_op(func,ALL); }
	/**
	 \~english @brief Loop over cells in serial order.
	 @param[in] func Function that defines how each grid cell is processed.
	 @param[in] type Type of target cells. ACTIVE or ALL.
	 \~japanese @brief セルをシリアルに処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 @param[in] type ターゲットセルのタイプ. ACTIVE か ALL。
	 */
	inline void serial_op( std::function<void(int dim, int i, int j, int k, typename bitarray3::iterator& it)> func, bool type=ALL ) {
		for( int dim : DIMS3 ) {
			(*this)[dim].serial_op([&](int i, int j, int k, typename bitarray3::iterator& it) {
				func(dim,i,j,k,it);
			},type );
		}
	}
	/**
	 \~english @brief Loop over all the cells in serial order by read-only fashion.
	 @param[in] func Function that defines how each grid cell is processed.
	 \~japanese @brief 全てのセルをシリアルに読み込みのみで処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 */
	inline void const_serial_all( std::function<void(const typename bitarray3::const_iterator& it)> func ) const { const_serial_op(func,ALL); }
	/**
	 \~english @brief Loop over the cells in serial order by read-only fashion.
	 @param[in] func Function that defines how each grid cell is processed.
	 @param[in] type Type of target cells. ACTIVE or ALL.
	 \~japanese @brief セルをシリアルに読み込みのみで処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 @param[in] type ターゲットセルのタイプ. ACTIVE か ALL。
	 */
	inline void const_serial_op( std::function<void(const typename bitarray3::const_iterator& it)> func, bool type=ALL ) const {
		const_serial_op([func](int dim, int i, int j, int k, const typename bitarray3::const_iterator& it) {
			func(it);
		},type);
	}
	/**
	 \~english @brief Loop over all the active cells in serial order by read-only fashion.
	 @param[in] func Function that defines how each grid cell is processed.
	 \~japanese @brief アクティブなセルをシリアルに読み込みのみで処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 */
	inline void const_serial_actives( std::function<void(int dim, int i, int j, int k)> func ) const {
		const_serial_op([&](int dim, int i, int j, int k, const typename bitarray3::const_iterator& it){ func(dim,i,j,k); }, ACTIVES); }
	/**
	 \~english @brief Loop over all the cells in serial order by read-only fashion.
	 @param[in] func Function that defines how each grid cell is processed.
	 \~japanese @brief 全てのセルをシリアルに読み込みのみで処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 */
	inline void const_serial_all( std::function<void(int dim, int i, int j, int k, const typename bitarray3::const_iterator& it)> func ) const { const_serial_op(func,ALL); }
	/**
	 \~english @brief Loop over cells in serial order by read-only fashion.
	 @param[in] func Function that defines how each grid cell is processed.
	 @param[in] type Type of target cells. ACTIVE or ALL.
	 \~japanese @brief セルをシリアルに読み込みのみで処理する。
	 @param[in] func それぞれのセルを処理する関数。
	 @param[in] type ターゲットセルのタイプ. ACTIVE か ALL。
	 */
	inline void const_serial_op( std::function<void(int dim, int i, int j, int k, const typename bitarray3::const_iterator& it)> func, bool type=ALL ) const {
		for( int dim : DIMS3 ) {
			(*this)[dim].const_serial_op([&](int i, int j, int k, const typename bitarray3::const_iterator& it) {
				func(dim,i,j,k,it);
			},type);
		}
	}
	/**
	 \~english @brief Loop over all the active cells in serial order.
	 @param[in] func Function that defines how each grid cell is processed. Stop the loop if return true.
	 \~japanese @brief アクティブなセルをシリアルに処理する。
	 @param[in] func それぞれのセルを処理する関数。\c true を返すと、ループを中断する。
	 */
	inline void interruptible_serial_actives( std::function<bool(typename bitarray3::iterator& it)> func ) { serial_op(func,ACTIVES); }
	/**
	 \~english @brief Loop over all the cells in serial order.
	 @param[in] func Function that defines how each grid cell is processed. Stop the loop if return true.
	 \~japanese @brief 全てのセルをシリアルに処理する。
	 @param[in] func それぞれのセルを処理する関数。\c true を返すと、ループを中断する。
	 */
	inline void interruptible_serial_all( std::function<bool(typename bitarray3::iterator& it)> func ) { serial_op(func,ALL); }
	/**
	 \~english @brief Loop over cells in serial order.
	 @param[in] func Function that defines how each grid cell is processed. Stop the loop if return true.
	 @param[in] type Type of target cells. ACTIVE or ALL.
	 \~japanese @brief セルをシリアルに処理する。
	 @param[in] func それぞれのセルを処理する関数。\c true を返すと、ループを中断する。
	 @param[in] type ターゲットセルのタイプ. ACTIVE か ALL。
	 */
	inline void interruptible_serial_op( std::function<bool(typename bitarray3::iterator& it)> func, bool type=ALL ) {
		serial_op([func](int dim, int i, int j, int k, typename bitarray3::iterator& it) {
			func(it);
		},type);
	}
	/**
	 \~english @brief Loop over all the active cells in serial order.
	 @param[in] func Function that defines how each grid cell is processed. Stop the loop if return true.
	 \~japanese @brief アクティブなセルをシリアルに処理する。\c true を返すと、ループを中断する。
	 @param[in] func それぞれのセルを処理する関数。
	 */
	inline void interruptible_serial_actives( std::function<bool(int dim, int i, int j, int k, typename bitarray3::iterator& it)> func ) { serial_op(func,ACTIVES); }
	/**
	 \~english @brief Loop over all the cells in serial order.
	 @param[in] func Function that defines how each grid cell is processed. Stop the loop if return true.
	 \~japanese @brief 全てのセルをシリアルに処理する。\c true を返すと、ループを中断する。
	 @param[in] func それぞれのセルを処理する関数。
	 */
	inline void interruptible_serial_all( std::function<bool(int dim, int i, int j, int k, typename bitarray3::iterator& it)> func ) { serial_op(func,ALL); }
	/**
	 \~english @brief Loop over cells in serial order.
	 @param[in] func Function that defines how each grid cell is processed.
	 @param[in] type Type of target cells. ACTIVE or ALL. Stop the loop if return true.
	 \~japanese @brief セルをシリアルに処理する。
	 @param[in] func それぞれのセルを処理する関数。\c true を返すと、ループを中断する。
	 @param[in] type ターゲットセルのタイプ. ACTIVE か ALL。
	 */
	inline void interruptible_serial_op( std::function<bool(int dim, int i, int j, int k, typename bitarray3::iterator& it)> func, bool type=ALL ) {
		for( int dim : DIMS3 ) {
			(*this)[dim].serial_op([&](int i, int j, int k, typename bitarray3::iterator& it) {
				func(dim,i,j,k,it);
			},type);
		}
	}
	/**
	 \~english @brief Loop over all the cells in serial order by read-only fashion.
	 @param[in] func Function that defines how each grid cell is processed. Stop the loop if return true.
	 \~japanese @brief 全てのセルを読み込み可能に限定してシリアルに処理する。
	 @param[in] func それぞれのセルを処理する関数。\c true を返すと、ループを中断する。
	 */
	inline void interruptible_const_serial_all( std::function<bool(const typename bitarray3::const_iterator& it)> func ) const { const_serial_op(func,ALL); }
	/**
	 \~english @brief Loop over cells in serial order by read-only fashion.
	 @param[in] func Function that defines how each grid cell is processed. Stop the loop if return true.
	 @param[in] type Type of target cells. ACTIVE or ALL. Stop the loop if return true.
	 \~japanese @brief セルを読み込み可能に限定してシリアルに処理する。
	 @param[in] func それぞれのセルを処理する関数。\c true を返すと、ループを中断する。
	 @param[in] type ターゲットセルのタイプ. ACTIVE か ALL。
	 */
	inline void interruptible_const_serial_op( std::function<bool(const typename bitarray3::const_iterator& it)> func, bool type=ALL ) const {
		const_serial_op([func](int dim, int i, int j, int k, const typename bitarray3::const_iterator& it) {
			func(it);
		},type);
	}
	/**
	 \~english @brief Loop over all the active cells in serial order by read-only fashion.
	 @param[in] func Function that defines how each grid cell is processed. Stop the loop if return true.
	 \~japanese @brief アクティブなセルを読み込み可能に限定してシリアルに処理する。
	 @param[in] func それぞれのセルを処理する関数。\c true を返すと、ループを中断する。
	 */
	inline void interruptible_const_serial_actives( std::function<bool(int dim, int i, int j, int k)> func ) const {
		const_serial_op([&](int dim, int i, int j, int k, const typename bitarray3::const_iterator& it){ return func(dim,i,j,k); }, ACTIVES); }
	/**
	 \~english @brief Loop over all the cells in serial order by read-only fashion.
	 @param[in] func Function that defines how each grid cell is processed. Stop the loop if return true.
	 \~japanese @brief 全てのセルを読み込み可能に限定してシリアルに処理する。
	 @param[in] func それぞれのセルを処理する関数。\c true を返すと、ループを中断する。
	 */
	inline void interruptible_const_serial_all( std::function<bool(int dim, int i, int j, int k, const typename bitarray3::const_iterator& it)> func ) const { const_serial_op(func,ALL); }
	/**
	 \~english @brief Loop over cells in serial order by read-only fashion.
	 @param[in] func Function that defines how each grid cell is processed. Stop the loop if return true.
	 @param[in] type Type of target cells. ACTIVE or ALL. Stop the loop if return true.
	 \~japanese @brief セルを読み込み可能に限定してシリアルに処理する。
	 @param[in] func それぞれのセルを処理する関数。\c true を返すと、ループを中断する。
	 @param[in] type ターゲットセルのタイプ. ACTIVE か ALL。
	 */
	inline void interruptible_const_serial_op( std::function<bool(int dim, int i, int j, int k, const typename bitarray3::const_iterator& it)> func, bool type=ALL ) const {
		for( int dim : DIMS3 ) {
			(*this)[dim].const_serial_op([&](int i, int j, int k, const typename bitarray3::const_iterator& it) {
				func(dim,i,j,k,it);
			},type);
		}
	}
	/**
	 \~english @brief Dilate cells.
	 @param[in] func Function that specifies what value to assign on dilated cells.
	 @param[in] count Number of dilation count.
	 \~japanese @brief 拡張する。
	 @param[in] func 拡張されたセルにどのような値を与えるか指定する関数。
	 @param[in] count 拡張の回数。
	 */
	void dilate( std::function<void(int dim, int i, int j, int k, typename bitarray3::iterator& it, int thread_index)> func, int count=1 ) {
		while( count -- ) {
			m_parallel.for_each(DIM3,[&]( size_t dim ) {
				operator[](dim).dilate([&](int i, int j, int k, typename bitarray3::iterator& it, int thread_index) {
					func(dim,i,j,k,it,thread_index);
				});
			});
		}
	}
	/**
	 \~english @brief Dilate cells.
	 @param[in] func Function that specifies what value to assign on dilated cells.
	 @param[in] count Number of dilation count.
	 \~japanese @brief 拡張する。
	 @param[in] func 拡張されたセルにどのような値を与えるか指定する関数。
	 @param[in] count 拡張の回数。
	 */
	void dilate( std::function<void(int dim, int i, int j, int k, typename bitarray3::iterator& it)> func, int count=1 ) {
		dilate([&](int dim, int i, int j, int k, typename bitarray3::iterator& it, int thread_index) {
			func(dim,i,j,k,it);
		},count);
	}
	/**
	 \~english @brief Dilate cells.
	 @param[in] count Number of dilation count.
	 \~japanese @brief 拡張する。
	 @param[in] count 拡張の回数。
	 */
	void dilate( int count=1 ) {
		dilate([&](int dim, int i, int j, int k, typename bitarray3::iterator& it){ it.set(); },count);
	}
	/**
	 \~english @brief Erode cells.
	 @param[in] func Function that specifies whether to inactivate the cell.
	 @param[in] count Number of erode count.
	 \~japanese @brief 縮小する。
	 @param[in] func 拡張されたセルを非アクティブにするか指定する関数。
	 @param[in] count 縮小の回数。
	 */
	void erode( std::function<bool(int dim, int i, int j, int k, int thread_index)> func, int count=1 ) {
		while( count -- ) {
			m_parallel.for_each(DIM3,[&]( size_t dim ) {
				operator[](dim).erode([&](int i, int j, int k, int thread_index) {
					return func(dim,i,j,k,thread_index);
				});
			});
		}
	}
	/**
	 \~english @brief Erode cells.
	 @param[in] func Function that specifies whether to inactivate the cell.
	 @param[in] count Number of erode count.
	 \~japanese @brief 縮小する。
	 @param[in] func 拡張されたセルを非アクティブにするか指定する関数。
	 @param[in] count 縮小の回数。
	 */
	void erode( std::function<bool(int dim, int i, int j, int k)> func, int count=1 ) {
		erode([&](int dim, int i, int j, int k, int thread_index) {
			return func(dim,i,j,k);
		},count);
	}
	/**
	 \~english @brief Erode cells.
	 @param[in] count Number of erode count.
	 \~japanese @brief 縮小する。
	 @param[in] count 縮小の回数。
	 */
	void erode( int count=1 ) {
		return erode([&](int dim, int i, int j, int k){ return true; },count);
	}
	/**
	 \~english @brief Set the core name of module of this grid.
	 @param[in] Name of the core name.
	 \~japanese @brief グリッドのモジュールのコアネームを取得する。
	 @param[in] コアネームの名前。
	 */
	void set_core_name( std::string core_name ) {
		m_array_0.set_core_name(core_name);
		m_array_1.set_core_name(core_name);
		m_array_2.set_core_name(core_name);
	}
	/**
	 \~english @brief Get the core name of module of this grid.
	 @return Name of the core name.
	 \~japanese @brief グリッドのモジュールのコアネームを取得する。
	 @return コアネームの名前。
	 */
	std::string get_core_name() const {
		return m_array_0.get_core_name();
	}
	/// \~english @brief Collection of properties of this grid.
	/// \~japanese @brief このグリッドのプロパティー集。
	struct type3 {
		/**
		 \~english @brief Core name of the module.
		 \~japanese @brief モジュールのコアネーム。
		 */
		std::string core_name;
		/**
		 \~english @brief Shape of the grid.
		 \~japanese @brief 格子の形。
		 */
		shape3 shape;
		/**
		 \~english @brief Type for x dimensional grid.
		 \~japanese @brief x 次元のグリッドのタイプ。
		 */
		typename bitarray3::type3 type0;
		/**
		 \~english @brief Type for y dimensional grid.
		 \~japanese @brief y 次元のグリッドのタイプ。
		 */
		typename bitarray3::type3 type1;
		/**
		 \~english @brief Type for z dimensional grid.
		 \~japanese @brief z 次元のグリッドのタイプ。
		 */
		typename bitarray3::type3 type2;
		/**
		 \~english @brief Check equality.
		 @return \c true if equal \c false otherwise.
		 \~japanese @brief 同値の確認。
		 @return 同じなら \c true を、そうでなければ \c false を返す。
		 */
		bool operator==( const type3 &type ) const {
			return
				core_name == type.core_name && shape == type.shape &&
				type0 == type.type0 && type1 == type.type1 && type2 == type.type2;
		}
	};
	/**
	 \~english @brief Get the type of this grid.
	 @return Type of this grid.
	 \~japanese @brief このグリッドの type を取得する。
	 @return このグリッドの type。
	 */
	type3 type() const {
		return { get_core_name(), m_shape, m_array_0.type(), m_array_1.type(), m_array_2.type() };
	}
	/**
	 \~english @brief Set the type of this grid.
	 @param[in] type An instance of type to set.
	 \~japanese @brief グリッドの type を設定する。
	 @param[in] type セットする type のインスタンス。
	 */
	void set_type( const type3 &type ) {
		m_shape = type.shape;
		m_array_0.set_type(type.type0);
		m_array_1.set_type(type.type1);
		m_array_2.set_type(type.type2);
	}
private:
	parallel_driver m_parallel{this};
	bitarray3 m_array_0{this};
	bitarray3 m_array_1{this};
	bitarray3 m_array_2{this};
	shape3 m_shape;
	//
	virtual void initialize( const filestream &file ) override {
		file.r(m_shape);
	}
	virtual void serialize( const filestream &file ) const override {
		file.w(m_shape);
	}
};
//
SHKZ_END_NAMESPACE
//
#endif