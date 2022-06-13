/*
**	filesystem.h
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on Feb 6, 2017.
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
#ifndef SHKZ_FILESYSTEM_H
#define SHKZ_FILESYSTEM_H
//
#include <shiokaze/core/common.h>
#include <string>
#include <vector>
#include <functional>
//
SHKZ_BEGIN_NAMESPACE
//
/** @file */
/// \~english @brief Class that perform file system related tasks.
/// \~japanese @brief ファイルシステム関連の処理を行うクラス。
class filesystem {
public:
	/**
	 \~english @brief Check if the file to a path exists
	 @param[in] path Path to check.
	 @return \c true if exists \c false otherwise.
	 \~japanese @brief パスで指定されたファイルが存在するか取得する。
	 @param[in] path 確認するパス。
	 @return 存在すれば \c true 存在していなければ \c false が返る。
	 */
	static bool is_exist( std::string path );
	/**
	 \~english @brief Check if a path contains a parent directory.
	 @param[in] path Path to check.
	 @return \c true if exists \c false otherwise.
	 \~japanese @brief パスに親ディレクトリが存在するか確認する。
	 @param[in] path 確認するパス。
	 @return 存在すれば \c true 存在していなければ \c false が返る。
	 */
	static bool has_parent( std::string path );
	/**
	 \~english @brief Get the path to the parent directory.
	 @param[in] path Full path.
	 @return Path to the parent directory.
	 \~japanese @brief 親ディレクトリへのパスを得る。
	 @param[in] path 完全なパス。
	 @return 親ディレクトリへのパス。
	 */
	static std::string parent_path( std::string path );
	/**
	 \~english @brief Check if a path contains a root path.
	 @param[in] path Path to check.
	 @return \c true if exists \c false otherwise.
	 \~japanese @brief パスにルートパスが存在するか確認する。
	 @param[in] path 確認するパス。
	 @return 存在すれば \c true 存在していなければ \c false が返る。
	 */
	static bool has_root( std::string path );
	/**
	 \~english @brief Create directory to a path.
	 @return \c true if successful \c false otherwise.
	 \~japanese @brief ディレクトリを指定されたパスに作成する。
	 @return 成功すれば \c true 失敗すれば \c false が返る。
	 */
	static bool create_directory( std::string path );
	/**
	 \~english @brief Recursively create directory to a path.
	 @return \c true if successful \c false otherwise.
	 \~japanese @brief ディレクトリを指定されたパスに再帰的に作成する。
	 @return 成功すれば \c true 失敗すれば \c false が返る。
	 */
	static bool create_directories( std::string path );
	/**
	 \~english @brief Rename the path (equivalently moving a file/directory)
	 @param[in] old_path Old path.
	 @param[in] new_path New path.
	 \~japanese @brief パスを変更する (ファイル/ディレクトリの移動)
	 @param[in] old_path 古いパス。
	 @param[in] new_path 新しいパス。
	 */
	static void rename( std::string old_path, std::string new_path );
	/**
	 \~english @brief Delete the file to a path.
	 @param[in] path Path to the file to delete.
	 \~japanese @brief パスで指定されたファイルを削除する。
	 @param[in] path 削除するファイルへのパス。
	 */
	static void remove_file( std::string path );
	/**
	 \~english @brief Delete all the files in a directory.
	 @param[in] path Path to the directory.
	 \~japanese @brief パスで指定されたディレクトリの中身を全て削除する。
	 @param[in] path ディレクトリのパス。
	 */
	static void remove_dir_contents( std::string path );
	/**
	 \~english @brief Get a file path in the resource directory.
	 @param[in] directory Name of the directory
	 @param[in] name Name of the file.
	 \~japanese @brief リソースディレクトリからファイルのパスを取得する。
	 @param[in] directory ディレクトリの名前。
	 @param[in] name ファイルの名前。
	 */
	static std::string find_resource_path( std::string directory, std::string name );
	/**
	 \~english @brief Get a file path in the dynamic library directory.
	 @param[in] name Name of the dynamic library without "lib" prefix and ".so" suffix.
	 \~japanese @brief ライブラリのパスを取得する。
	 @param[in] name "lib" プレフィックスと ".so" サフィックスを除いたライブラリの名前。
	 */
	static std::string resolve_libname( std::string name );
};
/** @file */
/// \~english @brief Class that perform data writing and reading.
/// \~japanese @brief ファイルの書き込みと入力を行うクラス。
class filestream {
public:
	/**
	 \~english @brief Option for reading and writing.
	 \~japanese @brief 読み書きのオプション。
	 */
	enum OPTION { READ, WRITE };
	/**
	 \~english @brief Open a file.
	 @param[in] path Path to a file.
	 @param[in] OPTION Opening option.
	 \~japanese @brief ファイルを開く。
	 @param[in] path ファイルへのパス。
	 @param[in] OPTION ファイルを開く時のオプション。
	 */
	filestream( std::string path, OPTION );
	/**
	 \~english @brief ファイルストリームのデストラクタ。
	 \~japanese @brief Destructor.
	 */
	virtual ~filestream();
	/**
	 \~english @brief Write data.
	 @param[in] ptr Pointer to the data.
	 @param[in] size Byte length.
	 \~japanese @brief データを書き込む。
	 @param[in] ptr データの先頭のポインタ。
	 @param[in] size バイト長。
	 */
	void write( const void *ptr, size_t size ) const;
	/**
	 \~english @brief Read data.
	 @param[out] ptr Pointer to the data.
	 @param[in] size Byte length.
	 \~japanese @brief データを読み込む。
	 @param[out] ptr データの先頭のポインタ。
	 @param[in] size バイト長。
	 */
	void read( void *ptr, size_t size ) const;
	/**
	 \~english @brief Write an element.
	 @param[in] v Reference to an element.
	 \~japanese @brief 要素を書き込む。
	 @param[in] v 書き込む要素への参照。
	 */
	template<class T> void w( const T &v ) const {
		write(&v,sizeof(T));
	}
	/**
	 \~english @brief Read an element.
	 @param[out] v Reference to an element to read.
	 \~japanese @brief 要素を読み込む。
	 @param[out] v 読み込む要素への参照。
	 */
	template<class T> void r( T &v ) const {
		read(&v,sizeof(T));
	}
	/**
	 \~english @brief Write an element vector.
	 @param[in] v Reference to an element vector of vector.
	 @param[in] func Delegate function for writing data.
	 \~japanese @brief 要素を書き込む。
	 @param[in] v 書き込むベクトルのベクトル要素への参照。
	 @param[in] func 要素を書き込むための代理関数。
	 */
	template<class T> void write( const std::vector<std::vector<T> > &v, std::function<void(const filestream &file, const T& e)> func=nullptr ) const {
		const size_t size = v.size();
		w(size);
		for( size_t n=0; n<size; ++n ) {
			write(v[n],func);
		}
	}
	/**
	 \~english @brief Read an element vector.
	 @param[out] v Reference to an element vector of vector.
	 @param[in] func Delegate function for reading data.
	 \~japanese @brief 要素を書き込む。
	 @param[out] v 読み込むベクトルのベクトル要素への参照。
	 @param[in] func 要素を読み込むための代理関数。
	 */
	template<class T> void read( std::vector<std::vector<T> > &v, std::function<void(const filestream &file, T& e)> func=nullptr ) const {
		size_t size;
		r(size);
		v.resize(size);
		for( size_t n=0; n<size; ++n ) {
			read(v[n],func);
		}
	}
	/**
	 \~english @brief Write an element vector.
	 @param[in] v Reference to an element vector.
	 @param[in] func Delegate function for writing data.
	 \~japanese @brief 要素を書き込む。
	 @param[in] v 書き込むベクトル要素への参照。
	 @param[in] func 要素を書き込むための代理関数。
	 */
	template<class T> void write( const std::vector<T> &v, std::function<void(const filestream &file, const T& e)> func=nullptr ) const {
		const size_t size = v.size();
		w(size);
		if( func ) {
			for( size_t n=0; n<size; ++n ) func(*this,v[n]);
		} else {
			write(v.data(),size*sizeof(T));
		}
	}
	/**
	 \~english @brief Read an element vector.
	 @param[out] v Reference to an element vector.
	 @param[in] func Delegate function for reading data.
	 \~japanese @brief 要素を書き込む。
	 @param[out] v 読み込むベクトル要素への参照。
	 @param[in] func 要素を読み込むための代理関数。
	 */
	template<class T> void read( std::vector<T> &v, std::function<void(const filestream &file, T& e)> func=nullptr ) const {
		size_t size;
		r(size);
		v.resize(size);
		if( func ) {
			for( size_t n=0; n<size; ++n ) func(*this,v[n]);
		} else {
			read(v.data(),size*sizeof(T));
		}
	}
	/**
	 \~english @brief Write a string.
	 @param[in] str Reference to an element.
	 \~japanese @brief 文字列を書き込む。
	 @param[in] str 書き込む文字列への参照。
	 */
	void write( const std::string &str ) const;
	/**
	 \~english @brief Read a string.
	 @param[out] str Reference to an output string.
	 \~japanese @brief 文字列を読み込む。
	 @param[out] str 読み込まれた文字列。
	 */
	void read( std::string &str ) const;
	/**
	 \~english @brief Write boundary for debugging purpose.
	 \~japanese @brief デバッグ目的のための境界情報を書き込む。
	 */
	void write_boundary() const {
		const unsigned char boundary ('B');
		w(boundary);
	}
	/**
	 \~english @brief Read boundary for debugging purpose and return if it was a boundary.
	 @return \true if boundary is read \false otherwise.
	 \~japanese @brief デバッグ目的のための境界情報を読み込んで境界が読み込まれたどうが確認する。
	 @return もし境界が読み込まれたら \true を、そうでなければ \false を返す。
	 */
	bool check_boundary() const {
		unsigned char boundary;
		r(boundary);
		return boundary == 'B';
	}
	//
private:
	void *m_gzfp {nullptr};
	bool m_writable {false};
	bool m_readable {false};
};
//
template<> void filestream::write<bool>( const std::vector<bool> &v, std::function<void(const filestream &file, const bool& e)> func ) const;
template<> void filestream::read<bool>( std::vector<bool> &v, std::function<void(const filestream &file, bool& e)> func ) const;
//
SHKZ_END_NAMESPACE
//
#endif