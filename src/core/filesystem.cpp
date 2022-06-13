/*
**	filesystem.cpp
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
#include <shiokaze/core/filesystem.h>
#include <boost/filesystem.hpp>
#include <iostream>
#include <zlib.h>
//
SHKZ_USING_NAMESPACE
//
bool filesystem::is_exist( std::string path ) {
	return boost::filesystem::exists(path);
}
//
bool filesystem::has_parent( std::string path ) {
	const auto p = boost::filesystem::path(path);
	return p.has_parent_path();
}
//
std::string filesystem::parent_path( std::string path ) {
	const auto p = boost::filesystem::path(path);
	return p.parent_path().string();
}
//
bool filesystem::has_root( std::string path ) {
	const auto p = boost::filesystem::path(path);
	return p.has_root_path();
}
//
bool filesystem::create_directory( std::string path ) {
	return boost::filesystem::create_directory(path);
}
//
bool filesystem::create_directories( std::string path ) {
	return boost::filesystem::create_directories(path);
}
//
void filesystem::rename( std::string old_path, std::string new_path ) {
	boost::filesystem::rename(old_path,new_path);
}
//
void filesystem::remove_file( std::string path ) {
	 boost::filesystem::remove(path);
}
//
void filesystem::remove_dir_contents( std::string path ) {
	 boost::filesystem::remove_all(path);
}
//
std::string filesystem::find_resource_path( std::string directory, std::string name ) {
	//
	std::string path = "./resources/"+ directory + "/" + name;
	if(filesystem::is_exist(path)) return path;
	else {
		std::string private_path = "./private-resources/"+ directory + "/" + name;
		if(filesystem::is_exist(private_path)) return private_path;
		else {
			std::cout << "Could not locate the resource path=" << directory << "/" << name << std::endl;
			exit(0);
		}
	}
	return std::string();
}
//
std::string filesystem::resolve_libname( std::string name ) {
	//
#ifdef __APPLE__
	static const std::string extension ("dylib");
#else
	static const std::string extension ("so");
#endif
	std::string libname = std::string("lib")+name+SHKZ_SUFFIX+"."+extension;
	return libname;
}
//
filestream::filestream( std::string path, OPTION option ) {
	//
	if( option == READ ) {
		m_readable = true;
	} else if( option == WRITE ) {
		m_writable = true;
	}
	if( m_writable ) {
		m_gzfp = (void *)gzopen(path.c_str(),"wb");
	} else if( m_readable ) {
		m_gzfp = (void *)gzopen(path.c_str(),"rb");
	}
	assert(m_gzfp);
}
//
filestream::~filestream() {
	gzclose(reinterpret_cast<gzFile>(m_gzfp));
}
//
void filestream::write( const void *ptr, size_t size ) const {
	assert(m_writable);
	gzwrite(reinterpret_cast<gzFile>(m_gzfp),ptr,size);
}
//
template<> void filestream::write<bool>( const std::vector<bool> &v, std::function<void(const filestream &file, const bool& e)> func ) const {
	size_t size = v.size();
	size_t size8 = (size+7)/8;
	w(size);
	for( size_t n=0; n<size8; ++n ) {
		unsigned char byte (0);
		for( char i=0; i<8; ++i ) {
			size_t idx (8*n+i);
			if( idx < size) {
				if( v[idx] ) byte = byte | (0x01U << i);
			}
		}
		w(byte);
	}
}
//
void filestream::read( void *ptr, size_t size ) const {
	assert(m_readable);
	gzread(reinterpret_cast<gzFile>(m_gzfp),ptr,size);
}
//
template<> void filestream::read<bool>( std::vector<bool> &v, std::function<void(const filestream &file, bool& e)> func ) const {
	size_t size;
	r(size);
	v.resize(size);
	size_t size8 = (size+7)/8;
	for( size_t n=0; n<size8; ++n ) {
		unsigned char byte (0);
		r(byte);
		for( char i=0; i<8; ++i ) {
			size_t idx (8*n+i);
			if( idx < size) {
				v[idx] = byte & (0x01U << i);
			}
		}
	}
}
//
void filestream::write( const std::string &str ) const {
	//
	assert(m_writable);
	size_t size = str.size();
	write(&size,sizeof(size_t));
	for( size_t n=0; n<size; ++n ) write(str.c_str(),size*sizeof(unsigned char));
}
//
void filestream::read( std::string &str ) const {
	//
	assert(m_readable);
	size_t size;
	read(&size,sizeof(size_t));
	std::vector<unsigned char> buffer;
	buffer.resize(size+1);
	buffer[size] = 0;
	char *ptr = reinterpret_cast<char *>(buffer.data());
	read(ptr,size);
	str = std::string((const char *)ptr);
}
//