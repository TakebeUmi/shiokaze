/*
**	tiledarray3.cpp
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on Feb 7, 2018.
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
#ifndef SHKZ_TILEDARRAY3_H
#define SHKZ_TILEDARRAY3_H
//
#include <vector>
#include <cmath>
#include <cassert>
#include <limits>
#include <shiokaze/array/array_core3.h>
#include "bitcount/bitcount.h"
#include "dilate3.h"
//
SHKZ_BEGIN_NAMESPACE
//
class tiledarray3 : public array_core3 {
public:
	//
	LONG_NAME("Tiled Array 3D")
	ARGUMENT_NAME("TiledArray")
	//
	tiledarray3 () = default;
	virtual ~tiledarray3() {
		dealloc();
	}
	//
protected:
	//
	virtual void configure( configuration &config ) override {
		config.get_unsigned("TileSize",m_Z,"Tile size per dimension");
		assert( m_Z*m_Z*m_Z <= std::numeric_limits<unsigned short>::max());
		assert( m_Z <= std::numeric_limits<unsigned char>::max());
	}
	//
	virtual void initialize( unsigned nx, unsigned ny, unsigned nz, unsigned element_bytes ) override {
		//
		assert( element_bytes <= std::numeric_limits<unsigned char>::max() );
		dealloc();
		//
		m_nx = nx;
		m_ny = ny;
		m_nz = nz;
		//
		m_bx = std::ceil(m_nx/(double)m_Z);
		m_by = std::ceil(m_ny/(double)m_Z);
		m_bz = std::ceil(m_nz/(double)m_Z);
		m_element_bytes = element_bytes;
		m_plane = m_bx*m_by;
		//
		m_tiles.resize(m_bx*m_by*m_bz);
	}
	//
	virtual void get( unsigned &nx, unsigned &ny, unsigned &nz, unsigned &element_bytes ) const override {
		nx = m_nx;
		ny = m_ny;
		nz = m_nz;
		element_bytes = m_element_bytes;
	}
	//
	virtual size_t count( const parallel_driver &parallel ) const override {
		std::vector<size_t> total_slots(parallel.get_thread_num());
		parallel.for_each(m_bx*m_by*m_bz,[&]( size_t n, int thread_index ) {
			if( m_tiles[n] ) {
				total_slots[thread_index] += m_tiles[n]->count();
			}
		});
		size_t total (0);
		for( const auto &e : total_slots ) total += e;
		return total;
	}
	//
	virtual void copy( const array_core3 &array, std::function<void(void *target, const void *src)> copy_func, const parallel_driver &parallel ) override {
		//
		dealloc();
		unsigned nx, ny, nz, element_bytes;
		array.get(nx,ny,nz,element_bytes);
		//
		auto mate_array = dynamic_cast<const tiledarray3 *>(&array);
		if( mate_array ) {
			m_Z = mate_array->m_Z;
			initialize(nx,ny,nz,element_bytes);
			m_fill_mask = mate_array->m_fill_mask;
			for( size_t n=0; n<m_bx*m_by*m_bz; ++n ) {
				if( mate_array->m_tiles[n] ) {
					m_tiles[n] = new chunk3(*mate_array->m_tiles[n],copy_func);
					if( block_filled(n)) m_tiles[n]->fill_all();
				}
			}
		} else {
			initialize(nx,ny,nz,element_bytes);
			array.const_serial_actives([&](int i, int j, int k, const void *src_ptr, const bool& filled) {
				set(i,j,k,[&](void *dst_ptr, bool &active) {
					copy_func(dst_ptr,src_ptr);
					active = true;
				});
				return false;
			});
			if( element_bytes ) {
				array.const_serial_inside([&](int i, int j, int k, const void *src_ptr, const bool &active) {
					unsigned bi = i / m_Z;
					unsigned bj = j / m_Z;
					unsigned bk = k / m_Z;
					size_t n = encode(bi,bj,bk);
					if( m_fill_mask.empty()) m_fill_mask.resize(m_bx*m_by*m_bz);
					if( ! m_tiles[n] ) m_fill_mask[n] = true;
					else m_tiles[n]->set_filled(i-bi*m_Z,j-bj*m_Z,k-bk*m_Z);
					return false;
				});
			}
		}
	}
	//
	void dealloc() {
		for( size_t n=0; n<m_bx*m_by*m_bz; ++n ) {
			if( m_tiles[n] ) {
				delete m_tiles[n];
				m_tiles[n] = nullptr;
			}
		}
		m_fill_mask.resize(0);
	}
	bool block_filled( size_t n ) const {
		return m_fill_mask.empty() ? false : m_fill_mask[n];
	}
	//
	bool check_bound( int i, int j, int k ) const {
		if( i >= 0 && j >= 0 && k >= 0 && i < m_nx && j < m_ny && k < m_nz ) {
			return true;
		} else {
			printf( "Out of bounds (i=%d,j=%d,k=%d), (w=%d,h=%d,d=%d)\n", i, j, k, m_nx, m_ny, m_nz );
			return false;
		}
	}
	//
	virtual void set( int i, int j, int k, std::function<void(void *value_ptr, bool &active)> func ) override {
		//
#if SHKZ_DEBUG
		assert(check_bound(i,j,k));
#endif
		unsigned bi = i / m_Z;
		unsigned bj = j / m_Z;
		unsigned bk = k / m_Z;
		int oi = bi*m_Z;
		int oj = bj*m_Z;
		int ok = bk*m_Z;
		size_t n = encode(bi,bj,bk);
		//
		if( ! m_tiles[n] ) {
			bool active (false);
			unsigned char buffer[m_element_bytes ? m_element_bytes : 1];
			func(m_element_bytes ? buffer : nullptr,active);
			if( active ) {
				unsigned Zx = std::min(m_nx-oi,m_Z);
				unsigned Zy = std::min(m_ny-oj,m_Z);
				unsigned Zz = std::min(m_nz-ok,m_Z);
				m_tiles[n] = new chunk3(oi,oj,ok,Zx,Zy,Zz,m_element_bytes);
				if( block_filled(n)) m_tiles[n]->fill_all();
				m_tiles[n]->set(i-oi,j-oj,k-ok,buffer);
			}
		} else {
			m_tiles[n]->set(i-oi,j-oj,k-ok,func);
			if( m_tiles[n]->deletable()) {
				delete m_tiles[n];
				m_tiles[n] = nullptr;
			}
		}
	}
	//
	virtual const void * operator()( int i, int j, int k, bool &filled ) const override {
		//
#if SHKZ_DEBUG
		assert(check_bound(i,j,k));
#endif
		unsigned bi = i / m_Z;
		unsigned bj = j / m_Z;
		unsigned bk = k / m_Z;
		size_t n = encode(bi,bj,bk);
		//
		filled = false;
		if( ! m_tiles[n] ) {
			filled = block_filled(n);
			return nullptr;
		} else {
			return m_tiles[n]->get(i-bi*m_Z,j-bj*m_Z,k-bk*m_Z,&filled);
		}
	}
	//
	virtual void dilate( std::function<void(int i, int j, int k, void *value_ptr, bool &active, const bool &filled, int thread_index)> func, const parallel_driver &parallel ) override {
		dilate3::dilate<size_t>(this,func,parallel);
	}
	//
	virtual void flood_fill( std::function<bool(void *value_ptr)> inside_func, const parallel_driver &parallel ) override {
		//
		if( ! m_element_bytes ) return;
		//
		parallel.for_each(m_bx*m_by*m_bz,[&]( size_t n ) {
			if( m_tiles[n] ) {
				m_tiles[n]->flood_fill(inside_func);
			}
		});
		//
		m_fill_mask.clear();
		m_fill_mask.resize(m_bx*m_by*m_bz,false);
		std::stack<size_t> start_queue;
		//
		for( size_t n=0; n<m_bx*m_by*m_bz; ++n ) {
			if( m_tiles[n] ) {
				int bi, bj, bk;
				decode(n,bi,bj,bk);
				for( int dim : DIMS3 ) for( int dir=-1; dir<=1; dir+=2 ) {
					int ni(bi+dir*(dim==0)), nj(bj+dir*(dim==1)), nk(bk+dir*(dim==2));
					if( ! shape3(m_bx,m_by,m_bz).out_of_bounds(ni,nj,nk)) {
						size_t m = encode(ni,nj,nk);
						if( ! m_tiles[m] ) {
							if( m_tiles[n]->filled(
								(m_Z-1)*(dir==1)*(dim==0),
								(m_Z-1)*(dir==1)*(dim==1),
								(m_Z-1)*(dir==1)*(dim==2))
							) {
								if( ! m_fill_mask[m] ) {
									start_queue.push(m);
									m_fill_mask[m] = true;
								}
							}
						}
					}
				}
			}
		}
		//
		std::stack<vec3i> queue;
		auto markable = [&]( const vec3i &ni ) {
			if( ! shape3(m_bx,m_by,m_bz).out_of_bounds(ni)) {
				size_t n = encode(ni[0],ni[1],ni[2]);
				return m_fill_mask[n] == false && ( ! m_tiles[n] || m_tiles[n]->count_filled() == 0 );
			} else {
				return false;
			}
		};
		//
		while( ! start_queue.empty()) {
			size_t n = start_queue.top();
			start_queue.pop();
			int i, j, k; decode(n,i,j,k);
			vec3i pi(i,j,k);
			queue.push(pi);
			while( ! queue.empty()) {
				vec3i qi = queue.top();
				size_t m = encode(qi[0],qi[1],qi[2]);
				m_fill_mask[m] = true;
				if( m_tiles[m] && ! m_tiles[m]->count_filled()) {
					m_tiles[m]->fill_all();
				}
				queue.pop();
				for( int dim : DIMS3 ) for( int dir=-1; dir<=1; dir+=2 ) {
					vec3i ni = qi+dir*vec3i(dim==0,dim==1,dim==2);
					if( markable(ni)) queue.push(ni);
				}
			}
		}
		//
#ifdef SHKZ_DEBUG
		parallel.for_each(m_bx*m_by,[&]( size_t n ) {
			if( m_tiles[n] ) assert(m_tiles[n]->debug_verify_active_count());
		});
#endif
	}
	//
	virtual void const_parallel_inside ( std::function<void(int i, int j, int k, const void *value_ptr, const bool &active, int thread_index )> func, const parallel_driver &parallel ) const override {
		//
		parallel.for_each(m_bx*m_by*m_bz,[&]( size_t n, int thread_index ) {
			if( m_tiles[n] ) {
				m_tiles[n]->const_loop_inside([&]( int i, int j, int k, const void *value_ptr, const bool &active ) {
					func(i,j,k,value_ptr,active,thread_index);
					return false;
				});
			} else if( block_filled(n) ) {
				int bi, bj, bk;
				decode(n,bi,bj,bk);
				int oi = m_Z*bi;
				int oj = m_Z*bj;
				int ok = m_Z*bk;
				unsigned Zx = std::min(m_nx-oi,m_Z);
				unsigned Zy = std::min(m_ny-oj,m_Z);
				unsigned Zz = std::min(m_nz-ok,m_Z);
				for( int kk=0; kk<Zz; ++kk ) for( int jj=0; jj<Zy; ++jj ) for( int ii=0; ii<Zx; ++ii ) {
					func(oi+ii,oj+jj,ok+kk,nullptr,false,thread_index);
				}
			}
		});
	}
	virtual void const_serial_inside ( std::function<bool(int i, int j, int k, const void *value_ptr, const bool &active )> func ) const override {
		//
		for( size_t n=0; n<m_bx*m_by*m_bz; ++n ) {
			if( m_tiles[n] ) {
				if(m_tiles[n]->const_loop_inside(func)) break;
			} else if( block_filled(n) ) {
				int bi, bj, bk;
				decode(n,bi,bj,bk);
				int oi = m_Z*bi;
				int oj = m_Z*bj;
				int ok = m_Z*bk;
				unsigned Zx = std::min(m_nx-oi,m_Z);
				unsigned Zy = std::min(m_ny-oj,m_Z);
				unsigned Zz = std::min(m_nz-ok,m_Z);
				//
				bool do_break (false);
				for( int kk=0; kk<Zz; ++kk ) for( int jj=0; jj<Zy; ++jj ) for( int ii=0; ii<Zx; ++ii ) {
					if( func(oi+ii,oj+jj,ok+kk,nullptr,false)) {
						do_break = true;
						goto loop_escape;
					}
				}
loop_escape:
				if( do_break ) break;
			}
		}
	}
	//
	void parallel_loop_actives_body ( int bi, int bj, int bk, std::function<void(int i, int j, int k, void *value_ptr, bool &active, const bool &filled )> func ) {
		loop_actives_body(bi,bj,bk,[&](int i, int j, int k, void *value_ptr, bool &active, const bool &filled ) {
			func(i,j,k,value_ptr,active,filled);
			return false;
		});
	}
	bool loop_actives_body ( int bi, int bj, int bk, std::function<bool(int i, int j, int k, void *value_ptr, bool &active, const bool &filled )> func ) {
		size_t n = encode(bi,bj,bk);
		if( m_tiles[n] ) {
			bool result = m_tiles[n]->loop_actives(func);
			if( m_tiles[n]->deletable()) {
				delete m_tiles[n];
				m_tiles[n] = nullptr;
			}
			if( result ) return true;
		}
		return false;
	}
	//
	void parallel_const_loop_actives_body ( int bi, int bj, int bk, std::function<void(int i, int j, int k, const void *value_ptr, const bool &filled )> func ) const {
		const_loop_actives_body(bi,bj,bk,[&](int i, int j, int k, const void *value_ptr, const bool &filled) {
			func(i,j,k,value_ptr,filled);
			return false;
		});
	}
	bool const_loop_actives_body ( int bi, int bj, int bk, std::function<bool(int i, int j, int k, const void *value_ptr, const bool &filled )> func ) const {
		size_t n = encode(bi,bj,bk);
		if( m_tiles[n] ) {
			if(m_tiles[n]->const_loop_actives(func)) return true;
		}
		return false;
	}
	//
	void parallel_loop_all_body ( int bi, int bj, int bk, std::function<void(int i, int j, int k, void *value_ptr, bool &active, const bool &filled )> func ) {
		loop_all_body(bi,bj,bk,[&](int i, int j, int k, void *value_ptr, bool &active, const bool &filled ) {
			func(i,j,k,value_ptr,active,filled);
			return false;
		});
	}
	bool loop_all_body ( int bi, int bj, int bk, std::function<bool(int i, int j, int k, void *value_ptr, bool &active, const bool &filled )> func ) {
		size_t n = encode(bi,bj,bk);
		if( m_tiles[n] ) {
			bool result = m_tiles[n]->loop_all(func);
			if( m_tiles[n]->deletable()) {
				delete m_tiles[n];
				m_tiles[n] = nullptr;
			}
			if( result ) return true;
		} else {
			unsigned char buffer[m_element_bytes ? m_element_bytes : 1];
			int oi = bi*m_Z;
			int oj = bj*m_Z;
			int ok = bk*m_Z;
			unsigned Zx = std::min(m_Z,m_nx-oi);
			unsigned Zy = std::min(m_Z,m_ny-oj);
			unsigned Zz = std::min(m_Z,m_nz-ok);
			for( int kk=0; kk<Zz; ++kk ) for( int jj=0; jj<Zy; ++jj ) for( int ii=0; ii<Zx; ++ii ) {
				bool active (false);
				int i = oi+ii;
				int j = oj+jj;
				int k = ok+kk;
				func(i,j,k,m_element_bytes ? buffer : nullptr,active,block_filled(n));
				if( active ) {
					if( ! m_tiles[n] ) {
						unsigned Zx = std::min(m_nx-oi,m_Z);
						unsigned Zy = std::min(m_ny-oj,m_Z);
						unsigned Zz = std::min(m_nz-ok,m_Z);
						m_tiles[n] = new chunk3(oi,oj,ok,Zx,Zy,Zz,m_element_bytes);
						if( block_filled(n)) m_tiles[n]->fill_all();
					}
					m_tiles[n]->set(ii,jj,kk,buffer);
				}
			}
		}
#ifdef SHKZ_DEBUG
		if( m_tiles[n] ) assert(m_tiles[n]->debug_verify_active_count());
#endif
		return false;
	}
	//
	void parallel_const_loop_all_body ( int bi, int bj, int bk, std::function<void(int i, int j, int k, const void *value_ptr, const bool &active, const bool &filled )> func ) const {
		const_loop_all_body(bi,bj,bk,[&](int i, int j, int k, const void *value_ptr, const bool &active, const bool &filled) {
			func(i,j,k,value_ptr,active,filled);
			return false;
		});
	}
	bool const_loop_all_body ( int bi, int bj, int bk, std::function<bool(int i, int j, int k, const void *value_ptr, const bool &active, const bool &filled)> func ) const {
		size_t n = encode(bi,bj,bk);
		if( m_tiles[n] ) {
			if(m_tiles[n]->const_loop_all(func)) return true;
		} else {
			int oi = bi*m_Z;
			int oj = bj*m_Z;
			int ok = bk*m_Z;
			unsigned Zx = std::min(m_Z,m_nx-oi);
			unsigned Zy = std::min(m_Z,m_ny-oj);
			unsigned Zz = std::min(m_Z,m_nz-ok);
			for( int kk=0; kk<Zz; ++kk ) for( int jj=0; jj<Zy; ++jj ) for( int ii=0; ii<Zx; ++ii ) {
				bool active (false);
				func(oi+ii,oj+jj,ok+kk,nullptr,active,block_filled(n));
			}
		}
		return false;
	}
	//
	virtual void parallel_actives ( std::function<void(int i, int j, int k, void *value_ptr, bool &active, const bool &filled, int thread_index )> func, const parallel_driver &parallel ) override {
		parallel.for_each(shape3(m_bx,m_by,m_bz),[&](int bi, int bj, int bk, int thread_index) {
			parallel_loop_actives_body(bi,bj,bk,[&](int i, int j, int k, void *value_ptr, bool &active, const bool &filled){
				func(i,j,k,value_ptr,active,filled,thread_index);
			});
		});
	}
	virtual void serial_actives ( std::function<bool(int i, int j, int k, void *value_ptr, bool &active, const bool &filled )> func ) override {
		for( int bk=0; bk<m_bz; ++bk ) for( int bj=0; bj<m_by; ++bj ) for( int bi=0; bi<m_bx; ++bi ) if(loop_actives_body(bi,bj,bk,func)) goto serial_actives_end;
serial_actives_end: ;
	}
	//
	virtual void const_parallel_actives ( std::function<void(int i, int j, int k, const void *value_ptr, const bool &filled, int thread_index )> func, const parallel_driver &parallel ) const override {
		parallel.for_each(shape3(m_bx,m_by,m_bz),[&](int bi, int bj, int bk, int thread_index) {
			parallel_const_loop_actives_body(bi,bj,bk,[&](int i, int j, int k, const void *value_ptr, const bool &filled) {
				func(i,j,k,value_ptr,filled,thread_index);
			});
		});
	}
	virtual void const_serial_actives ( std::function<bool(int i, int j, int k, const void *value_ptr, const bool &filled )> func ) const override {
		for( int bk=0; bk<m_bz; ++bk ) for( int bj=0; bj<m_by; ++bj ) for( int bi=0; bi<m_bx; ++bi ) if(const_loop_actives_body(bi,bj,bk,func)) goto const_serial_actives_end;
const_serial_actives_end: ;
	}
	//
	virtual void parallel_all ( std::function<void(int i, int j, int k, void *value_ptr, bool &active, const bool &filled, int thread_index )> func, const parallel_driver &parallel ) override {
		parallel.for_each(shape3(m_bx,m_by,m_bz),[&](int bi, int bj, int bk, int thread_index) {
			parallel_loop_all_body(bi,bj,bk,[&](int i, int j, int k, void *value_ptr, bool &active, const bool &filled) {
				func(i,j,k,value_ptr,active,filled,thread_index);
			});
		});
	}
	virtual void serial_all ( std::function<bool(int i, int j, int k, void *value_ptr, bool &active, const bool &filled )> func ) override {
		for( int bk=0; bk<m_bz; ++bk ) for( int bj=0; bj<m_by; ++bj ) for( int bi=0; bi<m_bx; ++bi ) if(loop_all_body(bi,bj,bk,func)) goto serial_all_end;
serial_all_end: ;
	}
	//
	virtual void const_parallel_all ( std::function<void(int i, int j, int k, const void *value_ptr, const bool &active, const bool &filled, int thread_index )> func, const parallel_driver &parallel ) const override {
		parallel.for_each(shape3(m_bx,m_by,m_bz),[&](int bi, int bj, int bk, int thread_index) {
			parallel_const_loop_all_body(bi,bj,bk,[&](int i, int j, int k, const void *value_ptr, const bool &active, const bool &filled) {
				func(i,j,k,value_ptr,active,filled,thread_index);
			});
		});
	}
	virtual void const_serial_all ( std::function<bool(int i, int j, int k, const void *value_ptr, const bool &active, const bool &filled )> func ) const override {
		for( int bk=0; bk<m_bz; ++bk ) for( int bj=0; bj<m_by; ++bj ) for( int bi=0; bi<m_bx; ++bi ) if(const_loop_all_body(bi,bj,bk,func)) goto const_serial_all_end;
const_serial_all_end: ;
	}
	//
private:
	//
	struct chunk3 {
		//
		chunk3 ( int oi, int oj, int ok, unsigned Zx, unsigned Zy, unsigned Zz, unsigned element_bytes ) : m_oi(oi), m_oj(oj), m_ok(ok), m_Zx(Zx), m_Zy(Zy), m_Zz(Zz), m_element_bytes(element_bytes) {
			//
			if( m_element_bytes ) {
				m_buffer = new unsigned char [m_Zx*m_Zy*m_Zz*m_element_bytes];
			}
			m_bit_mask_size = std::ceil(m_Zx*m_Zy*m_Zz/8.0);
			m_bit_mask = new unsigned char [m_bit_mask_size];
			std::memset(m_buffer,0,m_element_bytes*m_Zx*m_Zy*m_Zz);
			std::memset(m_bit_mask,0,m_bit_mask_size);
			m_num_active = 0;
		}
		~chunk3 () {
			if( m_buffer ) delete [] m_buffer;
			delete [] m_bit_mask;
			if( m_fill_mask ) delete [] m_fill_mask;
		}
		chunk3( const chunk3 &instance, std::function<void(void *target, const void *src)> copy_func ) {
			//
			m_num_active = instance.m_num_active;
			m_Zx = instance.m_Zx;
			m_Zy = instance.m_Zy;
			m_Zz = instance.m_Zz;
			m_oi = instance.m_oi;
			m_oj = instance.m_oj;
			m_ok = instance.m_ok;
			m_element_bytes = instance.m_element_bytes;
			//
			m_bit_mask_size = instance.m_bit_mask_size;
			m_bit_mask = new unsigned char [m_bit_mask_size];
			std::memcpy(m_bit_mask,instance.m_bit_mask,m_bit_mask_size);
			if( instance.m_fill_mask ) {
				m_fill_mask = new unsigned char [m_bit_mask_size];
				std::memcpy(m_fill_mask,instance.m_fill_mask,m_bit_mask_size);
			}
			//
			if( m_element_bytes ) {
				size_t size = m_Zx*m_Zy*m_Zz*m_element_bytes;
				m_buffer = new unsigned char [size];
			}
			//
			for( int kk=0; kk<m_Zz; ++kk ) for( int jj=0; jj<m_Zy; ++jj ) for( int ii=0; ii<m_Zx; ++ii ) {
				size_t n = encode(ii,jj,kk);
				unsigned char &mask = *(m_bit_mask+(n>>3));
				if( mask ) {
					if((mask >> (n&7)) & 1U) {
						size_t offset = n*m_element_bytes;
						if( m_element_bytes ) copy_func(m_buffer+offset,instance.m_buffer+offset);
					}
				}
			}
			//
		}
		bool debug_verify_active_count() const {
			//
			size_t verify_count (0);
			for( size_t n8=0; n8<m_bit_mask_size; ++n8 ) {
				if( m_bit_mask[n8] == 0xFF ) verify_count += 8;
				else if( m_bit_mask[n8] ) {
					for( size_t n=8*n8; n<8*(n8+1); ++n ) {
						if( (*(m_bit_mask+(n>>3)) >> (n&7)) & 1U ) ++ verify_count;
					}
				}
			}
			if(verify_count != m_num_active) {
				printf( "===== Verification failed! =====\n");
				printf( "verify_count = %d, m_num_active = %d\n", (int)verify_count, (int)m_num_active );
				printf( "================================\n");
				return false;
			}
			return true;
		}
		void alloc_fill( unsigned char with_value ) {
			m_fill_mask = new unsigned char [m_bit_mask_size];
			std::memset(m_fill_mask,with_value,m_bit_mask_size);
		}
		size_t count() const {
			return bitcount::count(m_bit_mask,m_bit_mask_size);
		}
		size_t count_filled() const {
			return m_fill_mask ? bitcount::count(m_fill_mask,m_bit_mask_size) : 0;
		}
		void fill_all() {
			if( ! m_fill_mask ) alloc_fill(0xFF);
			else std::memset(m_fill_mask,0xFF,m_bit_mask_size);
		}
		void set( int bi, int bj, int bk, const void *value_ptr ) {
			//
			set(bi,bj,bk,[&](void *target_ptr, bool &active) {
				if( value_ptr ) {
					std::memcpy(target_ptr,value_ptr,m_element_bytes);
					active = true;
				} else {
					active = false;
				}
			});
		}
		void set( int bi, int bj, int bk, std::function<void(void *value_ptr, bool &active)> func ) {
			//
			size_t n = encode(bi,bj,bk);
			unsigned char &mask = *(m_bit_mask+(n>>3));
			bool active = (mask >> (n&7)) & 1U;
			unsigned char *ptr = m_buffer ? m_buffer+n*m_element_bytes : nullptr;
			//
			if( active ) {
				func(ptr,active);
				if( ! active ) {
					assert(m_num_active);
					m_num_active --;
					mask &= ~(1UL << (n&7));
				}
			} else {
				func(ptr,active);
				if( active ) {
					m_num_active ++;
					mask |= 1UL << (n&7);
				}
			}
		}
		void set_filled( int bi, int bj, int bk ) {
			//
			if( ! m_fill_mask ) alloc_fill(0);
			size_t n = encode(bi,bj,bk);
			*(m_fill_mask+(n>>3)) |= 1UL << (n&7);
		}
		void flood_fill( std::function<bool(void *value_ptr)> inside_func ) {
			//
			if( ! m_fill_mask ) alloc_fill(0);
			else std::memset(m_fill_mask,0,m_bit_mask_size);
			//
			std::stack<vec3i> queue;
			shape3 local_shape = shape3(m_Zx,m_Zy,m_Zz);
			shape3 global_shape = shape3(m_oi+m_Zx,m_oj+m_Zy,m_ok+m_Zz);
			auto markable = [&]( vec3i pi, bool default_result ) {
				if( ! local_shape.out_of_bounds(pi) && ! global_shape.out_of_bounds(pi+vec3i(m_oi,m_oj,m_ok))) {
					auto pass_fill_mask = [&]( size_t n ) {
						return ! ((*(m_fill_mask+(n>>3)) >> (n&7)) & 1U);
					};
					const size_t n = encode(pi[0],pi[1],pi[2]);
					if( (*(m_bit_mask+(n>>3)) >> (n&7)) & 1U ) {
						return inside_func(m_buffer ? m_buffer+n*m_element_bytes : nullptr) && pass_fill_mask(n);
					} else {
						return default_result && pass_fill_mask(n);
					}
				} else {
					return false;
				}
			};
			auto mark = [&]( size_t n ) {
				*(m_fill_mask+(n>>3)) |= 1UL << (n&7);
			};
			size_t count = local_shape.count();
			for( size_t n8=0; n8<m_bit_mask_size; ++n8 ) {
				if( *(m_bit_mask+n8) ) {
					for( size_t n=8*n8; n < 8*(n8+1); ++n ) if ( n < count ) {
						int bi, bj, bk; decode(n,bi,bj,bk);
						vec3i pi(bi,bj,bk);
						if( markable(pi,false)) {
							if( (*(m_bit_mask+n8) >> (n&7)) & 1U ) {
								queue.push(pi);
								while(! queue.empty()) {
									vec3i qi = queue.top();
									mark(encode(qi[0],qi[1],qi[2]));
									queue.pop();
									for( int dim : DIMS3 ) for( int dir=-1; dir<=1; dir+=2 ) {
										vec3i ni = qi+dir*vec3i(dim==0,dim==1,dim==2);
										if( markable(ni,true)) queue.push(ni);
									}
								}
							}
						}
					}
				}
			}
			//
		}
		bool const_loop_inside( std::function<bool(int i, int j, int k, const void *value_ptr, const bool &active )> func ) const {
			//
			if( m_fill_mask ) {
				for( int kk=0; kk<m_Zz; ++kk ) for( int jj=0; jj<m_Zy; ++jj ) for( int ii=0; ii<m_Zx; ++ii ) {
					size_t n = encode(ii,jj,kk);
					unsigned char &mask = *(m_fill_mask+(n>>3));
					if( (mask >> (n&7)) & 1U ) {
						bool active = (*(m_bit_mask+(n>>3)) >> (n&7)) & 1U;
						int i = m_oi+ii;
						int j = m_oj+jj;
						int k = m_ok+kk;
						if(func(i,j,k,active ? (m_buffer ? m_buffer+n*m_element_bytes : nullptr) : nullptr,active)) return true;
					}
				}
			}
			return false;
		}
		bool loop_actives( std::function<bool(int i, int j, int k, void *value_ptr, bool &active, const bool &filled )> func ) {
			//
			for( int kk=0; kk<m_Zz; ++kk ) for( int jj=0; jj<m_Zy; ++jj ) for( int ii=0; ii<m_Zx; ++ii ) {
				size_t n = encode(ii,jj,kk);
				unsigned char &mask = *(m_bit_mask+(n>>3));
				if( mask ) {
					bool active = (mask >> (n&7)) & 1U;
					if( active ) {
						int i = m_oi+ii;
						int j = m_oj+jj;
						int k = m_ok+kk;
						bool result = func(i,j,k,m_buffer ? m_buffer+n*m_element_bytes : nullptr,active,filled(n));
						if( ! active ) {
							assert(m_num_active);
							m_num_active --;
							mask &= ~(1UL << (n&7));
						}
						if( result ) return true;
					}
				}
			}
#ifdef SHKZ_DEBUG
			assert(debug_verify_active_count());
#endif
			return false;
		}
		bool const_loop_actives( std::function<bool(int i, int j, int k, const void *value_ptr, const bool &filled )> func ) const {
			//
			for( int kk=0; kk<m_Zz; ++kk ) for( int jj=0; jj<m_Zy; ++jj ) for( int ii=0; ii<m_Zx; ++ii ) {
				size_t n = encode(ii,jj,kk);
				unsigned char &mask = *(m_bit_mask+(n>>3));
				if( mask ) {
					if( (mask >> (n&7)) & 1U ) {
						int i = m_oi+ii;
						int j = m_oj+jj;
						int k = m_ok+kk;
						if(func(i,j,k,m_buffer ? m_buffer+n*m_element_bytes : nullptr,filled(n))) return true;
					}
				}
			}
			return false;
		}
		bool loop_all( std::function<bool(int i, int j, int k, void *value_ptr, bool &active, const bool &filled )> func ) {
			//
			for( int kk=0; kk<m_Zz; ++kk ) for( int jj=0; jj<m_Zy; ++jj ) for( int ii=0; ii<m_Zx; ++ii ) {
				size_t n = encode(ii,jj,kk);
				unsigned char &mask = *(m_bit_mask+(n>>3));
				bool active = (mask >> (n&7)) & 1U;
				bool new_active (active);
				int i = m_oi+ii;
				int j = m_oj+jj;
				int k = m_ok+kk;
				bool result = func(i,j,k,m_buffer ? m_buffer+n*m_element_bytes : nullptr,new_active,filled(n));
				if( new_active != active ) {
					if( new_active ) {
						m_num_active ++;
						mask |= 1UL << (n&7);
					} else {
						assert(m_num_active);
						m_num_active --;
						mask &= ~(1UL << (n&7));
					}
				}
				if( result ) return true;
			}
#ifdef SHKZ_DEBUG
			assert(debug_verify_active_count());
#endif
			return false;
		}
		bool const_loop_all( std::function<bool(int i, int j, int k, const void *value_ptr, const bool &active, const bool &filled )> func ) const {
			//
			for( int kk=0; kk<m_Zz; ++kk ) for( int jj=0; jj<m_Zy; ++jj ) for( int ii=0; ii<m_Zx; ++ii ) {
				size_t n = encode(ii,jj,kk);
				const unsigned char &mask = *(m_bit_mask+(n>>3));
				bool active = (mask >> (n&7)) & 1U;
				int i = m_oi+ii;
				int j = m_oj+jj;
				int k = m_ok+kk;
				if(func(i,j,k,m_buffer ? m_buffer+n*m_element_bytes : nullptr,active,filled(n))) return true;
			}
			return false;
		}
		const void* get( int bi, int bj, int bk, bool *_filled=nullptr ) const {
			//
			size_t n = encode(bi,bj,bk);
			unsigned char &mask = *(m_bit_mask+(n>>3));
			if( _filled ) *_filled = filled(n);
			static char tmp;
			if( (mask >> (n&7)) & 1U ) return m_buffer ? m_buffer+n*m_element_bytes : (void *)&tmp;
			else return nullptr;
		}
		bool filled( size_t n ) const {
			if( m_fill_mask ) {
				return (*(m_fill_mask+(n>>3)) >> (n&7)) & 1U;
			} else {
				return false;
			}
		}
		bool filled( int bi, int bj, int bk ) const {
			return filled(encode(bi,bj,bk));
		}
		bool deletable () const {
			return m_num_active == 0;
		}
		//
		unsigned short m_num_active {0};
		unsigned m_oi {0}, m_oj {0}, m_ok {0};
		unsigned char m_Zx {0}, m_Zy {0}, m_Zz {0};
		unsigned char m_element_bytes {0};
		unsigned short m_bit_mask_size {0};
		unsigned char *m_buffer {nullptr};
		unsigned char *m_bit_mask {nullptr};
		unsigned char *m_fill_mask {nullptr};
		//
		size_t encode ( int i, int j, int k ) const {
			return i + j * m_Zx + k * (m_Zx * m_Zy);
		}
		void decode ( size_t n, int &i, int &j, int &k ) const {
			unsigned plane_num = m_Zx * m_Zy;
			i = (n % plane_num) % m_Zx;
			j = (n % plane_num) / m_Zx;
			k = n / plane_num;
		};
	};
	//
	std::vector<chunk3 *> m_tiles;
	std::vector<bool> m_fill_mask;
	unsigned m_nx {0}, m_ny {0}, m_nz {0}, m_bx {0}, m_by {0}, m_bz {0}, m_element_bytes {0};
	unsigned m_Z {16};
	size_t m_plane {0};
	//
	size_t encode ( int bi, int bj, int bk ) const {
		return bi + bj * m_bx + bk * m_plane;
	};
	void decode ( size_t n, int &bi, int &bj, int &bk ) const {
		bi = (n % m_plane) % m_bx;
		bj = (n % m_plane) / m_bx;
		bk = n / m_plane;
	};
};
//
extern "C" module * create_instance() {
	return new tiledarray3();
}
//
extern "C" const char *license() {
	return "BSD-{2,3}-Clause";
}
//
SHKZ_END_NAMESPACE
//
#endif