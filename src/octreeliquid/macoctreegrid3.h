/*
**	macoctreegrid3.h
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on November 11, 2019. All rights reserved.
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
#ifndef SHKZ_OCTREEGRID3_H
#define SHKZ_OCTREEGRID3_H
//
#include <shiokaze/math/vec.h>
#include <shiokaze/array/bitarray3.h>
#include <shiokaze/array/array3.h>
#include <shiokaze/array/macarray3.h>
#include <shiokaze/visualizer/gridvisualizer3_interface.h>
#include <shiokaze/utility/macutility3_interface.h>
#include <shiokaze/utility/meshutility3_interface.h>
#include <shiokaze/graphics/graphics_engine.h>
#include <vector>
#include <memory>
//
SHKZ_BEGIN_NAMESPACE
//
namespace macotreeliquid3_namespace {
	//
	using uint_type = uint32_t;
	//
	struct cell_id3 {
		//
		char depth;
		vec3i pi;
		uint_type index;
		bool operator<(const cell_id3& rhs) const { return index < rhs.index; }
		bool operator==(const cell_id3& rhs) const { return index == rhs.index; }
	};
	//
	struct face_id3 {
		//
		char depth;
		char dim;
		vec3i pi;
		uint_type index;
		bool operator<(const face_id3& rhs) const { return index < rhs.index; }
		bool operator==(const face_id3& rhs) const { return index == rhs.index; }
	};
	//
	struct interpolation_data3 {
		vec3d p;
		double coef[4];
		uint_type indices[4];
		unsigned char count {0};
	};
	//
	const static unsigned char shrink_left = 0x01 << 0;
	const static unsigned char shrink_right = 0x01 << 1;
	const static unsigned char shrink_top = 0x01 << 2;
	const static unsigned char shrink_bottom = 0x01 << 3;
	const static unsigned char shrink_front = 0x01 << 4;
	const static unsigned char shrink_back = 0x01 << 5;
	//
	struct ghost_cell3 {
		//
		std::vector<interpolation_data3> data;
		vec3d p;
		std::vector<std::pair<uint_type,double> > combined;
		unsigned char shrink_info {0};
		//
		void initialize( const filestream &file ) {
			file.r(data);
			file.r(p);
			size_t size;
			file.r(size);
			combined.resize(size);
			for( unsigned n=0; n<size; ++n ) {
				file.r(combined[n].first);
				file.r(combined[n].second);
			}
			file.r(shrink_info);
		}
		void serialize( const filestream &file ) const {
			file.w(data);
			file.w(p);
			size_t size = combined.size();
			file.w(size);
			for( unsigned n=0; n<size; ++n ) {
				file.w(combined[n].first);
				file.w(combined[n].second);
			}
			file.w(shrink_info);
		}
	};
	struct ghost_face3 : public ghost_cell3 {};
	//
	struct layer3 : public recursive_configurable {
		//
		shape3 shape;
		double dx;
		//
		bitarray3 fill_flags{this};
		array3<uint_type> active_cells{this};
		macarray3<uint_type> active_faces{this};
		array3<uint_type> ghost_cell_indices{this};
		macarray3<uint_type> ghost_face_indices{this};
		//
		macutility3_driver macutility{this,"macutility3"};
		//
		virtual void post_initialize( bool initialized_from_file ) override {
			//
			if( ! initialized_from_file ) {
				fill_flags.initialize(shape);
				active_cells.initialize(shape);
				active_faces.initialize(shape);
				ghost_cell_indices.initialize(shape);
				ghost_face_indices.initialize(shape);
			}
		}
		//
	private:
		//
		virtual void initialize( const filestream &file ) override {
			file.r(shape);
			file.r(dx);
		}
		virtual void serialize( const filestream &file ) const override {
			file.w(shape);
			file.w(dx);
		}
	};
	//
	struct flux_boundary_condition3 {
		//
		Real velocity[DIM3][2];
		//
		flux_boundary_condition3 () {
			for( int dim : DIMS3 ) {
				velocity[dim][0] = velocity[dim][1] = 0.0;
			}
		}
		//
		bool has_flux () const {
			bool result (false);
			for( int dim : DIMS3 ) if( velocity[dim][0] || velocity[dim][1] ) {
				result = true;
				break;
			}
			return result;
		}
	};
	//
	struct grid3 : public recursive_configurable, public credit {
		//
		LONG_NAME("MAC Octee Grid 3D")
		ARGUMENT_NAME("macoctreegrid")
		//
		grid3 ( recursive_configurable *parent ) {
			if( parent ) parent->add_child(this);
			else setup_now();
		}
		virtual void configure( configuration &config ) override;
		void copy( const grid3 &grid );
		//
		std::vector<std::shared_ptr<layer3> > layers;
		std::vector<Real> levelset;
		std::vector<Real> velocity;
		std::vector<Real> area;
		std::vector<bool> solid_cell;
		std::vector<std::vector<ghost_cell3> > ghost_cells;
		std::vector<std::vector<ghost_face3> > ghost_faces;
		flux_boundary_condition3 flux_boundary_condition;
		//
		uint_type cell_count;
		uint_type face_count;
		uint_type ghost_cell_reference_count;
		uint_type ghost_face_reference_count;
		//
		parallel_driver parallel{this};
		meshutility3_driver meshutility{this,"meshutility3"};
		//
		std::vector<uint_type> cell_map, face_map;
		uint_type valid_cell_count, valid_face_count;
		//
		virtual void initialize( const filestream &file ) override {
			clear();
			size_t size;
			file.r(size);
			for( unsigned n=0; n<size; ++n ) {
				auto layer = std::shared_ptr<layer3>(new layer3());
				layer->setup_now(file);
				layers.push_back(layer);
			}
			//
			std::function<void(const filestream &file, ghost_cell3& e)> func0 = []( const filestream &file, auto &e ) { e.initialize(file); };
			std::function<void(const filestream &file, ghost_face3& e)> func1 = []( const filestream &file, auto &e ) { e.initialize(file); };
			//
			file.read(levelset);
			file.read(velocity);
			file.read(area);
			file.read(solid_cell);
			file.read(ghost_cells,func0);
			file.read(ghost_faces,func1);
			file.r(flux_boundary_condition);
			file.r(cell_count);
			file.r(face_count);
			file.r(ghost_cell_reference_count);
			file.r(ghost_face_reference_count);
			file.read(cell_map);
			file.read(face_map);
			file.r(valid_cell_count);
			file.r(valid_face_count);
		}
		virtual void serialize( const filestream &file ) const override {
			size_t size = layers.size();
			file.w(size);
			for( unsigned n=0; n<size; ++n ) {
				layers[n]->recursive_serialize(file);
			}
			//
			std::function<void(const filestream &file, const ghost_cell3& e)> func0 = []( const filestream &file, const auto &e ) { e.serialize(file); };
			std::function<void(const filestream &file, const ghost_face3& e)> func1 = []( const filestream &file, const auto &e ) { e.serialize(file); };
			//
			file.write(levelset);
			file.write(velocity);
			file.write(area);
			file.write(solid_cell);
			file.write(ghost_cells,func0);
			file.write(ghost_faces,func1);
			file.w(flux_boundary_condition);
			file.w(cell_count);
			file.w(face_count);
			file.w(ghost_cell_reference_count);
			file.w(ghost_face_reference_count);
			file.write(cell_map);
			file.write(face_map);
			file.w(valid_cell_count);
			file.w(valid_face_count);
		}
		//
		struct Parameters {
			unsigned change_sparse_array_resolution {64};
			std::string sparse_grid_module {"tiledarray3"};
			int adaptivity_type {0};
			double padding_for_remeshing {1.0};
			bool steep_adapvitiy {true};
			unsigned dilate_count {3};
			double eps {1e-2};
			double precision_eps {1e-8};
			bool first_order {false};
			bool clamp_order {true};
			double clamp_fluid_eps {0.0};
			double clamp_solid_eps {1e-2};
			double MLS_eps {1e-2};
			bool MLS_constant_diag {false};
			bool accurate_interpolation {false};
			bool use_inverse {false};
			unsigned surftens_smooth_count {0};
			bool simple_redistance {true};
			bool normalize_laplacian {false};
			unsigned pde_update_count {0};
			unsigned erosion_count {0};
			bool debug {false};
		};
		Parameters param;
		//
		vec3d get_cell_position( const cell_id3 &cell_id ) const {
			const auto &layer = layers[cell_id.depth];
			return layer->dx*cell_id.pi.cell();
		}
		//
		vec3d get_face_position( const face_id3 &face_id ) const {
			const auto &layer = layers[face_id.depth];
			return layer->dx*face_id.pi.face(face_id.dim);
		}
		//
		double get_cell_dx( const cell_id3 &cell_id ) const {
			return layers[cell_id.depth]->dx;
		}
		//
		double get_face_dx( const face_id3 &face_id ) const {
			return layers[face_id.depth]->dx;
		}
		//
		uint_type get_cell_index( const cell_id3 &cell_id ) const {
			return layers[cell_id.depth]->active_cells(cell_id.pi);
		}
		//
		uint_type get_face_index( const face_id3 &face_id ) const {
			return layers[face_id.depth]->active_faces[face_id.dim](face_id.pi);
		}
		//
		void clear();
		void add_layer( const shape3 &shape, double dx );
		void activate_cells( std::function<bool(char depth, const vec3d &p)> func );
		void activate_cells( std::function<double(const vec3d &p)> fluid, std::function<double(const vec3d &p)> solid );
		void assign_levelset( std::function<double( const vec3d &p )> fluid, std::function<double( const vec3d &p )> solid );
		void set_velocity( std::function<double( const vec3d &p, char dim )> func );
		void set_flux_boundary_condition( const flux_boundary_condition3 &boundary_cond );
		void balance_layers();
		void assign_indices();
		void iterate_cell_neighbors( const cell_id3 &cell_id, std::function<void( char dim, const cell_id3 &cell_id )> func ) const;
		void iterate_face_neighbors( const face_id3 &face_id, std::function<void( const face_id3 &face_id )> func ) const;
		void iterate_face_neighbors( const cell_id3 &cell_id, char dim, std::function<void( const face_id3 &face_id )> func ) const;
		void iterate_active_cells( const std::function<void( const cell_id3 &cell_id, int thread_index )> func ) const;
		void iterate_active_faces( const std::function<void( const face_id3 &face_id, int thread_index )> func ) const;
		void serial_iterate_active_cells( const std::function<void( const cell_id3 &cell_id )> func ) const;
		void serial_iterate_active_faces( const std::function<void( const face_id3 &face_id )> func ) const;
		//
		struct gradient_info3 {
			Real levelset;
			Real dx;
			Real rho;
			Real area;
			bool t_junction;
			bool cross_interface;
			bool compromised;
		};
		//
		size_t get_compromised_gradient_count() const;
		size_t get_surface_T_junction_count() const;
		void get_gradient( const face_id3 &face_id, std::function<void( const cell_id3 &cell_id, double value, const gradient_info3 &info )> func ) const;
		void get_scaled_gradient( const face_id3 &face_id, std::function<void( const cell_id3 &cell_id, double value, const gradient_info3 &info )> func ) const;
		void get_divergence( const cell_id3 &cell_id, std::function<void( const face_id3 &face_id, double value0, double value1 )> func ) const;
		void get_unmofidied_divergence( const cell_id3 &cell_id, std::function<void( const face_id3 &face_id, double value )> func ) const;
		void add_surfacetense_force( double coeff, double dt );
		void extrapolate( std::function<double(const vec3d &p)> solid_func );
		void extrapolate_toward_solid( std::function<double(const vec3d &p)> solid_func );
		//
		bool find_cell( const vec3d &p, cell_id3 &cell_id ) const;
		bool find_face( const vec3d &p, char dim, face_id3 &face_id ) const;
		//
		double sample_cell( const vec3d &p, const std::vector<Real> &cell_values, Real *min_max_values=nullptr ) const;
		double sample_levelset( const vec3d &p, Real *min_max_values=nullptr ) const;
		double sample_face( const vec3d &p, char dim, const std::vector<Real> &face_values, Real *min_max_values=nullptr ) const;
		double sample_velocity( const vec3d &p, char dim, Real *min_max_values=nullptr ) const;
		vec3d sample_velocity( const vec3d &p ) const;
		double sample_solid_face_velocity( const face_id3 &face_id, std::function<vec3d( const vec3d &p )> func ) const;
		//
		bool try_interp_cell( const vec3d &p, const std::vector<Real> &cell_values, double &value, Real *min_max_values=nullptr ) const;
		bool try_interp_face( const vec3d &p, char dim, const std::vector<Real> &face_values, double &value, Real *min_max_values=nullptr ) const;
		//
		struct point_info3 {
			double dx;
			unsigned char shrink_info {0};
			Real value;
			double weight;
			vec3d position;
			uint_type index;
		};
		//
		void gather_neighbor_faces( const vec3d &p, char dim, std::vector<point_info3> &points, const std::vector<Real> &values ) const;
		void gather_neighbor_cells( const vec3d &p, std::vector<point_info3> &points, const std::vector<Real> &values ) const;
		//
		bool reconstruct_fluid( array3<Real> &fluid, std::function<bool(int i, int j, int k, const Real &levelset_value )> test_func ) const;
		bool reconstruct_velocity( macarray3<Real> &velocity, std::function<bool(int dim, int i, int j, int k)> test_func ) const;
		//
		void compute_cell_map ();
		void compute_face_map ();
		void clear_map ();
		//
		double get_finest_dx() const;
		bool is_surface_cell( const cell_id3 &cell_id ) const;
		bool is_T_junction_cell( const cell_id3 &cell_id ) const;
		double get_volume() const;
		double get_cell_volume( const cell_id3 &cell_id ) const;
		vec3d get_upwind_gradient( const cell_id3 &cell_id, const std::vector<Real> &levelset ) const;
		//
		void draw_grid( graphics_engine &g, double slice_z, bool fill_checkboard=false ) const;
	};
};
//
SHKZ_END_NAMESPACE
//
#endif
//