/*
**	macoctreemesher3.cpp
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
#include "macoctreemesher3.h"
#include <Eigen/Dense>
#include <shiokaze/utility/utility.h>
//
SHKZ_USING_NAMESPACE
using namespace macotreeliquid3_namespace;
//
void macoctreemesher3::configure( configuration &config ) {
	//
	configuration::auto_group group(config,*this);
	//
	config.get_bool("SnapVertices",m_param.snap_vertices,"Whether to snap vertices onto surfaces");
	config.get_bool("DisplaceSolidEmbedded",m_param.displace_solid_embedded_vertices,"Whether to move solid embedded vertices");
	config.get_double("DisplaceSolidEmbeddedRate",m_param.displace_solid_embedded_rate,"Displace rate for solid embedded vertices");
	config.get_bool("ProjectVertices",m_param.project_vertices,"Whether to project vertices on surfaces");
	config.get_double("MeshProximity",m_param.mesh_proximity,"Mesh merging proximity");
	config.get_bool("RemoveFace",m_param.remove_face,"Whether to remove skinny faces");
}
//
void macoctreemesher3::generate_mesh( const grid3 &grid, const std::vector<Real> &cell_values, double offset, std::function<double(const vec3d &p)> solid_func, std::vector<vec3d> &vertices, std::vector<std::vector<size_t> > &faces, bool enclose ) const {
	//
	assert(grid.layers.size());
	shape3 h_shape = grid.layers[0]->shape;
	double h_dx = grid.layers[0]->dx;
	if( ! solid_func ) solid_func = []( const vec3d &p ) { return 1.0; };
	//
	vertices.clear();
	faces.clear();
	//
	struct face_info {
		vec3d pos;
		vec3d normal;
		bool boundary;
	};
	//
	struct dc_layer {
		macarray3<face_info> active_faces;
	};
	//
	std::vector<std::shared_ptr<dc_layer> > dc_layers;
	//
	// Allocate a set of layers
	for( char depth=0; depth<grid.layers.size(); ++depth ) {
		//
		double dx = grid.layers[depth]->dx;
		shape3 shape = grid.layers[depth]->shape;
		auto layer = std::shared_ptr<dc_layer>(new dc_layer());
		//
		layer->active_faces.initialize(shape);
		//
		dc_layers.push_back(layer);
	}
	//
	// Find faces that should be activated
	for( char depth=0; depth<grid.layers.size(); ++depth ) {
		//
		grid.layers[depth]->active_cells.const_serial_actives([&]( int i, int j, int k, const auto &it) {
			auto &layer = *dc_layers[depth];
			for( char dim : DIMS3 ) {
				layer.active_faces[dim].set(i,j,k,{vec3d(),vec3d(),false});
				layer.active_faces[dim].set(i+(dim==0),j+(dim==1),k+(dim==2),{vec3d(),vec3d(),false});
			}
		});
	}
	//
	for( char depth=1; depth<grid.layers.size(); ++depth ) {
		//
		auto &coarse_layer = *dc_layers[depth];
		auto &fine_layer = *dc_layers[depth-1];
		//
		coarse_layer.active_faces.parallel_actives([&]( char dim, int i, int j, int k, auto &it ) {
			if( fine_layer.active_faces[dim].active(2*i,2*j,2*k)) {
				it.set_off();
			}
		});
	}
	//
	for( char depth=0; depth<grid.layers.size(); ++depth ) {
		//
		auto &layer = *dc_layers[depth];
		double dx = grid.layers[depth]->dx;
		layer.active_faces.parallel_actives([&]( char dim, int i, int j, int k, auto &it ) {
			//
			vec3d corners_pos[4];
			for( int jj=0; jj<2; ++jj ) for( int ii=0; ii<2; ++ii ) {
				int idx = ii+2*jj;
				if( dim == 0 ) {
					corners_pos[idx] = dx * vec3d(i,j+ii,k+jj);
				} else if( dim == 1 ) {
					corners_pos[idx] = dx * vec3d(i+ii,j,k+jj);
				} else if( dim == 2 ) {
					corners_pos[idx] = dx * vec3d(i+ii,j+jj,k);
				}
			}
			int inside_count (0);
			for( int n=0; n<4; ++n ) {
				 if( solid_func(corners_pos[n]) < -0.5*dx ) inside_count ++;
			}
			//
			if( inside_count == 4 ) it.set_off();
			else {
				//
				auto test_cell = [&]( const vec3i pi, double &phi, vec3d &pos, bool &boundary_face ) {
					//
					const auto &active_cells = grid.layers[depth]->active_cells;
					if( active_cells.shape().out_of_bounds(pi)) {
						if( enclose ) {
							boundary_face = true;
						} else {
							it.set_off();
						}
					} else {
						//
						if( active_cells.active(pi)) {
							uint_type idx = active_cells(pi);
							phi = cell_values[idx]+offset;
							pos = grid.get_cell_position({depth,pi,idx});
						} else if( depth < grid.layers.size()-1 ) {
							const auto &coarse_active_cells = grid.layers[depth+1]->active_cells;
							if( coarse_active_cells.active(pi/2)) {
								uint_type idx = coarse_active_cells(pi/2);
								phi = cell_values[idx]+offset;
								pos = grid.get_cell_position({(char)(depth+1),pi/2,idx});
							}
						}
					}
				};
				//
				bool boundary_back (false), boundary_front (false);
				double phi_back (1.0), phi_front (1.0);
				vec3d pos_back, pos_front;
				//
				test_cell(vec3i(i-(dim==0),j-(dim==1),k-(dim==2)),phi_back,pos_back,boundary_back);
				test_cell(vec3i(i,j,k),phi_front,pos_front,boundary_front);
				//
				if( it.active()) {
					if( boundary_back || boundary_front ) {
						if( enclose ) {
							if( phi_back < 0.0 || phi_front < 0.0 ) {
								vec3d pos = dx*vec3i(i,j,k).face(dim);
								vec3d normal = boundary_back ? -vec3d(dim==0,dim==1,dim==2) : vec3d(dim==0,dim==1,dim==2);
								it.set({pos,normal,true});
							}
							else it.set_off();
						}
					} else {
						if( phi_back * phi_front > 0.0 ) {
							it.set_off();
						} else {
							// Compute edge vertex
							double fraction = utility::fraction(phi_back,phi_front);
							vec3d p;
							if( phi_back < 0.0 ) {
								p = fraction*pos_front+(1.0-fraction)*pos_back;
							} else {
								p = fraction*pos_back+(1.0-fraction)*pos_front;
							}
							vec3d normal = phi_back < 0.0 ? vec3d(dim==0,dim==1,dim==2) : -vec3d(dim==0,dim==1,dim==2);
							it.set({p,normal,false});
						}
					}
				}
			}
		});
	}
	//
	struct vertex_info {
		size_t index;
		std::vector<vec3d> intersection_positions;
		std::vector<vec3d> face_positions;
	};
	//
	// Find vertices
	shape3 nodal_shape = h_shape+shape3(1,1,1);
	std::map<size_t,vertex_info> map_vertices;
	size_t index (0);
	//
	std::vector<vec3d> raw_vertices;
	std::vector<Real> vertex_dx;
	//
	for( char depth=0; depth<grid.layers.size(); ++depth ) {
		//
		auto &layer = *dc_layers[depth];
		double dx = grid.layers[depth]->dx;
		unsigned ratio = dx / h_dx;

		layer.active_faces.const_serial_actives([&]( char dim, int i, int j, int k, const auto &it ) {
			//
			vec3i corners[4];
			for( int jj=0; jj<2; ++jj ) for( int tmp_ii=0; tmp_ii<2; ++tmp_ii ) {
				int ii = jj == 0 ? tmp_ii : 1-tmp_ii;
				int idx = tmp_ii+2*jj;
				if( dim == 0 ) {
					corners[idx] = ratio * vec3i(i,j+ii,k+jj);
				} else if( dim == 1 ) {
					corners[idx] = ratio * vec3i(i+ii,j,k+jj);
				} else if( dim == 2 ) {
					corners[idx] = ratio * vec3i(i+ii,j+jj,k);
				}
			}
			//
			const double eps = 1e-3;
			vec3d fpos = dx*vec3i(i,j,k).face(dim);
			bool boundary_face = fpos == it().pos;
			//
			for( int n=0; n<4; ++n ) {
				size_t coord = nodal_shape.encode(corners[n]);
				vec3d default_pos = h_dx*vec3d(corners[n]);
				if( map_vertices.find(coord) == map_vertices.end()) {
					vec3d raw_pos = h_dx*vec3d(corners[n]);
					if( boundary_face ) {
						map_vertices[coord] = { index++, {}, {} };
					} else {
						map_vertices[coord] = { index++, {it().pos}, {fpos}};
					}
					vertices.push_back(raw_pos);
					raw_vertices.push_back(raw_pos);
					vertex_dx.push_back(dx);
				} else {
					if( ! boundary_face ) {
						map_vertices[coord].intersection_positions.push_back(it().pos);
						map_vertices[coord].face_positions.push_back(fpos);
					}
				}
			}
			//
			if( depth > 0 ) {
				if( dim == 0 ) {
					corners[0] = vec3i(ratio * vec3d(i,j+0.5,k));
					corners[1] = vec3i(ratio * vec3d(i,j+1,k+0.5));
					corners[2] = vec3i(ratio * vec3d(i,j+0.5,k+1));
					corners[3] = vec3i(ratio * vec3d(i,j,k+0.5));
				} else if( dim == 1 ) {
					corners[0] = vec3i(ratio * vec3d(i+0.5,j,k));
					corners[1] = vec3i(ratio * vec3d(i+1,j,k+0.5));
					corners[2] = vec3i(ratio * vec3d(i+0.5,j,k+1));
					corners[3] = vec3i(ratio * vec3d(i,j,k+0.5));
				} else if( dim == 2 ) {
					corners[0] = vec3i(ratio * vec3d(i+0.5,j,k));
					corners[1] = vec3i(ratio * vec3d(i+1,j+0.5,k));
					corners[2] = vec3i(ratio * vec3d(i+0.5,j+1,k));
					corners[3] = vec3i(ratio * vec3d(i,j+0.5,k));
				}
				//
				for( int n=0; n<4; ++n ) {
					size_t coord = nodal_shape.encode(corners[n]);
					vec3d default_pos = h_dx*vec3d(corners[n]);
					if( map_vertices.find(coord) != map_vertices.end()) {
						if( ! boundary_face ) {
							map_vertices[coord].intersection_positions.push_back(it().pos);
							map_vertices[coord].face_positions.push_back(fpos);
						}
					}
				}
			}
		});
	}
	//
	if( m_param.snap_vertices ) {
		//
		for( auto it=map_vertices.begin(); it!=map_vertices.end(); it++ ) {
			//
			size_t idx = it->second.index;
			double dx = vertex_dx[idx];
			//
			bool snaps[DIM3][2];
			for( char dim : DIMS3 ) {
				snaps[dim][0] = vertices[idx][dim] < 0.5 * h_dx;
				snaps[dim][1] = vertices[idx][dim] > h_dx * (h_shape[dim]-0.5);
			}
			//
			bool failed (false);
			double sum (0.0);
			vec3d average_pos;
			for( const auto &p : it->second.intersection_positions ) {
				double dist2 = (p-raw_vertices[idx]).norm2();
				if( dist2 ) {
					double w = 1.0 / sqrtf(dist2);
					average_pos += w * p;
					sum += w;
				}
			}
			if( sum ) {
				average_pos = average_pos / sum;
			}
			//
			double fsum (0.0);
			vec3d average_fpos;
			for( const auto &p : it->second.face_positions ) {
				average_fpos += p;
				fsum += 1.0;
			}
			//
			if( fsum ) {
				average_fpos = average_fpos / fsum;
				vertices[idx] = average_fpos;
			}
			if( m_param.project_vertices ) {
				try {
					const auto &positions_save = it->second.intersection_positions;
					//
					if( it->second.intersection_positions.size() >= 3 ) {
						//
						const auto &positions_save = it->second.intersection_positions;
						double area (0.0);
						int fi, fj, fk;
						for( int i=0; i<positions_save.size(); ++i ) for( int j=i+1; j<positions_save.size(); ++j ) for( int k=j+1; k<positions_save.size(); ++k ) {
							vec3d r0 = positions_save[j] - positions_save[i];
							vec3d r1 = positions_save[k] - positions_save[i];
							double a2 = (r0 ^ r1).norm2();
							if( area < a2 ) {
								area = a2; fi = i; fj = j; fk = k;
							}
						}
						std::vector<vec3d> positions;
						positions.push_back(positions_save[fi]);
						positions.push_back(positions_save[fj]);
						positions.push_back(positions_save[fk]);
						//
						vec3d origin = vertices[idx];
						//
						Eigen::MatrixXd A(3,3);
						Eigen::VectorXd b(3);
						for( unsigned row=0; row<3; ++row ) {
							A(row,0) = (positions[row][0]-origin[0]) / dx;
							A(row,1) = (positions[row][1]-origin[1]) / dx;
							A(row,2) = (positions[row][2]-origin[2]) / dx;
						}
						for( unsigned row=0; row<3; ++row ) {
							b(row) = 1.0;
						}
						Eigen::VectorXd x;
						if( A.determinant() ) {
							x = A.inverse() * b;
						} else throw 0;
						double norm = x.squaredNorm();
						if( norm ) {
							vertices[idx] += dx * vec3d(x(0),x(1),x(2)) / norm;
						} else throw 2;
					} else throw 3;
				} catch(...) {
					if( fsum ) {
						vertices[idx] = average_pos;
					}
				}
			}
			//
			for( char dim : DIMS3 ) {
				if( snaps[dim][0] ) vertices[idx][dim] = 0.0;
				if( snaps[dim][1] ) vertices[idx][dim] = h_dx * h_shape[dim];
			}
		}
	}
	//
	if( m_param.displace_solid_embedded_vertices ) {
		//
		const double t = m_param.displace_solid_embedded_rate;
		m_parallel.for_each(vertices.size(),[&]( size_t n ) {
			//
			double solid_levelset = solid_func(vertices[n]);
			double gap (0.5*h_dx);
			if( solid_levelset < -gap ) {
				//
				vec3d grad;
				for( int dim : DIMS3 ) {
					grad[dim] = solid_func(vertices[n]+0.5*h_dx*vec3d(dim==0,dim==1,dim==2)) -
								solid_func(vertices[n]-0.5*h_dx*vec3d(dim==0,dim==1,dim==2));
				}
				grad.normalize();
				grad[1] = 0.0;
				vec3d new_position = vertices[n]-(solid_levelset+gap) * grad;
				vertices[n] = (1.0-t)*vertices[n]+t*new_position;
			}
		});
	}
	//
	// Patch faces on them
	std::vector<vec3d> face_normals;
	std::vector<bool> face_flags;
	for( char depth=0; depth<grid.layers.size(); ++depth ) {
		//
		auto &layer = *dc_layers[depth];
		double dx = grid.layers[depth]->dx;
		unsigned ratio = dx / h_dx;

		layer.active_faces.const_serial_actives([&]( char dim, int i, int j, int k, const auto &it ) {
			//
			if( depth == 0 ) {
				vec3i corners[4];
				for( int jj=0; jj<2; ++jj ) for( int tmp_ii=0; tmp_ii<2; ++tmp_ii ) {
					int ii = jj == 0 ? tmp_ii : 1-tmp_ii;
					int idx = tmp_ii+2*jj;
					if( dim == 0 ) {
						corners[idx] = ratio * vec3i(i,j+ii,k+jj);
					} else if( dim == 1 ) {
						corners[idx] = ratio * vec3i(i+ii,j,k+jj);
					} else if( dim == 2 ) {
						corners[idx] = ratio * vec3i(i+ii,j+jj,k);
					}
				}
				//
				std::vector<size_t> face;
				for( int n=0; n<4; ++n ) {
					size_t coord = nodal_shape.encode(corners[n]);
					assert( map_vertices.find(coord) != map_vertices.end());
					face.push_back(map_vertices[coord].index);
				}
				faces.push_back(face);
				face_normals.push_back(it().normal);
				face_flags.push_back(it().boundary);
			} else {
				vec3i corners[8];
				if( dim == 0 ) {
					corners[0] = vec3i(ratio * vec3d(i,j,k));
					corners[1] = vec3i(ratio * vec3d(i,j+0.5,k));
					corners[2] = vec3i(ratio * vec3d(i,j+1,k));
					corners[3] = vec3i(ratio * vec3d(i,j+1,k+0.5));
					corners[4] = vec3i(ratio * vec3d(i,j+1,k+1));
					corners[5] = vec3i(ratio * vec3d(i,j+0.5,k+1));
					corners[6] = vec3i(ratio * vec3d(i,j,k+1));
					corners[7] = vec3i(ratio * vec3d(i,j,k+0.5));
				} else if( dim == 1 ) {
					corners[0] = vec3i(ratio * vec3d(i,j,k));
					corners[1] = vec3i(ratio * vec3d(i+0.5,j,k));
					corners[2] = vec3i(ratio * vec3d(i+1,j,k));
					corners[3] = vec3i(ratio * vec3d(i+1,j,k+0.5));
					corners[4] = vec3i(ratio * vec3d(i+1,j,k+1));
					corners[5] = vec3i(ratio * vec3d(i+0.5,j,k+1));
					corners[6] = vec3i(ratio * vec3d(i,j,k+1));
					corners[7] = vec3i(ratio * vec3d(i,j,k+0.5));
				} else if( dim == 2 ) {
					corners[0] = vec3i(ratio * vec3d(i,j,k));
					corners[1] = vec3i(ratio * vec3d(i+0.5,j,k));
					corners[2] = vec3i(ratio * vec3d(i+1,j,k));
					corners[3] = vec3i(ratio * vec3d(i+1,j+0.5,k));
					corners[4] = vec3i(ratio * vec3d(i+1,j+1,k));
					corners[5] = vec3i(ratio * vec3d(i+0.5,j+1,k));
					corners[6] = vec3i(ratio * vec3d(i,j+1,k));
					corners[7] = vec3i(ratio * vec3d(i,j+0.5,k));
				}
				int count (0);
				for( int n=0; n<8; ++n ) {
					size_t coord = nodal_shape.encode(corners[n]);
					if( map_vertices.find(coord) != map_vertices.end()) {
						++ count;
					}
				}
				if( count == 4 ) {
					std::vector<size_t> face;
					face.push_back(map_vertices[nodal_shape.encode(corners[0])].index);
					face.push_back(map_vertices[nodal_shape.encode(corners[2])].index);
					face.push_back(map_vertices[nodal_shape.encode(corners[4])].index);
					face.push_back(map_vertices[nodal_shape.encode(corners[6])].index);
					faces.push_back(face);
					face_normals.push_back(it().normal);
					face_flags.push_back(it().boundary);
				} else {
					//
					vec3d center = (
						vertices[map_vertices[nodal_shape.encode(corners[0])].index]+
						vertices[map_vertices[nodal_shape.encode(corners[2])].index]+
						vertices[map_vertices[nodal_shape.encode(corners[4])].index]+
						vertices[map_vertices[nodal_shape.encode(corners[6])].index] ) / 4.0;
					vec3d raw_center = dx*vec3i(i,j,k).face(dim);
					//
					vertices.push_back(center);
					raw_vertices.push_back(raw_center);
					vertex_dx.push_back(dx);
					//
					size_t coord0 = nodal_shape.encode(corners[0]);
					for( int n=1; n<9; ++n ) {
						size_t coord1 = nodal_shape.encode(corners[n%8]);
						if( map_vertices.find(coord1) != map_vertices.end()) {
							std::vector<size_t> face;
							face.push_back(index);
							face.push_back(map_vertices[coord0].index);
							face.push_back(map_vertices[coord1].index);
							faces.push_back(face);
							face_normals.push_back(it().normal);
							face_flags.push_back(it().boundary);
							coord0 = coord1;
						}
					}
					index ++;
				}
			}
		});
	}
	//
	// FLIP face normal
	m_parallel.for_each(faces.size(),[&]( size_t n ) {
		//
		vec3d vec0 = raw_vertices[faces[n][1]]-raw_vertices[faces[n][0]];
		vec3d vec1 = raw_vertices[faces[n][2]]-raw_vertices[faces[n][0]];
		vec3d normal = vec1 ^ vec0;
		vec3d gradient = face_normals[n];
		//
		if( gradient * normal < 0.0 ) {
			std::vector<size_t> save_face = faces[n];
			for( int i=0; i<faces[n].size(); ++i ) {
				faces[n][i] = save_face[faces[n].size()-i-1];
			}
		}
	});
	//
	// Remove skinny faces
	if( m_param.remove_face ) {
		//
		std::vector<std::vector<size_t> > connections(vertices.size());
		for( size_t n=0; n<faces.size(); ++n ) {
			if( face_flags[n] ) continue;
			for( size_t i=0; i<faces[n].size(); ++i ) {
				int j = (i+1) % faces[n].size();
				auto &ci = connections[faces[n][i]];
				if( find(ci.begin(),ci.end(),faces[n][j]) == ci.end()) {
					ci.push_back(faces[n][j]);
				}
				auto &cj = connections[faces[n][j]];
				if( find(cj.begin(),cj.end(),faces[n][i]) == cj.end()) {
					cj.push_back(faces[n][i]);
				}
			}
		}
		//
		std::vector<size_t> trimmed_vertices(vertices.size());
		std::vector<bool> trimmed_flags(vertices.size());
		for( size_t n=0; n<vertices.size(); ++n ) {
			trimmed_vertices[n] = n;
			for( size_t j=0; j<connections[n].size(); ++j ) {
				size_t m = connections[n][j];
				if( n < m && ! trimmed_flags[m] ) {
					double dx = std::max(vertex_dx[n],vertex_dx[m]);
					if( (vertices[n]-vertices[m]).len() < m_param.mesh_proximity * dx ) {
						trimmed_vertices[n] = m;
						vertices[m] = 0.5 * (vertices[n]+vertices[m]);
						trimmed_flags[m] = true;
						break;
					}
				}
			}
		}
		//
		size_t removed_count = 0;
		for( size_t n=0; n<trimmed_vertices.size(); ++n ) {
			if( trimmed_vertices[n] != n ) removed_count ++;
		}
		//
		if( removed_count ) {
			m_parallel.for_each(faces.size(),[&]( size_t n ) {
				for( int i=0; i<faces[n].size(); ++i ) {
					faces[n][i] = trimmed_vertices[faces[n][i]];
				}
			});
			//
			std::vector<std::vector<size_t> > face_save (faces);
			faces.clear();
			for( auto f : face_save ) {
				char duplicate_count (0);
				for( int i=0; i<f.size(); ++i ) for( int j=i+1; j<f.size(); ++j ) {
					if( f[i] == f[j] ) duplicate_count ++;
				}
				if( f.size() - duplicate_count > 2 ) {
					faces.push_back(f);
				}
			}
		}
		//
		// Re-order vertices
		std::vector<bool> touched (vertices.size(),false);
		for( size_t n=0; n<faces.size(); ++n ) {
			for( int i=0; i<faces[n].size(); ++i ) {
				touched[faces[n][i]] = true;
			}
		}
		size_t reorder_index (0);
		std::vector<size_t> reoderered_vertices (vertices.size(),0);
		for( size_t n=0; n<vertices.size(); ++n ) {
			if( touched[n] ) reoderered_vertices[n] = ++ reorder_index;
		}
		//
		for( size_t n=0; n<faces.size(); ++n ) {
			for( int i=0; i<faces[n].size(); ++i ) {
				faces[n][i] = reoderered_vertices[faces[n][i]]-1;
			}
		}
		//
		std::vector<vec3d> vertices_save (vertices);
		vertices.clear();
		for( size_t n=0; n<reoderered_vertices.size(); ++n ) {
			if( reoderered_vertices[n] ) {
				vertices.push_back(vertices_save[n]);
			}
		}
	}
}
//