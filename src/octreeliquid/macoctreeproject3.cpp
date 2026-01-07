/*
**	macoctreeproject3.cpp
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
#include "macoctreeproject3.h"
#include <shiokaze/math/RCMatrix_utility.h>
#include <shiokaze/core/console.h>
#include <shiokaze/core/timer.h>
#include <shiokaze/utility/utility.h>
#include <numeric>
#include <cassert>
#include <fstream>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseCholesky>
#include <eigen3/Eigen/IterativeLinearSolvers>
//
SHKZ_USING_NAMESPACE
using namespace macotreeliquid3_namespace;
//
bool macoctreeproject3::solveWithICCG (Eigen::SparseMatrix<double> &A, Eigen::VectorXd &b, Eigen::VectorXd &x) {
	// ラーメンの具材（行列）の準備（Aは正定値対称行列であることを想定）

    // 1. 不完全コレスキー分解
    Eigen::IncompleteCholesky<double> ic;
    ic.compute(A);

    if (ic.info() != Eigen::Success) {
        std::cerr << "前処理のIC分解に失敗...\n";
        return false;
    }

    // 2. 共役勾配法 (CG法) のソルバ
    // 前処理子として、ic (Incomplete Cholesky) を指定！
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<double>> cg;
    //cg.setPreconditioner(ic); // 前処理子を設定
	int max_iterations = 1000;
	double tol = 1e-6;
    cg.setMaxIterations(max_iterations);
    cg.setTolerance(tol);
    cg.compute(A);  // 行列を設定

    // 3. 解を求める！
    x = cg.solve(b);

	if (cg.iterations() < max_iterations) {
		console::dump("  PCG converged in %d iterations with residual %e\n", cg.iterations(), cg.error());
	} else {
		console::dump("  PCG did not converge in %d iterations, residual %e\n", cg.iterations(), cg.error());
	}
	return true;
}
void is_valid_double(double value) {
       if(std::isnan(value) || std::isinf(value)) console::dump("Invalid double value: %f\n", value);
}
void macoctreeproject3::configure( configuration &config ) {
	//
	configuration::auto_group group(config,*this);
	//
	config.get_bool("VolumeCorrection",m_param.volume_correction,"Whether to perform volume correction");
	config.get_bool("VolumeCorrectSkipSurface",m_param.volume_correct_skip_surfaces,"Skip near surfaces for volume correct");
	config.set_default_bool("ForceGlobalResidual",true);
	config.get_bool("CheckSymmetric",m_param.check_symmetric,"Check matrix symmetircity");
	config.get_bool("CheckPositiveDiag",m_param.check_positive_diag,"Check the all the diag elements are positive");
	config.get_bool("FixDivergence",m_param.fix_divergence,"Fix divergence due to the moving solid");
	config.get_bool("DebugAssemble",m_param.debug_assemble,"Debug mode for assembling matrix");
	config.get_double("VolumeRecoverRatio",m_param.volume_recover_ratio,"Volume recover ratio for volume correction");
	config.get_double("MinimalVolumeCorrectTimeStep",m_param.minimal_volume_correct_timestep,"Minimal time step size for volume correction");
	config.get_double("NP",m_param.np,"Constant np for volume correction");
	config.get_bool("RemoveOneDegreesOfFreedom",m_param.remove_one_degrees_of_freedom,"Remove one degrees of freedom");
	config.get_integer("SpecialBoundaryCondition",m_param.special_boundary_condition,"Type of special boundary condition");
}
//
void macoctreeproject3::post_initialize( bool initialized_from_file ) {
	//
	m_initial_volume = 0.0;
	m_y_prev = 0.0;
}
void macoctreeproject3::add_element_symm(uint_type row, uint_type col, double value) {
	if (m_matrix.Lhs->get(col, row)) {
		m_matrix.Lhs->add_to_element(row, col, m_matrix.Lhs->get(col, row));
	} else {
		m_matrix.Lhs->add_to_element(row, col, value);
	}
}
//
void macoctreeproject3::assemble_matrix( grid3 &grid ) {
	//
	scoped_timer timer(this);
	timer.tick(); console::dump( ">>> Assembling linear system...\n" );
	//
	grid.compute_cell_map();
	grid.compute_face_map();
	//
	uint_type cell_count = grid.valid_cell_count;
	uint_type face_count = grid.valid_face_count;
	//
	if( m_matrix.allocated ) {
		m_matrix.Lhs->initialize(cell_count,cell_count);
		if( m_param.debug_assemble ) {
			m_matrix.G->initialize(face_count,cell_count);
			m_matrix.D->initialize(cell_count,face_count);
		}
	} else {
		m_matrix.Lhs = m_factory->allocate_matrix(cell_count,cell_count);
		if( m_param.debug_assemble ) {
			m_matrix.G = m_factory->allocate_matrix(face_count,cell_count);
			m_matrix.D = m_factory->allocate_matrix(cell_count,face_count);
		}
		m_matrix.allocated = true;
	}
	//
	if( m_param.debug_assemble ) {
		//
		// Scaled gradient matrix
		grid.iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
			if( grid.face_map[face_id.index] ) {
				uint_type row = grid.face_map[face_id.index]-1;
				grid.get_scaled_gradient(face_id,[&]( const cell_id3 &cell_id, double value, const grid3::gradient_info3 &info ) {
					if( value ) {
						assert( grid.cell_map[cell_id.index] );
						uint_type column = grid.cell_map[cell_id.index]-1;
						m_matrix.G->add_to_element(row,column,value);
					}
				});
			}
		});
		//
		// Divergence matrix
		grid.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
			if( grid.cell_map[cell_id.index] ) {
				uint_type row = grid.cell_map[cell_id.index]-1;
				grid.get_divergence(cell_id,[&]( const face_id3 &face_id, double value0, double value1 ) {
					if( value0 && grid.face_map[face_id.index] ) {
						uint_type column = grid.face_map[face_id.index]-1;
						m_matrix.D->add_to_element(row,column,value0);
					}
				});
			}
		});
		//
		m_matrix.D->multiply(m_matrix.G.get(),m_matrix.Lhs.get());
		RCMatrix_utility<size_t,double>::report(m_matrix.Lhs.get(),"Lhs");
		//
	} else {
		//
		// Directily assemble Lhs matrix
		grid.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
			if( grid.cell_map[cell_id.index] ) {
				uint_type row0 = grid.cell_map[cell_id.index]-1;
				double diag (0.0);
				int count = 0;
				grid.get_divergence(cell_id,[&]( const face_id3 &face_id, double value0, double value1 ) {
					if( grid.face_map[face_id.index] ) {
						if( value0 ) {
							grid.get_scaled_gradient(face_id,[&]( const cell_id3 &cell_neigh_id, double value, const grid3::gradient_info3 &info ) {
								if( value ) {
									uint_type row1 = grid.cell_map[cell_neigh_id.index]-1;
									const double v = value0*value;

									if( row0 == row1 ) {
										diag += v;
										count++;
										// console::dump("value0: %f, value: %f v: %f row0: %zu, row1: %zu\n", value0, value, v, row0, row1);
									}
									else {
										m_matrix.Lhs->add_to_element(row0,row1,v);
									}
								}
							});
						}
					}
				});
				// if (diag != 0.0) console::dump("Adding Lhs(%zu,%zu) += %f\n", row0, row0, diag);
				m_matrix.Lhs->add_to_element(row0,row0,diag);
				//console::dump("Row %zu: diag=%f, count=%d\n", row0, diag, count);
			}
		});
		RCMatrix_utility<size_t,double>::report(m_matrix.Lhs.get(),"Lhs");
	}

	//
	if( m_param.check_symmetric ) {
		const double symm_error = RCMatrix_utility<size_t,double>::symmetricity_error(m_matrix.Lhs.get());
		assert( symm_error == 0.0 );
	}
	if( m_param.check_positive_diag ) {
		const double min_diag = RCMatrix_utility<size_t,double>::min_diag(m_matrix.Lhs.get());
		assert( min_diag > 0.0 );
	}
	m_matrix.assembled = true;
	//
	// Gather information
	std::vector<unsigned> num_t_junction_bucket(m_parallel.get_thread_num(),0);
	std::vector<unsigned> num_compromised_bucket(m_parallel.get_thread_num(),0);
	grid.iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
		//
		bool t_junction (false);
		bool compromised (false);
		//
		if( grid.face_map[face_id.index] ) {
			uint_type row = grid.face_map[face_id.index]-1;
			grid.get_scaled_gradient(face_id,[&]( const cell_id3 &cell_id, double value, const grid3::gradient_info3 &info ) {
				if( info.t_junction && info.cross_interface ) {
					t_junction = true;
					if( info.compromised ) compromised = true;
				}
			});
			if( t_junction ) {
				num_t_junction_bucket[tid] ++;
				if( compromised ) num_compromised_bucket[tid] ++;
			}
		}
	});
	//
	unsigned num_t_junction = std::accumulate(num_t_junction_bucket.begin(),num_t_junction_bucket.end(),0);
	unsigned num_compromised = std::accumulate(num_compromised_bucket.begin(),num_compromised_bucket.end(),0);
	double compromised_ratio = num_t_junction ? num_compromised/(double)num_t_junction : 0.0;
	console::dump( "num_t_junction = %u, num_compromised=%u (compromised_ratio=%.2f%%)\n", num_t_junction, num_compromised, 100.0*compromised_ratio );
	//
	console::write("num_proj_t_junction",num_t_junction);
	console::write("num_proj_compromised",num_compromised);
	console::write("num_proj_compromised_ratio",compromised_ratio);
	//
	console::dump( "<<< ...Done. Took %s.\n", timer.stock("assemble_matrix").c_str());
}
//圧力点p1,p2,p3,p4のAを構築するためのMatrix上のidと、それぞれのm_grid上のcell_id3を受け取り、LHS行列を構築する
void macoctreeproject3::assemble_Lhs(grid3 &grid,std::vector<int> p, cell_id3 top1, cell_id3 bottom1, cell_id3 top2, cell_id3 bottom2) {
	
	//top/bottomは圧力点のcell_id3
	double upper_virtual_face_x = (grid.upper_face_xposition(top1) < grid.upper_face_xposition(top2)) ? grid.upper_face_xposition(top1) : grid.upper_face_xposition(top2);
	double lower_virtual_face_x = (grid.lower_face_xposition(bottom1) > grid.lower_face_xposition(bottom2)) ? grid.lower_face_xposition(bottom1) : grid.lower_face_xposition(bottom2);
	double center_virtual_face_x = 0.5 * (upper_virtual_face_x + lower_virtual_face_x);
	double distance1 = grid.get_cell_position(top1)[0] - grid.get_cell_position(bottom1)[0];
	double distance2 = grid.get_cell_position(top2)[0] - grid.get_cell_position(bottom2)[0];
	double t = (grid.get_cell_position(top1)[0] - center_virtual_face_x) / distance1;
	double s = (grid.get_cell_position(top2)[0] - center_virtual_face_x) / distance2;
	std::vector<double> g = {1-t, t, -(1-s), -s};
	for (int j = 0; j<4; j++) {
		for (int k = 0; k<4; k++) {
			m_matrix.Lhs->add_to_element(p[j], p[k], g[j]*g[k]);
		}
	}
}

void macoctreeproject3::assemble_lhs_yz(grid3 &grid, std::vector<uint_type> p, std::vector<cell_id3> cell_ids, cell_id3 cell_id_this, char dim, double dx, int dir) {
	//j0<j1, j2<j3と仮定。jaはminimal faceの下、jbは上
	if (cell_ids[0].index != cell_ids[1].index && cell_ids[2].index != cell_ids[3].index) {
		std::vector<int> j(4);
		for (int i = 0; i < 4; i++) {
			j[i] = cell_ids[i].pi[0];
		}
		int J = cell_id_this.pi[0];
		int ja = std::max(j[0], j[2]);
		int jb = std::min(j[1], j[3]);
		//virtual faceの位置
		double scale_to_upper = dx * (double)(J-j[0]) / (j[1]-j[0]);
		double scale_to_lower = dx * (double)(j[1]-J) / (j[1]-j[0]);
		std::vector<double> coeff(4);
		//coeffは論文から符号反転がなされてこの値になる。
		//∇⋅(∇p)=−∇⋅uの右辺のマイナスが邪魔なので両辺にマイナスをかけて左辺のマイナスはAに含めてしまおう！という発想
		coeff[0] = (double) (j[1] - J) / (double)(j[1] - j[0]);
		coeff[1] = (double) (J - j[0]) / (double)(j[1] - j[0]);
		coeff[2] = - (double) (j[3] - J) / (double)(j[3] - j[2]);
		coeff[3] = - (double) (J - j[2]) / (double)(j[3] - j[2]);
		std::vector<uint_type> cell_maps(4);
		for (int i = 0; i < 4; i++) {
			cell_maps[i] = grid.cell_map[cell_ids[i].index];
		}
		//LHSの構築
		//水平方向成分
		for (int i = 0; i < 4; i++) {
			if (cell_maps[i]) {
				if (cell_maps[i] > grid.valid_cell_count) console::dump("cell_map[%zu] = %u\n", cell_ids[i].index, cell_maps[i]);
				assert(cell_maps[i] <= grid.valid_cell_count);
				if (std::isinf(scale_to_upper) || std::isinf(scale_to_lower)) {
					console::dump("1.Assemble_lhs_yz: scale is inf. p[%d]: %zu, i0=%d, i1=%d, i2=%d, i3=%d, J=%d\n", i, p[i], j[0], j[1], j[2], j[3], J);
				}
				double ratio = 1.0 / (jb - ja + 1);
				// console::dump("1.Assemble_lhs_yz: p[1]=%zu, p[0]=%zu, p[%d]=%zu, scale_to_upper=%f, scale_to_lower=%f, coeff=%f\n", p[1], p[0], i, p[i], scale_to_upper, scale_to_lower, coeff[i]);
				m_matrix.Lhs->add_to_element(p[1], p[i], ratio * scale_to_upper *coeff[i]);
				m_matrix.Lhs->add_to_element(p[0], p[i], ratio * scale_to_lower *coeff[i]);
				//仮に定数だったとしたら収束はどうなる？
				// A.coeffRef(p[1], p[i]) += ratio * scale_to_upper *coeff[i];
				// A.coeffRef(p[0], p[i]) += ratio * scale_to_lower *coeff[i];
			}
		}
	}
	else if (cell_ids[0].index != cell_ids[1].index) { //対セルがuniform j1!=j2 j3==j4
		std::vector<int> j(3);
		for (int i = 0; i < 2; i++) {
			j[i] = cell_ids[i].pi[0];
		}
		int J = cell_id_this.pi[0];
		//virtual faceの位置
		double scale_to_upper = dx * (double)(J-j[0]) / (j[1]-j[0]);
		double scale_to_lower = dx * (double)(j[1]-J) / (j[1]-j[0]);
		std::vector<double> coeff(3);
		coeff[0] = (double) (j[1] - J) / (double)(j[1] - j[0]);
		coeff[1] = (double) (J - j[0]) / (double)(j[1] - j[0]);
		coeff[2] = - 1.0;
		//LHSの構築
		//水平方向成分
		for (int i = 0; i < 3; i++) {
				if (std::isinf(scale_to_upper) || std::isinf(scale_to_lower)) {
					console::dump("2.Assemble_lhs_yz: scale is inf. p[%d]: %zu, i0=%d, i1=%d, J=%d\n", i, p[i], j[0], j[1], J);
				}
				//console::dump("2.Assemble_lhs_yz: p[1]=%zu, p[%d]=%zu, scale_to_upper=%f, coeff=%f\n", p[1], i, p[i], scale_to_upper, coeff[i]);
				m_matrix.Lhs->add_to_element(p[1], p[i],  scale_to_upper *coeff[i]);
				m_matrix.Lhs->add_to_element(p[0], p[i],  scale_to_lower *coeff[i]);
				// A.coeffRef(p[1], p[i]) += scale_to_upper *coeff[i];
				// A.coeffRef(p[0], p[i]) += scale_to_lower *coeff[i];
		}
	}
	else if (cell_ids[2].index != cell_ids[3].index) { //自セルがuniform j3!=j4 j1==j2
		std::vector<int> j(3);
		for (int i = 1; i < 4; i++) {
			j[i-1] = cell_ids[i].pi[0];
		}
		int J = cell_id_this.pi[0];
		//virtual faceの位置
		double scale_to_upper = dx;
		double scale_to_lower = dx;
		std::vector<double> coeff(3);
		coeff[0] = 1.0;
		coeff[1] = - (double) (j[2] - J) / (double)(j[2] - j[1]);
		coeff[2] = - (double) (J - j[1]) / (double)(j[2] - j[1]);
		//LHSの構築
		//水平方向成分
		for (int i = 0; i < 3; i++) {
				if (std::isinf(scale_to_upper) || std::isinf(scale_to_lower)) {
					console::dump("3.Assemble_lhs_yz: scale is inf. p[%d]: %zu, i1=%d, i2=%d, J=%d\n", i, p[i+1], j[0], j[1], J);
				}
				// console::dump("3.Assemble_lhs_yz: p[0]=%zu, p[%d]=%zu, scale_to_lower=%f, coeff=%f\n", p[0], i, p[i+1], scale_to_lower, coeff[i]);
				m_matrix.Lhs->add_to_element(p[1], p[i+1],  scale_to_lower *coeff[i]);
				// A.coeffRef(p[1], p[i+1]) += scale_to_lower *coeff[i];
		}
	}
	else { //両方uniform
		double scale = dx;
		//LHSの構築
		//水平方向成分
		if (std::isinf(scale)) {
			console::dump("4.Assemble_lhs_yz: scale is inf. p[0]: %zu\n", p[0]);
		}
		//console::dump("4.Assemble_lhs_yz: p[0]=%zu, p[2]=%zu, scale=%f\n", p[0], p[2], scale);
		m_matrix.Lhs->add_to_element(p[0], p[0], scale);
		m_matrix.Lhs->add_to_element(p[0], p[2], -scale);
		// A.coeffRef(p[0], p[0]) += scale;
		// A.coeffRef(p[0], p[2]) += -scale;
	}
}

void macoctreeproject3::assemble_lhs_x(grid3 &grid, std::vector<uint_type> p, std::vector<cell_id3> cell_ids, cell_id3 cell_id_this, char dim, double dx, int dir) {
	//i1<i2, i3<i4と仮定。iaはminimal faceの左、ibは右
	if (p[0] != p[1]) {
		std::vector<int> i(4);
		for (int k = 0; k < 4; k++) {
			i[k] = cell_ids[k].pi[0];
		}
		int I = cell_id_this.pi[1];
		int ia = std::max(i[0], i[2]);
		int ib = std::min(i[1], i[3]);
		double ratio = 1.0 / (ib - ia + 1);
		//virtual faceの位置
		double scale = dx / (i[1]-i[0]);
		//LHSの構築
		//鉛直方向成分
		if (std::isinf(scale)) {
			console::dump("5.Assemble_lhs_x: scale is inf. pi0=(%d,%d,%d), pi1=(%d,%d,%d), p[0]=%zu, p[1]=%zu, i[0]=%d, i[1]=%d\n", cell_ids[0].pi[0], cell_ids[0].pi[1], cell_ids[0].pi[2], cell_ids[1].pi[0], cell_ids[1].pi[1], cell_ids[1].pi[2], p[0], p[1], i[0], i[1]);
		}
		// if (!std::isinf(scale)) {
			m_matrix.Lhs->add_to_element(p[0], p[0], scale * ratio);
			m_matrix.Lhs->add_to_element(p[0], p[1], -scale * ratio);
			m_matrix.Lhs->add_to_element(p[1], p[0], -scale * ratio);
			m_matrix.Lhs->add_to_element(p[1], p[1], scale * ratio);
			// A.coeffRef(p[0], p[0]) += scale * ratio;
			// A.coeffRef(p[0], p[1]) += -scale * ratio;
			// A.coeffRef(p[1], p[0]) += -scale * ratio;
			// A.coeffRef(p[1], p[1]) += scale * ratio;
		// }
	}
	else {
		double scale = dx;
		//LHSの構築
		//鉛直方向成分
		if (std::isinf(scale)) {
			console::dump("6.Assemble_lhs_x: scale is inf. p[0]: %zu\n", p[0]);
		}
		// if (!std::isinf(scale)) {
			m_matrix.Lhs->add_to_element(p[0], p[0], scale);
			// A.coeffRef(p[0], p[0]) += scale;
		// }
		cell_id3 neighbor_cell_id;
		neighbor_cell_id.depth = cell_id_this.depth;
		neighbor_cell_id.pi = cell_id_this.pi + dir * vec3i(dim==0, dim==1, dim==2);
		neighbor_cell_id.index = grid.layers[neighbor_cell_id.depth]->active_cells(neighbor_cell_id.pi);
		if (grid.cell_map[neighbor_cell_id.index]) {
			uint_type neighbor_p = grid.cell_map[neighbor_cell_id.index]-1;
			m_matrix.Lhs->add_to_element(p[0], neighbor_p, -scale);
			// A.coeffRef(p[0], neighbor_p) += -scale;
		}
	}

}

bool macoctreeproject3::is_in_same_merged_cell( grid3 &grid , cell_id3 cell_id1, cell_id3 cell_id2, UnionFind &uf) {
	//
	return uf.issame(cell_id1.index, cell_id2.index);
}

bool macoctreeproject3::should_Lhs_calc( grid3 &grid , cell_id3 cell_id, UnionFind &uf, char dim) {
	//
	if (cell_id.pi[0] == 0) return true; //左端
	cell_id3 left_cell_id;
	cell_id3 upper_cell_id;
	cell_id3 upper_left_cell_id;
	left_cell_id.depth = cell_id.depth;
	left_cell_id.pi = cell_id.pi + vec3i(-1,0,0);
	left_cell_id.index = grid.layers[left_cell_id.depth]->active_cells(left_cell_id.pi);
	const auto &layer = *grid.layers[cell_id.depth];
	if (layer.active_cells.safe_active(cell_id.pi+vec3i(dim==0,dim==1,dim==2))) {
		upper_cell_id.depth = cell_id.depth;
		upper_cell_id.pi = cell_id.pi + vec3i(dim==0,dim==1,dim==2);
		upper_cell_id.index = grid.layers[upper_cell_id.depth]->active_cells(upper_cell_id.pi);
		upper_left_cell_id.depth = cell_id.depth;
		upper_left_cell_id.pi = cell_id.pi + vec3i(dim==0,dim==1,dim==2) + vec3i(-1,0,0);
		upper_left_cell_id.index = grid.layers[upper_left_cell_id.depth]->active_cells(upper_left_cell_id.pi);
		if (uf.issame(cell_id.index, left_cell_id.index) && uf.issame(upper_cell_id.index, upper_left_cell_id.index)) {
			return false;
		}
		else return true;
	}

	else if (uf.issame(cell_id.index, left_cell_id.index)) {
		return false;
	}
	else return true;
}

//LHS（Ap=bのA）の構築
//cloud_ocに定義されているcell_ids_included_merged_cells（merged_cell, 均一セルを含むセルのcell_id）を用いる
void macoctreeproject3::compute_maps( grid3 &grid , UnionFind &uf) {
	//
	grid.compute_cell_map_merged(uf);
	grid.compute_face_map_merged(uf);
	console::dump("cell_map size: %zu\n", grid.valid_cell_count);
}

void macoctreeproject3::assemble_matrix_merged_cell( grid3 &grid, UnionFind &uf, bool use_Eigen, int all_cell_count, int merged_cell_count, std::vector<cell_id3_and_is_merged> &cell_ids_included_merged_cells, std::vector<int> &merged_cell_id, int step_count) {
	//
	scoped_timer timer(this);
	timer.tick(); console::dump( ">>> Assembling linear system...\n" );
	//
	grid.compute_cell_map_merged(uf);
	grid.compute_face_map_merged(uf);
	console::dump("cell_map size: %zu\n", grid.valid_cell_count);
	//
	uint_type cell_count = grid.valid_cell_count;
	uint_type face_count = grid.valid_face_count;
	//merged_cellの上下2点分を追加
	// A.resize((int)cell_count, (int)cell_count);
	// b.resize((int)cell_count);
	// b.setZero();
	// x.setZero();
	 //to avoid singular matrix when remove one degree of freedom
	//
	if( m_matrix.allocated ) {
		m_matrix.Lhs->initialize(cell_count,cell_count);
		// if( m_param.debug_assemble ) {
		// 	m_matrix.G->initialize(face_count,cell_count);
		// 	m_matrix.D->initialize(cell_count,face_count);
		// }
	} else {
		m_matrix.Lhs = m_factory->allocate_matrix(cell_count,cell_count);
		// if( m_param.debug_assemble ) {
		// 	m_matrix.G = m_factory->allocate_matrix(face_count,cell_count);
		// 	m_matrix.D = m_factory->allocate_matrix(cell_count,face_count);
		// }
		m_matrix.allocated = true;
	}
	// for (int i = 0; i < grid.valid_cell_count; i++) {
	// 	for (int j = 0; j < grid.valid_cell_count; j++) {
	// 		double val1 = m_matrix.Lhs->get(i, j);
	// 		console::dump("Lhs before assembly at (%d, %d): %f\n", i, j, val1);
	// 	}
	// }
	//
	if( m_param.debug_assemble ) {
		//
		// Scaled gradient matrix
		// grid.iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
		// 	if( grid.face_map[face_id.index] ) {
		// 		uint_type row = grid.face_map[face_id.index]-1;
		// 		grid.get_scaled_gradient(face_id,[&]( const cell_id3 &cell_id, double value, const grid3::gradient_info3 &info ) {
		// 			if( value ) {
		// 				assert( grid.cell_map[cell_id.index] );
		// 				uint_type column = grid.cell_map[cell_id.index]-1;
		// 				m_matrix.G->add_to_element(row,column,value);
		// 			}
		// 		});
		// 	}
		// });
		// //
		// // Divergence matrix
		// grid.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
		// 	if( grid.cell_map[cell_id.index] ) {
		// 		uint_type row = grid.cell_map[cell_id.index]-1;
		// 		grid.get_divergence(cell_id,[&]( const face_id3 &face_id, double value0, double value1 ) {
		// 			if( value0 && grid.face_map[face_id.index] ) {
		// 				uint_type column = grid.face_map[face_id.index]-1;
		// 				m_matrix.D->add_to_element(row,column,value0);
		// 			}
		// 		});
		// 	}
		// });
		// //
		// m_matrix.D->multiply(m_matrix.G.get(),m_matrix.Lhs.get());
		// RCMatrix_utility<size_t,double>::report(m_matrix.Lhs.get(),"Lhs");
		//
	} else {
		//
		// Directily assemble Lhs matrix
		//走査するセルを圧力点（均一セルのセル中心とmerged_cellの上下2点）の配列に従って走査する

		grid.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
			if( grid.cell_map[uf.top[uf.root(cell_id.index)].index] || grid.cell_map[uf.bottom[uf.root(cell_id.index)].index] ) {
			for (int dir = -1; dir <= 1; dir+=2) {
				//merged_cellのtop/bottomいずれかに対応する圧力点が存在する場合のみLHSを構築
				for( char dim : DIMS3 ) {
					std::vector<uint_type> p(4); //圧力点のmatrix上のidを格納する配列
					cell_id3 upper_cell_id;
					upper_cell_id.depth = cell_id.depth;
					upper_cell_id.pi = cell_id.pi + dir * vec3i(0, dim==1, dim==2);
					if (!grid.layers[upper_cell_id.depth]->active_cells.safe_active(upper_cell_id.pi)) {
						continue;
					}
					upper_cell_id.index = grid.layers[upper_cell_id.depth]->active_cells(upper_cell_id.pi);
					//p0はbottom1、p1はtop1、p2はbottom2、p3はtop2
					//topはbottom+merged_cell_countで表現
					//ここをcell_mapに依拠した
					
					p[0] = grid.cell_map[uf.bottom[uf.root(cell_id.index)].index]- 1;
					p[1] = grid.cell_map[uf.top[uf.root(cell_id.index)].index]- 1;
					p[2] = grid.cell_map[uf.bottom[uf.root(upper_cell_id.index)].index]- 1;
					p[3] = grid.cell_map[uf.top[uf.root(upper_cell_id.index)].index]- 1;
					
					std::vector<cell_id3> cell_ids(4);
					cell_ids[0] = uf.bottom[uf.root(cell_id.index)];
					cell_ids[1] = uf.top[uf.root(cell_id.index)];
					cell_ids[2] = uf.bottom[uf.root(upper_cell_id.index)];
					cell_ids[3] = uf.top[uf.root(upper_cell_id.index)];
					
					const auto &layer = *grid.layers[cell_id.depth];
					if (layer.active_cells.safe_active(upper_cell_id.pi) && (grid.cell_map[uf.top[uf.root(upper_cell_id.index)].index]) && (grid.cell_map[uf.bottom[uf.root(upper_cell_id.index)].index])) {
						if (dim != 0) {
							assemble_lhs_yz(grid, p, cell_ids, cell_id, dim, grid.layers[0]->dx, dir);
						}
						//p[0]はtop1、p[1]はbottom1、p[2]はtop2、p[3]はbottom2それぞれの(cell_map上の値)-1。これはそのままLHSの行列idとして使える
						else {
							assemble_lhs_x(grid, p, cell_ids, cell_id, dim, grid.layers[0]->dx, dir);
						}
					}
				}
			}
		}
		});
	}
	// 	grid.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
	// 		if( grid.cell_map[cell_id.index] ) {
	// 			uint_type row0 = grid.cell_map[cell_id.index]-1;
	// 			double diag (0.0);
	// 			int count = 0;
	// 			grid.get_divergence(cell_id,[&]( const face_id3 &face_id, double value0, double value1 ) {
	// 				if( grid.face_map[face_id.index] ) {
	// 					if( value0 ) {
	// 						grid.get_scaled_gradient(face_id,[&]( const cell_id3 &cell_neigh_id, double value, const grid3::gradient_info3 &info ) {
	// 							if( value ) {
	// 								uint_type row1 = grid.cell_map[cell_neigh_id.index]-1;
	// 								const double v = value0*value;

	// 								if( row0 == row1 ) {
	// 									diag += v;
	// 									count++;
	// 									console::dump("value0: %f, value: %f v: %f row0: %zu, row1: %zu\n", value0, value, v, row0, row1);
	// 								}
	// 								else {
	// 									m_matrix.Lhs->add_to_element(row0,row1,v);
	// 								}
	// 							}
	// 						});
	// 					}
	// 				}
	// 			});
	// 			// if (diag != 0.0) console::dump("Adding Lhs(%zu,%zu) += %f\n", row0, row0, diag);
	// 			m_matrix.Lhs->add_to_element(row0,row0,diag);
	// 			//console::dump("Row %zu: diag=%f, count=%d\n", row0, diag, count);
	// 		}
	// 	});
		RCMatrix_utility<size_t,double>::report(m_matrix.Lhs.get(),"Lhs");
		
		// 行列をtxt形式で出力
// 行列をtxt形式で出力（Step 40でのみ）
	{
		
		if(step_count == 40) {  // Step 40でのみ出力
			char filename[256];
			
			sprintf(filename, "/home/takebe/shiokaze/project_merged_cell/matrix_lhs_step_40.txt");
			std::ofstream file(filename);
			file << "Row Col Value\n";
			for( size_t i=0; i<m_matrix.Lhs->rows(); ++i ) {
				m_matrix.Lhs->const_for_each(i,[&]( size_t j, double v ) {
					file << i << " " << j << " " << v << "\n";
				});
			}
			file.close();
			console::dump("Matrix dumped to: %s (Step %d)\n", filename, step_count);
		}
	}
	// }
	//

	//ここは通常のLhsを作る際にいくつかの条件（対称性・正の対角成分など）をチェックして想定外の計算ミスをチェックするためのコード
	// for (int i = 0; i < grid.valid_cell_count; i++) {
	// 	for (int j = 0; j < grid.valid_cell_count; j++) {
	// 		double val1 = m_matrix->get(i, j);
	// 		double val2 = m_matrix.Lhs->get(j, i);
	// 		if (val1 != val2) {
	// 			console::dump("Lhs is not symmetric at (%d, %d): %f != %f, difference %f\n", i, j, val1, val2, val1 - val2);
	// 		}
	// 	}
	// }
	// auto Lhs_copy = m_matrix.Lhs;
	// m_matrix.Lhs->clear();
	// //再度Lhs_copyを参照して値を取り出すテスト
	double max_element = 0.0;
	double min_element = 1000000.0;
	for (int i = 0; i < grid.valid_cell_count; i++) {
		for (int j = i; j < grid.valid_cell_count; j++) {
			int number_in_row = 0;
			double val1 = m_matrix.Lhs->get(i, j);
			double val2 = m_matrix.Lhs->get(j, i);

			if (val1 != val2) {
				m_matrix.Lhs->add_to_element(j, i, -val2);
				m_matrix.Lhs->add_to_element(j, i, val1);
				// console::dump("After re-adding, Lhs is not symmetric at (%d, %d): %f != %f, difference %f\n", i, j, val1, val2, val1 - val2);
			}
			if (val1) {
				number_in_row++;
			}
			// if (i == j) {
			// 	m_matrix.Lhs->add_to_element(i, j, -val1 + 750.0); //対角成分はそのまま
			// }
			// if (i != j && val1 && val2) {
			// 	m_matrix.Lhs->add_to_element(i, j, -val1 -187.5);
			// 	m_matrix.Lhs->add_to_element(j, i, -val1 -187.5);
			// }

			// if (val1 && val2 && ((double)std::abs(i-j)/grid.valid_cell_count > 0.1) ) {
			// 	m_matrix.Lhs->add_to_element(i, j, -val1);
			// 	m_matrix.Lhs->add_to_element(j, i, -val1);
			// }
		}
	}
	    // ★ 対角成分から条件数を推定（簡易版）
    std::vector<double> diag(grid.valid_cell_count);
    double max_diag = 0.0;
    double min_diag = INFINITY;
    for(int i = 0; i < grid.valid_cell_count; i++) {
        diag[i] = m_matrix.Lhs->get(i, i);
        max_diag = std::max(max_diag, diag[i]);
        min_diag = std::min(min_diag, diag[i]);
    }
    double cond_number_estimate = max_diag / min_diag;
    double estimated_iterations = std::sqrt(cond_number_estimate);

	console::dump("Matrix diagonal analysis:\n");
    console::dump("  max_diag = %.2e, min_diag = %.2e\n", max_diag, min_diag);
    console::dump("  Estimated condition number κ ≈ %.2e\n", cond_number_estimate);
	console::dump("  Estimated PCG iterations ≈ %.0f\n", estimated_iterations);

	if( m_param.check_symmetric ) {
		const double symm_error = RCMatrix_utility<size_t,double>::symmetricity_error(m_matrix.Lhs.get());
		assert( symm_error == 0.0 );
	}
	if( m_param.check_positive_diag ) {
		const double min_diag = RCMatrix_utility<size_t,double>::min_diag(m_matrix.Lhs.get());
		assert( min_diag > 0.0 );
	}
	m_matrix.assembled = true;
	//
	//以下はt_junction等の記録を書き出すもの
	// Gather information
	// std::vector<unsigned> num_t_junction_bucket(m_parallel.get_thread_num(),0);
	// std::vector<unsigned> num_compromised_bucket(m_parallel.get_thread_num(),0);
	// grid.iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
	// 	//
	// 	bool t_junction (false);
	// 	bool compromised (false);
	// 	//
	// 	if( grid.face_map[face_id.index] ) {
	// 		uint_type row = grid.face_map[face_id.index]-1;
	// 		grid.get_scaled_gradient(face_id,[&]( const cell_id3 &cell_id, double value, const grid3::gradient_info3 &info ) {
	// 			if( info.t_junction && info.cross_interface ) {
	// 				t_junction = true;
	// 				if( info.compromised ) compromised = true;
	// 			}
	// 		});
	// 		if( t_junction ) {
	// 			num_t_junction_bucket[tid] ++;
	// 			if( compromised ) num_compromised_bucket[tid] ++;
	// 		}
	// 	}
	// });
	// //
	// unsigned num_t_junction = std::accumulate(num_t_junction_bucket.begin(),num_t_junction_bucket.end(),0);
	// unsigned num_compromised = std::accumulate(num_compromised_bucket.begin(),num_compromised_bucket.end(),0);
	// double compromised_ratio = num_t_junction ? num_compromised/(double)num_t_junction : 0.0;
	// console::dump( "num_t_junction = %u, num_compromised=%u (compromised_ratio=%.2f%%)\n", num_t_junction, num_compromised, 100.0*compromised_ratio );
	// //
	// console::write("num_proj_t_junction",num_t_junction);
	// console::write("num_proj_compromised",num_compromised);
	// console::write("num_proj_compromised_ratio",compromised_ratio);
	//
	console::dump( "<<< ...Done. Took %s.\n", timer.stock("assemble_matrix").c_str());
}
void macoctreeproject3::assemble_matrix_density( grid3 &grid ) {
	//
	scoped_timer timer(this);
	timer.tick(); console::dump( ">>> Assembling linear system...\n" );
	//
	grid.compute_cell_map();
	grid.compute_face_map();
	//
	uint_type cell_count = grid.valid_cell_count;
	uint_type face_count = grid.valid_face_count;
	//
	if( m_matrix.allocated ) {
		m_matrix.Lhs->initialize(cell_count,cell_count);
		if( m_param.debug_assemble ) {
			m_matrix.G->initialize(face_count,cell_count);
			m_matrix.D->initialize(cell_count,face_count);
		}
	} else {
		m_matrix.Lhs = m_factory->allocate_matrix(cell_count,cell_count);
		if( m_param.debug_assemble ) {
			m_matrix.G = m_factory->allocate_matrix(face_count,cell_count);
			m_matrix.D = m_factory->allocate_matrix(cell_count,face_count);
		}
		m_matrix.allocated = true;
	}

	//
	if( m_param.debug_assemble ) {
		//
		// Scaled gradient matrix
		grid.iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
			if( grid.face_map[face_id.index] ) {
				uint_type row = grid.face_map[face_id.index]-1;
				grid.get_scaled_gradient_density(face_id,[&]( const cell_id3 &cell_id, double value, const grid3::gradient_info3_scalar &info ) {
					if( value ) {
						assert( grid.cell_map[cell_id.index] );
						uint_type column = grid.cell_map[cell_id.index]-1;
						m_matrix.G->add_to_element(row,column,value);
					}
				});
			}
		});
		//
		// Divergence matrix
		grid.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
			if( grid.cell_map[cell_id.index] ) {
				uint_type row = grid.cell_map[cell_id.index]-1;
				grid.get_divergence(cell_id,[&]( const face_id3 &face_id, double value0, double value1 ) {
					if( value0 && grid.face_map[face_id.index] ) {
						uint_type column = grid.face_map[face_id.index]-1;
						m_matrix.D->add_to_element(row,column,value0);
					}
				});
			}
		});
		//
		m_matrix.D->multiply(m_matrix.G.get(),m_matrix.Lhs.get());
		RCMatrix_utility<size_t,double>::report(m_matrix.Lhs.get(),"Lhs");
		//
	} else {
		//
		// Directily assemble Lhs matrix
		grid.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
			if( grid.cell_map[cell_id.index] ) {
				uint_type row0 = grid.cell_map[cell_id.index]-1;
				double diag (0.0);
				grid.get_divergence(cell_id,[&]( const face_id3 &face_id, double value0, double value1 ) {
					if( grid.face_map[face_id.index] ) {
						if( value0 ) {
							grid.get_scaled_gradient(face_id,[&]( const cell_id3 &cell_neigh_id, double value, const grid3::gradient_info3 &info ) {
								if( value ) {
									uint_type row1 = grid.cell_map[cell_neigh_id.index]-1;
									const double v = value0*value;
									if( row0 == row1 ) diag += v;
									else m_matrix.Lhs->add_to_element(row0,row1,v);
								}
							});
						}
					}
				});
				m_matrix.Lhs->add_to_element(row0,row0,diag);
			}
		});
		RCMatrix_utility<size_t,double>::report(m_matrix.Lhs.get(),"Lhs");
	}
	//
	if( m_param.check_symmetric ) {
		const double symm_error = RCMatrix_utility<size_t,double>::symmetricity_error(m_matrix.Lhs.get());
		assert( symm_error == 0.0 );
	}
	if( m_param.check_positive_diag ) {
		const double min_diag = RCMatrix_utility<size_t,double>::min_diag(m_matrix.Lhs.get());
		assert( min_diag > 0.0 );
	}
	m_matrix.assembled = true;
	//
	// Gather information
	std::vector<unsigned> num_t_junction_bucket(m_parallel.get_thread_num(),0);
	std::vector<unsigned> num_compromised_bucket(m_parallel.get_thread_num(),0);
	grid.iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
		//
		bool t_junction (false);
		bool compromised (false);
		//
		if( grid.face_map[face_id.index] ) {
			uint_type row = grid.face_map[face_id.index]-1;
			grid.get_scaled_gradient(face_id,[&]( const cell_id3 &cell_id, double value, const grid3::gradient_info3 &info ) {
				if( info.t_junction && info.cross_interface ) {
					t_junction = true;
					if( info.compromised ) compromised = true;
				}
			});
			if( t_junction ) {
				num_t_junction_bucket[tid] ++;
				if( compromised ) num_compromised_bucket[tid] ++;
			}
		}
	});
	//
	unsigned num_t_junction = std::accumulate(num_t_junction_bucket.begin(),num_t_junction_bucket.end(),0);
	unsigned num_compromised = std::accumulate(num_compromised_bucket.begin(),num_compromised_bucket.end(),0);
	double compromised_ratio = num_t_junction ? num_compromised/(double)num_t_junction : 0.0;
	console::dump( "num_t_junction = %u, num_compromised=%u (compromised_ratio=%.2f%%)\n", num_t_junction, num_compromised, 100.0*compromised_ratio );
	//
	console::write("num_proj_t_junction",num_t_junction);
	console::write("num_proj_compromised",num_compromised);
	console::write("num_proj_compromised_ratio",compromised_ratio);
	//
	console::dump( "<<< ...Done. Took %s.\n", timer.stock("assemble_matrix").c_str());
}

void macoctreeproject3::assemble_matrix_qc( grid3 &grid , int step_count) {
	//
	scoped_timer timer(this);
	timer.tick(); console::dump( ">>> Assembling linear system...\n" );
	//
	grid.compute_cell_map_cloud();
	grid.compute_face_map_cloud();
	//
	uint_type cell_count = grid.valid_cell_count;
	console::dump("Valid cell count: %u\n", cell_count);
	uint_type face_count = grid.valid_face_count;
	//
	if( m_matrix.allocated ) {//行列がもう存在する場合
		m_matrix.Lhs->initialize(cell_count,cell_count);
		if( m_param.debug_assemble ) {
			m_matrix.G->initialize(face_count,cell_count);
			m_matrix.D->initialize(cell_count,face_count);
		}
	} else {//行列がまだない場合
		m_matrix.Lhs = m_factory->allocate_matrix(cell_count,cell_count);
		if( m_param.debug_assemble ) {
			m_matrix.G = m_factory->allocate_matrix(face_count,cell_count);
			m_matrix.D = m_factory->allocate_matrix(cell_count,face_count);
		}
		m_matrix.allocated = true;
	}
	//
	if( m_param.debug_assemble ) {
		//
		// Scaled gradient matrix
		console::dump("Assembling G matrix...\n");
		grid.iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
			if( grid.face_map[face_id.index] ) {
				uint_type row = grid.face_map[face_id.index]-1;
				grid.get_scaled_gradient_qc(face_id,[&]( const cell_id3 &cell_id, double value, const grid3::gradient_info3_cloud &info ) {
					if( value ) {
						assert( grid.cell_map[cell_id.index] );
						uint_type column = grid.cell_map[cell_id.index]-1;
						m_matrix.G->add_to_element(row,column,value);
					}
				});
			}
		});
		//
		// Divergence matrix
		console::dump("Assembling D matrix...\n");
		grid.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
			if( grid.cell_map[cell_id.index] ) {
				uint_type row = grid.cell_map[cell_id.index]-1;
				grid.get_divergence_qc(cell_id,[&]( const face_id3 &face_id, double value0, double value1 ) {
					if( value0 && grid.face_map[face_id.index] ) {
						uint_type column = grid.face_map[face_id.index]-1;
						m_matrix.D->add_to_element(row,column,value0);
					}
				});
			}
		});
		//
		m_matrix.D->multiply(m_matrix.G.get(),m_matrix.Lhs.get());
		RCMatrix_utility<size_t,double>::report(m_matrix.Lhs.get(),"Lhs");
		//
	} else {
		//
		// Directily assemble Lhs matrix
		console::dump("Assembling Lhs matrix...\n");
		grid.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
			if( grid.cell_map[cell_id.index] ) {
				uint_type row0 = grid.cell_map[cell_id.index]-1;
				//console::dump("Assembling row %zu...\n", row0);
				double diag (0.0);
				grid.get_divergence_qc(cell_id,[&]( const face_id3 &face_id, double value0, double value1 ) {
					if( grid.face_map[face_id.index] ) {
						if( value0 ) {
							grid.get_scaled_gradient_qc(face_id,[&]( const cell_id3 &cell_neigh_id, double value, const grid3::gradient_info3_cloud &info ) {
								if( value ) {
									uint_type row1 = grid.cell_map[cell_neigh_id.index]-1;
									const double v = value0*value;
									// if ( row0 == row1 && value0 != 0.0 && value != 0.0) {
									// 	console::dump("Divergence value: %f, Gradient value: %f\n", value0, value);
									// }
									//value0はt*area
									//console::dump("Adding Lhs(%zu,%zu) += %f\n", row0, row1, v);
									if( row0 == row1 ){diag += v; 
										// console::dump("value0: %f, value: %f v: %f, row0: %zu, row1: %zu\n", value0, value, v, row0, row1);
									}
									else {m_matrix.Lhs->add_to_element(row0,row1,v); 
										// console::dump("Adding Lhs(%zu,%zu) += %f\n", row0, row1, v);
									}
								}
							});
						}
					}
				});
				// if (diag != 0.0) console::dump("Adding Lhs(%zu,%zu) += %f\n", row0, row0, diag);
				m_matrix.Lhs->add_to_element(row0,row0,diag);
			}
		});
		RCMatrix_utility<size_t,double>::report(m_matrix.Lhs.get(),"Lhs");
	}
	
		{
		
		if(step_count == 40) {  // Step 40でのみ出力
			char filename[256];
			sprintf(filename, "/home/takebe/matrix_dumps/project/matrix_lhs_step_40.txt");
			std::ofstream file(filename);
			file << "Row Col Value\n";
			for( size_t i=0; i<m_matrix.Lhs->rows(); ++i ) {
				m_matrix.Lhs->const_for_each(i,[&]( size_t j, double v ) {
					file << i << " " << j << " " << v << "\n";
				});
			}
			file.close();
			console::dump("Matrix dumped to: %s (Step %d)\n", filename, step_count);
		}
	}
	//
	if( m_param.check_symmetric ) {
		const double symm_error = RCMatrix_utility<size_t,double>::symmetricity_error(m_matrix.Lhs.get());
		assert( symm_error == 0.0 );
	}
	if( m_param.check_positive_diag ) {
		const double min_diag = RCMatrix_utility<size_t,double>::min_diag(m_matrix.Lhs.get());
		console::dump("Min diag: %f\n", min_diag);
		assert( min_diag > 0.0 );
	}
	m_matrix.assembled = true;
	//
	// Gather information
	std::vector<unsigned> num_t_junction_bucket(m_parallel.get_thread_num(),0);
	std::vector<unsigned> num_compromised_bucket(m_parallel.get_thread_num(),0);
	grid.iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
		//
		bool t_junction (false);
		bool compromised (false);
		//
		if( grid.face_map[face_id.index] ) {
			uint_type row = grid.face_map[face_id.index]-1;
			grid.get_scaled_gradient_qc(face_id,[&]( const cell_id3 &cell_id, double value, const grid3::gradient_info3_cloud &info ) {
				if( info.t_junction && info.cross_interface ) {
					t_junction = true;
					if( info.compromised ) compromised = true;
				}
			});
			if( t_junction ) {
				num_t_junction_bucket[tid] ++;
				if( compromised ) num_compromised_bucket[tid] ++;
			}
		}
	});
	//
	unsigned num_t_junction = std::accumulate(num_t_junction_bucket.begin(),num_t_junction_bucket.end(),0);
	unsigned num_compromised = std::accumulate(num_compromised_bucket.begin(),num_compromised_bucket.end(),0);
	double compromised_ratio = num_t_junction ? num_compromised/(double)num_t_junction : 0.0;
	console::dump( "num_t_junction = %u, num_compromised=%u (compromised_ratio=%.2f%%)\n", num_t_junction, num_compromised, 100.0*compromised_ratio );
	//
	console::write("num_proj_t_junction",num_t_junction);
	console::write("num_proj_compromised",num_compromised);
	console::write("num_proj_compromised_ratio",compromised_ratio);
	//
	console::dump( "<<< ...Done. Took %s.\n", timer.stock("assemble_matrix").c_str());
}


void macoctreeproject3::assemble_matrix_qv( grid3 &grid ) {
	//
	scoped_timer timer(this);
	timer.tick(); console::dump( ">>> Assembling linear system...\n" );
	//
	grid.compute_cell_map_cloud();
	grid.compute_face_map_cloud();
	//
	uint_type cell_count = grid.valid_cell_count;
	uint_type face_count = grid.valid_face_count;
	//
	if( m_matrix.allocated ) {
		m_matrix.Lhs->initialize(cell_count,cell_count);
		if( m_param.debug_assemble ) {
			m_matrix.G->initialize(face_count,cell_count);
			m_matrix.D->initialize(cell_count,face_count);
		}
	} else {
		m_matrix.Lhs = m_factory->allocate_matrix(cell_count,cell_count);
		if( m_param.debug_assemble ) {
			m_matrix.G = m_factory->allocate_matrix(face_count,cell_count);
			m_matrix.D = m_factory->allocate_matrix(cell_count,face_count);
		}
		m_matrix.allocated = true;
	}
	//
	if( m_param.debug_assemble ) {
		//
		// Scaled gradient matrix
		grid.iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
			if( grid.face_map[face_id.index] ) {
				uint_type row = grid.face_map[face_id.index]-1;
				grid.get_scaled_gradient_qv(face_id,[&]( const cell_id3 &cell_id, double value, const grid3::gradient_info3_cloud &info ) {
					if( value ) {
						assert( grid.cell_map[cell_id.index] );
						uint_type column = grid.cell_map[cell_id.index]-1;
						m_matrix.G->add_to_element(row,column,value);
					}
				});
			}
		});
		//
		// Divergence matrix
		grid.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
			if( grid.cell_map[cell_id.index] ) {
				uint_type row = grid.cell_map[cell_id.index]-1;
				grid.get_divergence_qv(cell_id,[&]( const face_id3 &face_id, double value0, double value1 ) {
					if( value0 && grid.face_map[face_id.index] ) {
						uint_type column = grid.face_map[face_id.index]-1;
						m_matrix.D->add_to_element(row,column,value0);
					}
				});
			}
		});
		//
		m_matrix.D->multiply(m_matrix.G.get(),m_matrix.Lhs.get());
		RCMatrix_utility<size_t,double>::report(m_matrix.Lhs.get(),"Lhs");
		//
	} else {
		//
		// Directily assemble Lhs matrix
		grid.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
			if( grid.cell_map[cell_id.index] ) {
				uint_type row0 = grid.cell_map[cell_id.index]-1;
				double diag (0.0);
				grid.get_divergence_qv(cell_id,[&]( const face_id3 &face_id, double value0, double value1 ) {
					if( grid.face_map[face_id.index] ) {
						if( value0 ) {
							grid.get_scaled_gradient_qv(face_id,[&]( const cell_id3 &cell_neigh_id, double value, const grid3::gradient_info3_cloud &info ) {
								if( value ) {
									uint_type row1 = grid.cell_map[cell_neigh_id.index]-1;
									const double v = value0*value;
									if( row0 == row1 ) diag += v;
									else m_matrix.Lhs->add_to_element(row0,row1,v);
								}
							});
						}
					}
				});
				m_matrix.Lhs->add_to_element(row0,row0,diag);
			}
		});
		RCMatrix_utility<size_t,double>::report(m_matrix.Lhs.get(),"Lhs");
	}
	//
	if( m_param.check_symmetric ) {
		const double symm_error = RCMatrix_utility<size_t,double>::symmetricity_error(m_matrix.Lhs.get());
		assert( symm_error == 0.0 );
	}
	if( m_param.check_positive_diag ) {
		const double min_diag = RCMatrix_utility<size_t,double>::min_diag(m_matrix.Lhs.get());
		assert( min_diag > 0.0 );
	}
	m_matrix.assembled = true;
	//
	// Gather information
	std::vector<unsigned> num_t_junction_bucket(m_parallel.get_thread_num(),0);
	std::vector<unsigned> num_compromised_bucket(m_parallel.get_thread_num(),0);
	grid.iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
		//
		bool t_junction (false);
		bool compromised (false);
		//
		if( grid.face_map[face_id.index] ) {
			uint_type row = grid.face_map[face_id.index]-1;
			grid.get_scaled_gradient_qv(face_id,[&]( const cell_id3 &cell_id, double value, const grid3::gradient_info3_cloud &info ) {
				if( info.t_junction && info.cross_interface ) {
					t_junction = true;
					if( info.compromised ) compromised = true;
				}
			});
			if( t_junction ) {
				num_t_junction_bucket[tid] ++;
				if( compromised ) num_compromised_bucket[tid] ++;
			}
		}
	});
	//
	unsigned num_t_junction = std::accumulate(num_t_junction_bucket.begin(),num_t_junction_bucket.end(),0);
	unsigned num_compromised = std::accumulate(num_compromised_bucket.begin(),num_compromised_bucket.end(),0);
	double compromised_ratio = num_t_junction ? num_compromised/(double)num_t_junction : 0.0;
	console::dump( "num_t_junction = %u, num_compromised=%u (compromised_ratio=%.2f%%)\n", num_t_junction, num_compromised, 100.0*compromised_ratio );
	//
	console::write("num_proj_t_junction",num_t_junction);
	console::write("num_proj_compromised",num_compromised);
	console::write("num_proj_compromised_ratio",compromised_ratio);
	//
	console::dump( "<<< ...Done. Took %s.\n", timer.stock("assemble_matrix").c_str());
}

void macoctreeproject3::assemble_matrix_qr( grid3 &grid ) {
	//
	scoped_timer timer(this);
	timer.tick(); console::dump( ">>> Assembling linear system...\n" );
	//
	grid.compute_cell_map_cloud();
	grid.compute_face_map_cloud();
	//
	uint_type cell_count = grid.valid_cell_count;
	uint_type face_count = grid.valid_face_count;
	//
	if( m_matrix.allocated ) {
		m_matrix.Lhs->initialize(cell_count,cell_count);
		if( m_param.debug_assemble ) {
			m_matrix.G->initialize(face_count,cell_count);
			m_matrix.D->initialize(cell_count,face_count);
		}
	} else {
		m_matrix.Lhs = m_factory->allocate_matrix(cell_count,cell_count);
		if( m_param.debug_assemble ) {
			m_matrix.G = m_factory->allocate_matrix(face_count,cell_count);
			m_matrix.D = m_factory->allocate_matrix(cell_count,face_count);
		}
		m_matrix.allocated = true;
	}
	//
	if( m_param.debug_assemble ) {
		//
		// Scaled gradient matrix
		grid.iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
			if( grid.face_map[face_id.index] ) {
				uint_type row = grid.face_map[face_id.index]-1;
				grid.get_scaled_gradient_qr(face_id,[&]( const cell_id3 &cell_id, double value, const grid3::gradient_info3_cloud &info ) {
					if( value ) {
						assert( grid.cell_map[cell_id.index] );
						uint_type column = grid.cell_map[cell_id.index]-1;
						m_matrix.G->add_to_element(row,column,value);
					}
				});
			}
		});
		//
		// Divergence matrix
		grid.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
			if( grid.cell_map[cell_id.index] ) {
				uint_type row = grid.cell_map[cell_id.index]-1;
				grid.get_divergence_qr(cell_id,[&]( const face_id3 &face_id, double value0, double value1 ) {
					if( value0 && grid.face_map[face_id.index] ) {
						uint_type column = grid.face_map[face_id.index]-1;
						m_matrix.D->add_to_element(row,column,value0);
					}
				});
			}
		});
		//
		m_matrix.D->multiply(m_matrix.G.get(),m_matrix.Lhs.get());
		RCMatrix_utility<size_t,double>::report(m_matrix.Lhs.get(),"Lhs");
		//
	} else {
		//
		// Directily assemble Lhs matrix
		grid.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
			if( grid.cell_map[cell_id.index] ) {
				uint_type row0 = grid.cell_map[cell_id.index]-1;
				double diag (0.0);
				grid.get_divergence_qr(cell_id,[&]( const face_id3 &face_id, double value0, double value1 ) {
					if( grid.face_map[face_id.index] ) {
						if( value0 ) {
							grid.get_scaled_gradient_qr(face_id,[&]( const cell_id3 &cell_neigh_id, double value, const grid3::gradient_info3_cloud &info ) {
								if( value ) {
									uint_type row1 = grid.cell_map[cell_neigh_id.index]-1;
									const double v = value0*value;
									if( row0 == row1 ) diag += v;
									else m_matrix.Lhs->add_to_element(row0,row1,v);
								}
							});
						}
					}
				});
				m_matrix.Lhs->add_to_element(row0,row0,diag);
			}
		});
		RCMatrix_utility<size_t,double>::report(m_matrix.Lhs.get(),"Lhs");
	}
	//
	if( m_param.check_symmetric ) {
		const double symm_error = RCMatrix_utility<size_t,double>::symmetricity_error(m_matrix.Lhs.get());
		assert( symm_error == 0.0 );
	}
	if( m_param.check_positive_diag ) {
		const double min_diag = RCMatrix_utility<size_t,double>::min_diag(m_matrix.Lhs.get());
		assert( min_diag > 0.0 );
	}
	m_matrix.assembled = true;
	//
	// Gather information
	std::vector<unsigned> num_t_junction_bucket(m_parallel.get_thread_num(),0);
	std::vector<unsigned> num_compromised_bucket(m_parallel.get_thread_num(),0);
	grid.iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
		//
		bool t_junction (false);
		bool compromised (false);
		//
		if( grid.face_map[face_id.index] ) {
			uint_type row = grid.face_map[face_id.index]-1;
			grid.get_scaled_gradient_qr(face_id,[&]( const cell_id3 &cell_id, double value, const grid3::gradient_info3_cloud &info ) {
				if( info.t_junction && info.cross_interface ) {
					t_junction = true;
					if( info.compromised ) compromised = true;
				}
			});
			if( t_junction ) {
				num_t_junction_bucket[tid] ++;
				if( compromised ) num_compromised_bucket[tid] ++;
			}
		}
	});
	//
	unsigned num_t_junction = std::accumulate(num_t_junction_bucket.begin(),num_t_junction_bucket.end(),0);
	unsigned num_compromised = std::accumulate(num_compromised_bucket.begin(),num_compromised_bucket.end(),0);
	double compromised_ratio = num_t_junction ? num_compromised/(double)num_t_junction : 0.0;
	console::dump( "num_t_junction = %u, num_compromised=%u (compromised_ratio=%.2f%%)\n", num_t_junction, num_compromised, 100.0*compromised_ratio );
	//
	console::write("num_proj_t_junction",num_t_junction);
	console::write("num_proj_compromised",num_compromised);
	console::write("num_proj_compromised_ratio",compromised_ratio);
	//
	console::dump( "<<< ...Done. Took %s.\n", timer.stock("assemble_matrix").c_str());
}
//
void macoctreeproject3::clear_matrix() {
	//
	if( m_matrix.allocated ) {
		m_matrix.Lhs.reset();
		if( m_param.debug_assemble ) {
			m_matrix.G.reset();
			m_matrix.D.reset();
		}
		m_matrix.allocated = false;
	}
	m_matrix.assembled = false;
}
//
void macoctreeproject3::set_moving_solid( std::function<double(const vec3d &p)> moving_solid_func ) {
	m_moving_solid_func = moving_solid_func;
}
//
void macoctreeproject3::project( grid3 &grid, double dt, std::function<vec3d(const vec3d &p)> solid_velocity, std::vector<Real> *pressure ) {
	//
	std::vector<uint_type> regions;
	std::vector<Real> current_volumes;
	std::vector<Real> target_volumes;
	std::vector<Real> y_list;
	//
	project(grid,dt,0,regions,current_volumes,target_volumes,y_list,solid_velocity,pressure);
}
//
void macoctreeproject3::project( grid3 &grid, double dt,
					  size_t region_count,
					  const std::vector<uint_type> &regions,
					  const std::vector<Real> &current_volumes,
					  const std::vector<Real> &target_volumes,
					  std::vector<Real> &y_list,
					  std::function<vec3d(const vec3d &p)> solid_velocity,
					  std::vector<Real> *pressure_vector ) {
	//
	scoped_timer timer(this);
	timer.tick(); console::dump( ">>> Projection step...\n" );
	//
	double rhs_correct (0.0);
	std::vector<Real> rhs_corrects(region_count);
	//
	if( m_param.volume_correction ) {
		const double dt0 = std::max(m_param.minimal_volume_correct_timestep,dt);
		if( region_count ) {
			timer.tick(); console::dump( "Computing regional volume corrector..." );
			for( unsigned n=0; n<region_count; ++n ) {
				double x = (current_volumes[n]-target_volumes[n])/target_volumes[n];
				double y = y_list[n] + x*dt0; y_list[n] = y;
				double kp = -log(1.0-m_param.volume_recover_ratio)/(m_param.np*dt0);
				double ki = kp*kp/16.0;
				rhs_corrects[n] = -(kp*x+ki*y)/(x+1.0);

				assert( ! utility::is_nan(rhs_corrects[n]));
				if( rhs_corrects[n] ) rhs_correct = 1.0;
			}
			console::dump( "Done. Took %s.\n", timer.stock("regional_volume_compute").c_str());
		} else {
			//
			timer.tick(); console::dump( "Computing global volume corrector..." );
			double current_volume = grid.get_volume();//ここの値がlevelsetを参照しているためNanになっている
			if( ! m_initial_volume ) m_initial_volume = current_volume;
			double x = (current_volume-m_initial_volume)/m_initial_volume;
			double y = m_y_prev + x*dt0; m_y_prev = y;
			double kp = -log(1.0-m_param.volume_recover_ratio)/(m_param.np*dt0);
			double ki = kp*kp/16.0;
			rhs_correct = -(kp*x+ki*y)/(x+1.0);
				// assert( ! utility::is_nan(x));
				// assert( ! utility::is_nan(y));
				// assert( ! utility::is_nan(kp));
				// assert( ! utility::is_nan(ki));
			// console::dump("x=%e, y=%e, kp=%e, ki=%e\n",x,y,kp,ki);
			assert( ! utility::is_nan(rhs_correct));
			//
			const double volume_ratio = current_volume / m_initial_volume;
			console::dump("rhs_correct=%e\n",rhs_correct);
			console::write("current_volume",current_volume);
			console::write("volume_ratio",volume_ratio);
			console::dump( "Done. Took %s. Volume ratio=%.2e\n", volume_ratio, timer.stock("volume_compute").c_str());
		}
	}
	//
	auto is_surface_cell = [&]( const cell_id3 &cell_id ) {
		//
		if( m_param.volume_correct_skip_surfaces ) {
			bool touching_air (false);
			grid.iterate_cell_neighbors(cell_id,[&]( char dim, const cell_id3 &cell_n_id ) {
				if( grid.levelset[cell_n_id.index] > 0.0 ) touching_air = true;
			});
			return touching_air;
		} else {
			return false;
		}
	};
	//
	assert( m_matrix.assembled );
	//
	timer.tick(); console::dump( "Allocating vectors..." );
	auto rhs = m_factory->allocate_vector(grid.valid_cell_count);
	auto result = m_factory->allocate_vector(grid.valid_cell_count);

	console::dump( "Done. Took %s.\n", timer.stock("allocate_vector").c_str());
	//
	timer.tick(); console::dump( "Computing divergence..." );
	//
	grid.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
		if( grid.cell_map[cell_id.index] ) {
			uint_type row = grid.cell_map[cell_id.index]-1;
			if( rhs_correct ) {
				bool solid_cell (false);
				for( char dim : DIMS3 ) {
					grid.iterate_face_neighbors(cell_id,dim,[&]( const face_id3 &face_id ){
						if( grid.area[face_id.index] < 1.0 ) solid_cell = true;
					});
					if( solid_cell ) break;
				}
				if( rhs_correct && ! solid_cell && ! is_surface_cell(cell_id)) {
					const double dx = grid.get_cell_dx(cell_id);
					if( region_count ) {
						rhs->add(row,(dx*dx*dx)*rhs_corrects[regions[cell_id.index]-1]);
					} else {
						rhs->add(row,(dx*dx*dx)*rhs_correct);
					}
				}
			}
			//
			if( m_param.fix_divergence ) {
				double err (0.0);
				if( solid_velocity ) {
					grid.get_unmofidied_divergence(cell_id,[&]( const face_id3 &face_id, double value ) {
						err += value*grid.sample_solid_face_velocity(face_id,solid_velocity);
					});
				}
				bool boundary_flag (false);
				for( int dim : DIMS3 ) {
					boundary_flag = cell_id.pi[dim] == 0 || cell_id.pi[dim] == grid.layers[cell_id.depth]->shape[dim]-1;
					if( boundary_flag ) break;
				}
				if( ! boundary_flag ) rhs->add(row,-err);
			}
			//
			const vec3d p = grid.get_cell_position(cell_id);
			const double dx = grid.get_cell_dx(cell_id);
			//
			grid.get_divergence(cell_id,[&]( const face_id3 &face_id, double value0, double value1 ) {
				bool add_divergence (true);
				if( m_param.special_boundary_condition == 1 && cell_id.pi[0] == grid.layers[cell_id.depth]->shape[0]-1 && face_id.dim == 0 ) {
					add_divergence = false;
				}
				if( m_param.special_boundary_condition == 1 && cell_id.pi[0] == 0 && face_id.dim == 0 ) {
					if( m_moving_solid_func(p) < dx ) {
						add_divergence = false;
					}
				}
				if( add_divergence ) {
					if( solid_velocity && value1 ) {
						if( cell_id.pi[0] != grid.layers[cell_id.depth]->shape[0]-1 ) {
							rhs->add(row,value1*grid.sample_solid_face_velocity(face_id,solid_velocity));
						}
					}
					if( value0 ) {
						rhs->add(row,value0*grid.velocity[face_id.index]);
					}
				}
			});
			//
			if( grid.flux_boundary_condition.has_flux()) {
				const vec3d p = grid.get_cell_position(cell_id);
				for( int dim : DIMS3 ) {
					if( cell_id.pi[dim] == 0 ) {
						rhs->add(row,(dx*dx)*grid.flux_boundary_condition.velocity[dim][0]);
					} else if( cell_id.pi[dim] == grid.layers[cell_id.depth]->shape[dim]-1 ) {
						rhs->add(row,-(dx*dx)*grid.flux_boundary_condition.velocity[dim][1]);
					}
				}
			}
		}
	});
	console::dump( "Done. Took %s.\n", timer.stock("divergence_compute").c_str());
	//
	auto compute_vector_kind = [&]( std::vector<unsigned char> &result ) {
		unsigned total_kinds (0);
		result.resize(grid.valid_cell_count);
		grid.serial_iterate_active_cells([&]( const cell_id3 &cell_id ) {
			if( grid.cell_map[cell_id.index] ) {
				uint_type row = grid.cell_map[cell_id.index]-1;
				result[row] = cell_id.depth;
				total_kinds = std::max(total_kinds,(unsigned)result[row]);
			}
		});
		return total_kinds;
	};
	//
	timer.tick(); console::dump( "Counting and forming the vector norm kind..." );
	std::vector<unsigned char> vector_kind; unsigned total_kinds = compute_vector_kind(vector_kind);
	m_solver->register_vector_norm_kind(vector_kind);
	console::dump( "Done. Total=%u. Took %s.\n", total_kinds, timer.stock("count_vector_norm_kind").c_str());
	//
	timer.tick(); console::dump( "Solving the linear system..." );
	rhs->const_for_each([]( size_t row, double value ){
		assert( ! utility::is_nan(value));
	});
	//
	if( m_param.remove_one_degrees_of_freedom ) {
		m_matrix.Lhs->clear(0);
		m_matrix.Lhs->add_to_element(0,0,1.0);
		rhs->set(0,0);
		m_parallel.for_each(grid.valid_cell_count,[&]( size_t row ) {
			m_matrix.Lhs->for_each(row,[&]( size_t column, double& value ) {
				if( column == 0 ) value = 0.0;
			});
		});
	}
	auto status = m_solver->solve(m_matrix.Lhs.get(),rhs.get(),result.get());
	//
	auto build_residual_string = [&]( const std::vector<double> &rhs ) {
		std::string str;
		if( rhs.size() <= 4 ) {
			for( const auto &e : rhs ) {
				str += " " + console::format_str("%.2e",e);
			}
		} else {
			double max_resid (0.0);
			double min_resid (-1.0);
			unsigned slot_max (0), slot_min(0);
			for( int n=0; n<rhs.size(); ++n ) {
				const auto &e = rhs[n];
				if( e > max_resid ) {
					max_resid = e;
					slot_max = n;
				}
				if( min_resid < 0.0 ) min_resid = max_resid;
				else if( e < min_resid ) {
					min_resid = e;
					slot_min = n;
				}
			}
			str = console::format_str("(num=%u, min[%u]=%.2e, max[%u]=%.2e)",status.vector_reresid.size(),min_resid,slot_min,max_resid,slot_max);
		}
		return str;
	};
	for( int n=0; n<status.vector_reresid.size(); ++n ) {
		console::write(console::format_str("reresid%d",n),status.vector_reresid[n]);
		console::write(console::format_str("absresid%d",n),status.vector_absresid[n]);
	}
	console::dump( "Done. Reresid=%.2e. Iterations = %u. Took %s.\n", status.reresid, status.count, timer.stock("solve_linear_system").c_str());
	console::write("num_iteration_count",status.count);
	//
	if( ! status.vector_reresid.empty()) {
		std::string reresid_str = build_residual_string(status.vector_reresid);
		std::string absresid_str = build_residual_string(status.vector_absresid);
		console::dump( "Reresid=%s\n", reresid_str.c_str());
		console::dump( "Absresid=%s\n", absresid_str.c_str());
	}
	//
	timer.tick(); console::dump( "Updating velocity..." );
	for( uint_type n=0; n<grid.face_count; ++n ) {
		if( ! grid.face_map[n] ) grid.velocity[n] = 0.0;
	}
	grid.iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
		grid.get_scaled_gradient(face_id,[&]( const cell_id3 &cell_id, double value, const grid3::gradient_info3 &info ) {
			if( grid.face_map[face_id.index] && grid.cell_map[cell_id.index] && value ) {
				grid.velocity[face_id.index] -= value * result->at(grid.cell_map[cell_id.index]-1);
			}
		});
	});
	//
	if( pressure_vector ) {
		pressure_vector->resize(grid.cell_count);
		for( uint_type n=0; n<grid.cell_count; ++n ) {
			if( grid.cell_map[n] ) {
				(*pressure_vector)[n] = result->at(grid.cell_map[n]-1);
			} else {
				(*pressure_vector)[n] = 0.0;
			}
		}
	}
	//
	console::dump( "Done. Took %s.\n", timer.stock("update_velocity").c_str());
	//
	grid.clear_map();
	console::dump( "<<< Done. Took %s.\n", timer.stock("project").c_str());
}

void macoctreeproject3::project_simple( grid3 &grid, double dt, std::function<vec3d(const vec3d &p)> solid_velocity, std::vector<Real> *pressure ) {
	//
	std::vector<uint_type> regions;
	std::vector<Real> current_volumes;
	std::vector<Real> target_volumes;
	std::vector<Real> y_list;
	//
	project_simple(grid,dt,0,regions,current_volumes,target_volumes,y_list,solid_velocity,pressure);
}

void macoctreeproject3::project_simple( grid3 &grid, double dt,
					  size_t region_count,
					  const std::vector<uint_type> &regions,
					  const std::vector<Real> &current_volumes,
					  const std::vector<Real> &target_volumes,
					  std::vector<Real> &y_list,
					  std::function<vec3d(const vec3d &p)> solid_velocity,
					  std::vector<Real> *pressure_vector ) {
	//
	scoped_timer timer(this);
	timer.tick(); console::dump( ">>> Projection step...\n" );
	//
	double rhs_correct (0.0);
	std::vector<Real> rhs_corrects(region_count);
	//
	// if( m_param.volume_correction ) {
	// 	const double dt0 = std::max(m_param.minimal_volume_correct_timestep,dt);
	// 	if( region_count ) {
	// 		timer.tick(); console::dump( "Computing regional volume corrector..." );
	// 		for( unsigned n=0; n<region_count; ++n ) {
	// 			double x = (current_volumes[n]-target_volumes[n])/target_volumes[n];
	// 			double y = y_list[n] + x*dt0; y_list[n] = y;
	// 			double kp = -log(1.0-m_param.volume_recover_ratio)/(m_param.np*dt0);
	// 			double ki = kp*kp/16.0;
	// 			rhs_corrects[n] = -(kp*x+ki*y)/(x+1.0);

	// 			assert( ! utility::is_nan(rhs_corrects[n]));
	// 			if( rhs_corrects[n] ) rhs_correct = 1.0;
	// 		}
	// 		console::dump( "Done. Took %s.\n", timer.stock("regional_volume_compute").c_str());
	// 	} else {
	// 		//
	// 		timer.tick(); console::dump( "Computing global volume corrector..." );
	// 		double current_volume = grid.get_volume();//ここの値がlevelsetを参照しているためNanになっている
	// 		if( ! m_initial_volume ) m_initial_volume = current_volume;
	// 		double x = (current_volume-m_initial_volume)/m_initial_volume;
	// 		double y = m_y_prev + x*dt0; m_y_prev = y;
	// 		double kp = -log(1.0-m_param.volume_recover_ratio)/(m_param.np*dt0);
	// 		double ki = kp*kp/16.0;
	// 		rhs_correct = -(kp*x+ki*y)/(x+1.0);
	// 			// assert( ! utility::is_nan(x));
	// 			// assert( ! utility::is_nan(y));
	// 			// assert( ! utility::is_nan(kp));
	// 			// assert( ! utility::is_nan(ki));
	// 		console::dump("x=%e, y=%e, kp=%e, ki=%e\n",x,y,kp,ki);
	// 		assert( ! utility::is_nan(rhs_correct));
	// 		//
	// 		const double volume_ratio = current_volume / m_initial_volume;
	// 		console::write("current_volume",current_volume);
	// 		console::write("volume_ratio",volume_ratio);
	// 		console::dump( "Done. Took %s. Volume ratio=%.2e\n", volume_ratio, timer.stock("volume_compute").c_str());
	// 	}
	// }
	//
	auto is_surface_cell = [&]( const cell_id3 &cell_id ) {
		//
		if( m_param.volume_correct_skip_surfaces ) {
			bool touching_air (false);
			grid.iterate_cell_neighbors(cell_id,[&]( char dim, const cell_id3 &cell_n_id ) {
				if( grid.levelset[cell_n_id.index] > 0.0 ) touching_air = true;
			});
			return touching_air;
		} else {
			return false;
		}
	};
	//
	assert( m_matrix.assembled );
	//
	timer.tick(); console::dump( "Allocating vectors..." );
	auto rhs = m_factory->allocate_vector(grid.valid_cell_count);
	auto result = m_factory->allocate_vector(grid.valid_cell_count);
	console::dump( "Done. Took %s.\n", timer.stock("allocate_vector").c_str());
	//
	timer.tick(); console::dump( "Computing divergence..." );
	//
	grid.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
		if( grid.cell_map[cell_id.index] ) {
			uint_type row = grid.cell_map[cell_id.index]-1;
			// if( rhs_correct ) {
			// 	bool solid_cell (false);
			// 	for( char dim : DIMS3 ) {
			// 		grid.iterate_face_neighbors(cell_id,dim,[&]( const face_id3 &face_id ){
			// 			if( grid.area[face_id.index] < 1.0 ) solid_cell = true;
			// 		});
			// 		if( solid_cell ) break;
			// 	}
			// 	if( rhs_correct && ! solid_cell && ! is_surface_cell(cell_id)) {
			// 		const double dx = grid.get_cell_dx(cell_id);
			// 		if( region_count ) {
			// 			rhs->add(row,(dx*dx*dx)*rhs_corrects[regions[cell_id.index]-1]);
			// 		} else {
			// 			rhs->add(row,(dx*dx*dx)*rhs_correct);
			// 		}
			// 	}
			// }
			//rhs補正の部分を削除
			//
			if( m_param.fix_divergence ) {
				double err (0.0);
				if( solid_velocity ) {
					grid.get_unmofidied_divergence(cell_id,[&]( const face_id3 &face_id, double value ) {
						err += value*grid.sample_solid_face_velocity(face_id,solid_velocity);
					});
				}
				bool boundary_flag (false);
				for( int dim : DIMS3 ) {
					boundary_flag = cell_id.pi[dim] == 0 || cell_id.pi[dim] == grid.layers[cell_id.depth]->shape[dim]-1;
					if( boundary_flag ) break;
				}
				if( ! boundary_flag ) rhs->add(row,-err);
			}
			//
			const vec3d p = grid.get_cell_position(cell_id);
			const double dx = grid.get_cell_dx(cell_id);
			//
			grid.get_divergence(cell_id,[&]( const face_id3 &face_id, double value0, double value1 ) {
				bool add_divergence (true);
				if( m_param.special_boundary_condition == 1 && cell_id.pi[0] == grid.layers[cell_id.depth]->shape[0]-1 && face_id.dim == 0 ) {
					add_divergence = false;
				}
				if( m_param.special_boundary_condition == 1 && cell_id.pi[0] == 0 && face_id.dim == 0 ) {
					if( m_moving_solid_func(p) < dx ) {
						add_divergence = false;
					}
				}
				if( add_divergence ) {
					if( solid_velocity && value1 ) {
						//area=0(つまり固体のときはこちらが呼ばれ、固体界面の速度で色々計算する)
						if( cell_id.pi[0] != grid.layers[cell_id.depth]->shape[0]-1 ) {
							rhs->add(row,value1*grid.sample_solid_face_velocity(face_id,solid_velocity));
						}
					}
					if( value0 ) {
						//多分大体こっちが通常計算される
						rhs->add(row,value0*grid.velocity[face_id.index]);
					}
				}
			});
			//
			// if( grid.flux_boundary_condition.has_flux()) {
			// 	const vec3d p = grid.get_cell_position(cell_id);
			// 	for( int dim : DIMS3 ) {
			// 		if( cell_id.pi[dim] == 0 ) {
			// 			rhs->add(row,(dx*dx)*grid.flux_boundary_condition.velocity[dim][0]);
			// 		} else if( cell_id.pi[dim] == grid.layers[cell_id.depth]->shape[dim]-1 ) {
			// 			rhs->add(row,-(dx*dx)*grid.flux_boundary_condition.velocity[dim][1]);
			// 		}
			// 	}
			// }
		}
	});
	console::dump( "Done. Took %s.\n", timer.stock("divergence_compute").c_str());
	//
	auto compute_vector_kind = [&]( std::vector<unsigned char> &result ) {
		unsigned total_kinds (0);
		result.resize(grid.valid_cell_count);
		grid.serial_iterate_active_cells([&]( const cell_id3 &cell_id ) {
			if( grid.cell_map[cell_id.index] ) {
				uint_type row = grid.cell_map[cell_id.index]-1;
				result[row] = cell_id.depth;
				total_kinds = std::max(total_kinds,(unsigned)result[row]);
			}
		});
		return total_kinds;
	};
	//
	timer.tick(); console::dump( "Counting and forming the vector norm kind..." );
	std::vector<unsigned char> vector_kind; unsigned total_kinds = compute_vector_kind(vector_kind);
	m_solver->register_vector_norm_kind(vector_kind);
	console::dump( "Done. Total=%u. Took %s.\n", total_kinds, timer.stock("count_vector_norm_kind").c_str());
	//
	timer.tick(); console::dump( "Solving the linear system..." );
	rhs->const_for_each([]( size_t row, double value ){
		assert( ! utility::is_nan(value));
	});
	//
	if( m_param.remove_one_degrees_of_freedom ) {
		m_matrix.Lhs->clear(0);
		m_matrix.Lhs->add_to_element(0,0,1.0);
		rhs->set(0,0);
		m_parallel.for_each(grid.valid_cell_count,[&]( size_t row ) {
			m_matrix.Lhs->for_each(row,[&]( size_t column, double& value ) {
				if( column == 0 ) value = 0.0;
			});
		});
	}
	auto status = m_solver->solve(m_matrix.Lhs.get(),rhs.get(),result.get());
	//
	auto build_residual_string = [&]( const std::vector<double> &rhs ) {
		std::string str;
		if( rhs.size() <= 4 ) {
			for( const auto &e : rhs ) {
				str += " " + console::format_str("%.2e",e);
			}
		} else {
			double max_resid (0.0);
			double min_resid (-1.0);
			unsigned slot_max (0), slot_min(0);
			for( int n=0; n<rhs.size(); ++n ) {
				const auto &e = rhs[n];
				if( e > max_resid ) {
					max_resid = e;
					slot_max = n;
				}
				if( min_resid < 0.0 ) min_resid = max_resid;
				else if( e < min_resid ) {
					min_resid = e;
					slot_min = n;
				}
			}
			str = console::format_str("(num=%u, min[%u]=%.2e, max[%u]=%.2e)",status.vector_reresid.size(),min_resid,slot_min,max_resid,slot_max);
		}
		return str;
	};
	for( int n=0; n<status.vector_reresid.size(); ++n ) {
		console::write(console::format_str("reresid%d",n),status.vector_reresid[n]);
		console::write(console::format_str("absresid%d",n),status.vector_absresid[n]);
	}
	console::dump( "Done. Reresid=%.2e. Iterations = %u. Took %s.\n", status.reresid, status.count, timer.stock("solve_linear_system").c_str());
	console::write("num_iteration_count",status.count);
	//
	if( ! status.vector_reresid.empty()) {
		std::string reresid_str = build_residual_string(status.vector_reresid);
		std::string absresid_str = build_residual_string(status.vector_absresid);
		console::dump( "Reresid=%s\n", reresid_str.c_str());
		console::dump( "Absresid=%s\n", absresid_str.c_str());
	}
	//
	timer.tick(); console::dump( "Updating velocity..." );
	for( uint_type n=0; n<grid.face_count; ++n ) {
		if( ! grid.face_map[n] ) grid.velocity[n] = 0.0;
	}
	//速度場の更新
	grid.iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
		grid.get_scaled_gradient_qc(face_id,[&]( const cell_id3 &cell_id, double value, const grid3::gradient_info3_cloud &info ) {
			if( grid.face_map[face_id.index] && grid.cell_map[cell_id.index] /*&& value */ ) {
				grid.velocity[face_id.index] -= value * result->at(grid.cell_map[cell_id.index]-1);
				//grid.velocity[face_id.index] -= result->at(grid.cell_map[cell_id.index]-1);
			}//もしかしたらvalueの代わりにdt/rhoをかけるかもしれない
		});
	});
	//
	//圧力が入力としてある場合はresultにその値を格納
	if( pressure_vector ) {
		pressure_vector->resize(grid.cell_count);
		for( uint_type n=0; n<grid.cell_count; ++n ) {
			if( grid.cell_map[n] ) {
				(*pressure_vector)[n] = result->at(grid.cell_map[n]-1);
			} else {
				(*pressure_vector)[n] = 0.0;
			}
		}
	}
	//
	console::dump( "Done. Took %s.\n", timer.stock("update_velocity").c_str());
	//
	grid.clear_map();
	console::dump( "<<< Done. Took %s.\n", timer.stock("project").c_str());
}
//


void macoctreeproject3::project_density( grid3 &grid, double dt, std::function<vec3d(const vec3d &p)> solid_velocity, std::vector<Real> *pressure ) {
	//
	std::vector<uint_type> regions;
	std::vector<Real> current_volumes;
	std::vector<Real> target_volumes;
	std::vector<Real> y_list;
	//
	project_density(grid,dt,0,regions,current_volumes,target_volumes,y_list,solid_velocity,pressure);
}
//
void macoctreeproject3::project_density( grid3 &grid, double dt,
					  size_t region_count,
					  const std::vector<uint_type> &regions,
					  const std::vector<Real> &current_volumes,
					  const std::vector<Real> &target_volumes,
					  std::vector<Real> &y_list,
					  std::function<vec3d(const vec3d &p)> solid_velocity,
					  std::vector<Real> *pressure_vector ) {
	//
	scoped_timer timer(this);
	timer.tick(); console::dump( ">>> Projection step...\n" );
	//
	double rhs_correct (0.0);
	std::vector<Real> rhs_corrects(region_count);
	//
	if( m_param.volume_correction ) {
		const double dt0 = std::max(m_param.minimal_volume_correct_timestep,dt);
		if( region_count ) {
			timer.tick(); console::dump( "Computing regional volume corrector..." );
			for( unsigned n=0; n<region_count; ++n ) {
				double x = (current_volumes[n]-target_volumes[n])/target_volumes[n];
				double y = y_list[n] + x*dt0; y_list[n] = y;
				double kp = -log(1.0-m_param.volume_recover_ratio)/(m_param.np*dt0);
				double ki = kp*kp/16.0;
				rhs_corrects[n] = -(kp*x+ki*y)/(x+1.0);

				assert( ! utility::is_nan(rhs_corrects[n]));
				if( rhs_corrects[n] ) rhs_correct = 1.0;
			}
			console::dump( "Done. Took %s.\n", timer.stock("regional_volume_compute").c_str());
		} else {
			//
			timer.tick(); console::dump( "Computing global volume corrector..." );
			double current_volume = grid.get_volume_density();//ここの値がlevelsetを参照しているためNanになっている
			console::dump( "current volume = %e", current_volume);
			if( ! m_initial_volume ) m_initial_volume = current_volume;
			double x = (current_volume-m_initial_volume)/m_initial_volume;
			double y = m_y_prev + x*dt0; m_y_prev = y;
			double kp = -log(1.0-m_param.volume_recover_ratio)/(m_param.np*dt0);
			double ki = kp*kp/16.0;
			rhs_correct = -(kp*x+ki*y)/(x+1.0);
				// assert( ! utility::is_nan(x));
				// assert( ! utility::is_nan(y));
				// assert( ! utility::is_nan(kp));
				// assert( ! utility::is_nan(ki));
			// console::dump("x=%e, y=%e, kp=%e, ki=%e\n",x,y,kp,ki);
			assert( ! utility::is_nan(rhs_correct));
			//
			const double volume_ratio = current_volume / m_initial_volume;
			console::write("current_volume",current_volume);
			console::write("volume_ratio",volume_ratio);
			console::dump( "Done. Took %s. Volume ratio=%.2e\n", volume_ratio, timer.stock("volume_compute").c_str());
		}
	}
	//
	auto is_surface_cell = [&]( const cell_id3 &cell_id ) {
		//
		if( m_param.volume_correct_skip_surfaces ) {
			bool touching_air (false);
			grid.iterate_cell_neighbors(cell_id,[&]( char dim, const cell_id3 &cell_n_id ) {
				if( grid.levelset[cell_n_id.index] > 0.0 ) touching_air = true;
			});
			return touching_air;
		} else {
			return false;
		}
	};
	//
	assert( m_matrix.assembled );
	//
	timer.tick(); console::dump( "Allocating vectors..." );
	auto rhs = m_factory->allocate_vector(grid.valid_cell_count);
	auto result = m_factory->allocate_vector(grid.valid_cell_count);
	console::dump( "Done. Took %s.\n", timer.stock("allocate_vector").c_str());
	//
	timer.tick(); console::dump( "Computing divergence..." );
	//
	grid.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
		if( grid.cell_map[cell_id.index] ) {
			uint_type row = grid.cell_map[cell_id.index]-1;
			if( rhs_correct ) {
				bool solid_cell (false);
				for( char dim : DIMS3 ) {
					grid.iterate_face_neighbors(cell_id,dim,[&]( const face_id3 &face_id ){
						if( grid.area[face_id.index] < 1.0 ) solid_cell = true;
					});
					if( solid_cell ) break;
				}
				if( rhs_correct && ! solid_cell && ! is_surface_cell(cell_id)) {
					const double dx = grid.get_cell_dx(cell_id);
					if( region_count ) {
						rhs->add(row,(dx*dx*dx)*rhs_corrects[regions[cell_id.index]-1]);
					} else {
						rhs->add(row,(dx*dx*dx)*rhs_correct);
					}
				}
			}
			//
			if( m_param.fix_divergence ) {
				double err (0.0);
				if( solid_velocity ) {
					grid.get_unmofidied_divergence(cell_id,[&]( const face_id3 &face_id, double value ) {
						err += value*grid.sample_solid_face_velocity(face_id,solid_velocity);
					});
				}
				bool boundary_flag (false);
				for( int dim : DIMS3 ) {
					boundary_flag = cell_id.pi[dim] == 0 || cell_id.pi[dim] == grid.layers[cell_id.depth]->shape[dim]-1;
					if( boundary_flag ) break;
				}
				if( ! boundary_flag ) rhs->add(row,-err);
			}
			//
			const vec3d p = grid.get_cell_position(cell_id);
			const double dx = grid.get_cell_dx(cell_id);
			//
			grid.get_divergence(cell_id,[&]( const face_id3 &face_id, double value0, double value1 ) {
				bool add_divergence (true);
				if( m_param.special_boundary_condition == 1 && cell_id.pi[0] == grid.layers[cell_id.depth]->shape[0]-1 && face_id.dim == 0 ) {
					add_divergence = false;
				}
				if( m_param.special_boundary_condition == 1 && cell_id.pi[0] == 0 && face_id.dim == 0 ) {
					if( m_moving_solid_func(p) < dx ) {
						add_divergence = false;
					}
				}
				if( add_divergence ) {
					if( solid_velocity && value1 ) {
						if( cell_id.pi[0] != grid.layers[cell_id.depth]->shape[0]-1 ) {
							rhs->add(row,value1*grid.sample_solid_face_velocity(face_id,solid_velocity));
						}
					}
					if( value0 ) {
						rhs->add(row,value0*grid.velocity[face_id.index]);
					}
				}
			});
			//
			if( grid.flux_boundary_condition.has_flux()) {
				const vec3d p = grid.get_cell_position(cell_id);
				for( int dim : DIMS3 ) {
					if( cell_id.pi[dim] == 0 ) {
						rhs->add(row,(dx*dx)*grid.flux_boundary_condition.velocity[dim][0]);
					} else if( cell_id.pi[dim] == grid.layers[cell_id.depth]->shape[dim]-1 ) {
						rhs->add(row,-(dx*dx)*grid.flux_boundary_condition.velocity[dim][1]);
					}
				}
			}
		}
	});
	console::dump( "Done. Took %s.\n", timer.stock("divergence_compute").c_str());
	//
	auto compute_vector_kind = [&]( std::vector<unsigned char> &result ) {
		unsigned total_kinds (0);
		result.resize(grid.valid_cell_count);
		grid.serial_iterate_active_cells([&]( const cell_id3 &cell_id ) {
			if( grid.cell_map[cell_id.index] ) {
				uint_type row = grid.cell_map[cell_id.index]-1;
				result[row] = cell_id.depth;
				total_kinds = std::max(total_kinds,(unsigned)result[row]);
			}
		});
		return total_kinds;
	};
	//
	timer.tick(); console::dump( "Counting and forming the vector norm kind..." );
	std::vector<unsigned char> vector_kind; unsigned total_kinds = compute_vector_kind(vector_kind);
	m_solver->register_vector_norm_kind(vector_kind);
	console::dump( "Done. Total=%u. Took %s.\n", total_kinds, timer.stock("count_vector_norm_kind").c_str());
	//
	timer.tick(); console::dump( "Solving the linear system..." );
	rhs->const_for_each([]( size_t row, double value ){
		assert( ! utility::is_nan(value));
	});
	//
	if( m_param.remove_one_degrees_of_freedom ) {
		m_matrix.Lhs->clear(0);
		m_matrix.Lhs->add_to_element(0,0,1.0);
		rhs->set(0,0);
		m_parallel.for_each(grid.valid_cell_count,[&]( size_t row ) {
			m_matrix.Lhs->for_each(row,[&]( size_t column, double& value ) {
				if( column == 0 ) value = 0.0;
			});
		});
	}
	auto status = m_solver->solve(m_matrix.Lhs.get(),rhs.get(),result.get());
	//
	auto build_residual_string = [&]( const std::vector<double> &rhs ) {
		std::string str;
		if( rhs.size() <= 4 ) {
			for( const auto &e : rhs ) {
				str += " " + console::format_str("%.2e",e);
			}
		} else {
			double max_resid (0.0);
			double min_resid (-1.0);
			unsigned slot_max (0), slot_min(0);
			for( int n=0; n<rhs.size(); ++n ) {
				const auto &e = rhs[n];
				if( e > max_resid ) {
					max_resid = e;
					slot_max = n;
				}
				if( min_resid < 0.0 ) min_resid = max_resid;
				else if( e < min_resid ) {
					min_resid = e;
					slot_min = n;
				}
			}
			str = console::format_str("(num=%u, min[%u]=%.2e, max[%u]=%.2e)",status.vector_reresid.size(),min_resid,slot_min,max_resid,slot_max);
		}
		return str;
	};
	for( int n=0; n<status.vector_reresid.size(); ++n ) {
		console::write(console::format_str("reresid%d",n),status.vector_reresid[n]);
		console::write(console::format_str("absresid%d",n),status.vector_absresid[n]);
	}
	console::dump( "Done. Reresid=%.2e. Iterations = %u. Took %s.\n", status.reresid, status.count, timer.stock("solve_linear_system").c_str());
	console::write("num_iteration_count",status.count);
	//
	if( ! status.vector_reresid.empty()) {
		std::string reresid_str = build_residual_string(status.vector_reresid);
		std::string absresid_str = build_residual_string(status.vector_absresid);
		console::dump( "Reresid=%s\n", reresid_str.c_str());
		console::dump( "Absresid=%s\n", absresid_str.c_str());
	}
	//
	timer.tick(); console::dump( "Updating velocity..." );
	for( uint_type n=0; n<grid.face_count; ++n ) {
		if( ! grid.face_map[n] ) grid.velocity[n] = 0.0;
	}
	grid.iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
		grid.get_scaled_gradient_density(face_id,[&]( const cell_id3 &cell_id, double value, const grid3::gradient_info3_scalar &info ) {
			if( grid.face_map[face_id.index] && grid.cell_map[cell_id.index] && value ) {
				grid.velocity[face_id.index] -= value * result->at(grid.cell_map[cell_id.index]-1);
				console::dump("face_id=%u, cell_id=%u, pressure=%e\n",face_id.index,cell_id.index,result->at(grid.cell_map[cell_id.index]-1));
			}
		});
	});
	//
	if( pressure_vector ) {
		pressure_vector->resize(grid.cell_count);
		for( uint_type n=0; n<grid.cell_count; ++n ) {
			if( grid.cell_map[n] ) {
				(*pressure_vector)[n] = result->at(grid.cell_map[n]-1);
			} else {
				(*pressure_vector)[n] = 0.0;
			}
		}
	}
	//
	console::dump( "Done. Took %s.\n", timer.stock("update_velocity").c_str());
	//
	grid.clear_map();
	console::dump( "<<< Done. Took %s.\n", timer.stock("project").c_str());
}

void macoctreeproject3::project_cloud( grid3 &grid, double dt, int step_count, std::function<vec3d(const vec3d &p)> solid_velocity, std::vector<Real> *pressure ) {
	//
	std::vector<uint_type> regions;
	std::vector<Real> current_volumes;
	std::vector<Real> target_volumes;
	std::vector<Real> y_list;
	//
	project_cloud(grid,dt,step_count, 0,regions,current_volumes,target_volumes,y_list,solid_velocity,pressure);
}
//
void macoctreeproject3::project_cloud( grid3 &grid, double dt, int step_count,
					  size_t region_count,
					  const std::vector<uint_type> &regions,
					  const std::vector<Real> &current_volumes,
					  const std::vector<Real> &target_volumes,
					  std::vector<Real> &y_list,
					  std::function<vec3d(const vec3d &p)> solid_velocity,
					  std::vector<Real> *pressure_vector ) {
	//
	scoped_timer timer(this);
	timer.tick(); console::dump( ">>> Projection step...\n" );
	//
	double rhs_correct (0.0);
	std::vector<Real> rhs_corrects(region_count);
	//
	if( m_param.volume_correction ) {
		const double dt0 = std::max(m_param.minimal_volume_correct_timestep,dt);
		if( region_count ) {
			timer.tick(); console::dump( "Computing regional volume corrector..." );
			for( unsigned n=0; n<region_count; ++n ) {
				double x = (current_volumes[n]-target_volumes[n])/target_volumes[n];
				double y = y_list[n] + x*dt0; y_list[n] = y;
				double kp = -log(1.0-m_param.volume_recover_ratio)/(m_param.np*dt0);
				double ki = kp*kp/16.0;
				rhs_corrects[n] = -(kp*x+ki*y)/(x+1.0);

				assert( ! utility::is_nan(rhs_corrects[n]));
				if( rhs_corrects[n] ) rhs_correct = 1.0;
			}
			console::dump( "Done. Took %s.\n", timer.stock("regional_volume_compute").c_str());
		} else {
			//
			timer.tick(); console::dump( "Computing global volume corrector..." );
			double current_volume = grid.get_volume_density();//ここの値がlevelsetを参照しているためNanになっている
			console::dump( "current volume = %e", current_volume);
			if( ! m_initial_volume ) m_initial_volume = current_volume;
			double x = (current_volume-m_initial_volume)/m_initial_volume;
			double y = m_y_prev + x*dt0; m_y_prev = y;
			double kp = -log(1.0-m_param.volume_recover_ratio)/(m_param.np*dt0);
			double ki = kp*kp/16.0;
			rhs_correct = -(kp*x+ki*y)/(x+1.0);
				// assert( ! utility::is_nan(x));
				// assert( ! utility::is_nan(y));
				// assert( ! utility::is_nan(kp));
				// assert( ! utility::is_nan(ki));
			// console::dump("x=%e, y=%e, kp=%e, ki=%e\n",x,y,kp,ki);
			assert( ! utility::is_nan(rhs_correct));
			//
			const double volume_ratio = current_volume / m_initial_volume;
			console::dump("rhs_correct=%e\n",rhs_correct);
			console::write("current_volume",current_volume);
			console::write("volume_ratio",volume_ratio);
			console::dump( "Done. Took %s. Volume ratio=%.2e\n", volume_ratio, timer.stock("volume_compute").c_str());
		}
	}
	//
	auto is_surface_cell = [&]( const cell_id3 &cell_id ) {
		//
		if( m_param.volume_correct_skip_surfaces ) {
			bool touching_air (false);
			grid.iterate_cell_neighbors(cell_id,[&]( char dim, const cell_id3 &cell_n_id ) {
				if( grid.levelset[cell_n_id.index] > 0.0 ) touching_air = true;
			});
			return touching_air;
		} else {
			return false;
		}
	};
	//
	assert( m_matrix.assembled );
	//
	timer.tick(); console::dump( "Allocating vectors..." );
	auto rhs = m_factory->allocate_vector(grid.valid_cell_count);
	auto result = m_factory->allocate_vector(grid.valid_cell_count);
	console::dump( "Done. Took %s.\n", timer.stock("allocate_vector").c_str());
	//
	console::dump("vector size: %zu\n", grid.valid_cell_count);
	timer.tick(); console::dump( "Computing divergence..." );
	//
	grid.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
		if( grid.cell_map[cell_id.index] ) {
			uint_type row = grid.cell_map[cell_id.index]-1;
			if( rhs_correct ) {
				bool solid_cell (false);
				for( char dim : DIMS3 ) {
					grid.iterate_face_neighbors(cell_id,dim,[&]( const face_id3 &face_id ){
						if( grid.area[face_id.index] < 1.0 ) solid_cell = true;
					});
					if( solid_cell ) break;
				}
				if( rhs_correct && ! solid_cell && ! is_surface_cell(cell_id)) {
					const double dx = grid.get_cell_dx(cell_id);
					if( region_count ) {
						rhs->add(row,(dx*dx*dx)*rhs_corrects[regions[cell_id.index]-1]);
					} else {
						console::dump("rhs_correct=%e\n",rhs_correct);
						rhs->add(row,(dx*dx*dx)*rhs_correct);
					}
				}
			}
			//
			if( m_param.fix_divergence ) {
				double err (0.0);
				if( solid_velocity ) {
					grid.get_unmofidied_divergence(cell_id,[&]( const face_id3 &face_id, double value ) {
						err += value*grid.sample_solid_face_velocity(face_id,solid_velocity);
					});
				}
				bool boundary_flag (false);
				for( int dim : DIMS3 ) {
					boundary_flag = cell_id.pi[dim] == 0 || cell_id.pi[dim] == grid.layers[cell_id.depth]->shape[dim]-1;
					if( boundary_flag ) break;
				}
				if( ! boundary_flag ) rhs->add(row,-err);
			}
			//
			const vec3d p = grid.get_cell_position(cell_id);
			const double dx = grid.get_cell_dx(cell_id);
			//
			grid.get_divergence_qc(cell_id,[&]( const face_id3 &face_id, double value0, double value1 ) {
				bool add_divergence (true);
				if( m_param.special_boundary_condition == 1 && cell_id.pi[0] == grid.layers[cell_id.depth]->shape[0]-1 && face_id.dim == 0 ) {
					add_divergence = false;
				}
				if( m_param.special_boundary_condition == 1 && cell_id.pi[0] == 0 && face_id.dim == 0 ) {
					if( m_moving_solid_func(p) < dx ) {
						add_divergence = false;
					}
				}
				if( add_divergence ) {
					if( solid_velocity && value1 ) {
						if( cell_id.pi[0] != grid.layers[cell_id.depth]->shape[0]-1 ) {
							rhs->add(row,value1*grid.sample_solid_face_velocity(face_id,solid_velocity));
						}
					}
					if( value0 ) {
						//rhs->add(row,value0*grid.velocity[face_id.index]);
						//rhs->add(row,value0*grid.velocity[face_id.index]);
					}
				}
			});
			//
			if( grid.flux_boundary_condition.has_flux()) {
				const vec3d p = grid.get_cell_position(cell_id);
				for( int dim : DIMS3 ) {
					if( cell_id.pi[dim] == 0 ) {
						rhs->add(row,(dx*dx)*grid.flux_boundary_condition.velocity[dim][0]);
					} else if( cell_id.pi[dim] == grid.layers[cell_id.depth]->shape[dim]-1 ) {
						rhs->add(row,-(dx*dx)*grid.flux_boundary_condition.velocity[dim][1]);
					}
				}
			}
		}
	});
	console::dump( "Done. Took %s.\n", timer.stock("divergence_compute").c_str());
	//
	auto compute_vector_kind = [&]( std::vector<unsigned char> &result ) {
		unsigned total_kinds (0);
		result.resize(grid.valid_cell_count);
		grid.serial_iterate_active_cells([&]( const cell_id3 &cell_id ) {
			if( grid.cell_map[cell_id.index] ) {
				uint_type row = grid.cell_map[cell_id.index]-1;
				result[row] = cell_id.depth;
				total_kinds = std::max(total_kinds,(unsigned)result[row]);
			}
		});
		return total_kinds;
	};
	if (step_count == 40) {
			char filename_rhs[256];
			sprintf(filename_rhs, "/home/takebe/shiokaze/project/rhs_step_40.txt");
			std::ofstream file_rhs(filename_rhs);
			file_rhs << "Row RHS_Value\n";
            rhs->const_for_each([&]( size_t row, double value ) {
                file_rhs << row << " " << value << "\n";
            });
			file_rhs.close();
		}
	//
	timer.tick(); console::dump( "Counting and forming the vector norm kind..." );
	std::vector<unsigned char> vector_kind; unsigned total_kinds = compute_vector_kind(vector_kind);
	m_solver->register_vector_norm_kind(vector_kind);
	console::dump( "Done. Total=%u. Took %s.\n", total_kinds, timer.stock("count_vector_norm_kind").c_str());
	//
	timer.tick(); console::dump( "Solving the linear system..." );
	rhs->const_for_each([]( size_t row, double value ){
		assert( ! utility::is_nan(value));
	});
	//
	if( m_param.remove_one_degrees_of_freedom ) {
		m_matrix.Lhs->clear(0);
		m_matrix.Lhs->add_to_element(0,0,1.0);
		rhs->set(0,0);
		m_parallel.for_each(grid.valid_cell_count,[&]( size_t row ) {
			m_matrix.Lhs->for_each(row,[&]( size_t column, double& value ) {
				if( column == 0 ) value = 0.0;
			});
		});
	}

	double sum = 0.0;
	for (size_t n=0; n<grid.valid_cell_count; ++n ) {
		sum += rhs->at(n);
	}
	console::dump("RHS sum=%e\n",sum);
	// auto status = m_solver->solve(m_matrix.Lhs.get(),rhs.get(),result.get());
	// //
	// auto build_residual_string = [&]( const std::vector<double> &rhs ) {
	// 	std::string str;
	// 	if( rhs.size() <= 4 ) {
	// 		for( const auto &e : rhs ) {
	// 			str += " " + console::format_str("%.2e",e);
	// 		}
	// 	} else {
	// 		double max_resid (0.0);
	// 		double min_resid (-1.0);
	// 		unsigned slot_max (0), slot_min(0);
	// 		for( int n=0; n<rhs.size(); ++n ) {
	// 			const auto &e = rhs[n];
	// 			if( e > max_resid ) {
	// 				max_resid = e;
	// 				slot_max = n;
	// 			}
	// 			if( min_resid < 0.0 ) min_resid = max_resid;
	// 			else if( e < min_resid ) {
	// 				min_resid = e;
	// 				slot_min = n;
	// 			}
	// 		}
	// 		str = console::format_str("(num=%u, min[%u]=%.2e, max[%u]=%.2e)",status.vector_reresid.size(),min_resid,slot_min,max_resid,slot_max);
	// 	}
	// 	return str;
	// };
	// for( int n=0; n<status.vector_reresid.size(); ++n ) {
	// 	console::write(console::format_str("reresid%d",n),status.vector_reresid[n]);
	// 	console::write(console::format_str("absresid%d",n),status.vector_absresid[n]);
	// }
	// console::dump( "Done. Reresid=%.2e. Iterations = %u. Took %s.\n", status.reresid, status.count, timer.stock("solve_linear_system").c_str());
	// console::write("num_iteration_count",status.count);
	// //
	// if( ! status.vector_reresid.empty()) {
	// 	std::string reresid_str = build_residual_string(status.vector_reresid);
	// 	std::string absresid_str = build_residual_string(status.vector_absresid);
	// 	console::dump( "Reresid=%s\n", reresid_str.c_str());
	// 	console::dump( "Absresid=%s\n", absresid_str.c_str());
	// }
	//
	timer.tick(); console::dump( "Updating velocity..." );
	for( uint_type n=0; n<grid.face_count; ++n ) {
		if( ! grid.face_map[n] ) grid.velocity[n] = 0.0;
	}
	grid.iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
		grid.get_scaled_gradient_qc(face_id,[&]( const cell_id3 &cell_id, double value, const grid3::gradient_info3_cloud &info ) {
			if( grid.face_map[face_id.index] && grid.cell_map[cell_id.index] && value ) {
				grid.velocity[face_id.index] -= value * result->at(grid.cell_map[cell_id.index]-1);
				// console::dump("face_id=%u, cell_id=%u, pressure=%e\n",face_id.index,cell_id.index,result->at(grid.cell_map[cell_id.index]-1));
			}
		});
	});
	//
	if( pressure_vector ) {
		pressure_vector->resize(grid.cell_count);
		for( uint_type n=0; n<grid.cell_count; ++n ) {
			if( grid.cell_map[n] ) {
				(*pressure_vector)[n] = result->at(grid.cell_map[n]-1);
			} else {
				(*pressure_vector)[n] = 0.0;
			}
		}
	}
	//
	console::dump( "Done. Took %s.\n", timer.stock("update_velocity").c_str());
	//
	grid.clear_map();
	console::dump( "<<< Done. Took %s.\n", timer.stock("project").c_str());
}
void macoctreeproject3::project_merged_cell( grid3 &grid, double dt, bool use_Eigen, int step_count,UnionFind &uf, int all_cell_count,
					  int merged_cell_count,
					  std::vector<cell_id3_and_is_merged> &cell_ids_included_merged_cells,
					  std::vector<int> &merged_cell_id, std::function<vec3d(const vec3d &p)> solid_velocity, std::vector<Real> *pressure) {
	//
	std::vector<uint_type> regions;
	std::vector<Real> current_volumes;
	std::vector<Real> target_volumes;
	std::vector<Real> y_list;
	//
	project_merged_cell(grid,dt,use_Eigen,step_count,0,regions,current_volumes,target_volumes,y_list, uf, all_cell_count, merged_cell_count, cell_ids_included_merged_cells, merged_cell_id, solid_velocity,pressure);
}

template<typename RHS_type>
void macoctreeproject3::assemble_RHS(grid3 &grid, std::vector<uint_type> p, double dx, cell_id3 top1, cell_id3 bottom1, cell_id3 top2, cell_id3 bottom2, char dim, bool upper, RHS_type rhs, UnionFind uf) {
	if (dim == 1 || dim == 2) {
		//p[1]は今回top1, p[0]はbottom1に対応
		//LHSと現状逆になっていることに注意
		cell_id3 p_top = (grid.get_cell_position(top1)[0] < grid.get_cell_position(top2)[0]) ? top1 : top2;
		cell_id3 p_bottom = (grid.get_cell_position(bottom1)[0] > grid.get_cell_position(bottom2)[0]) ? bottom1 : bottom2;
		int n = (p_top.pi - p_bottom.pi)[0];
		double u_avg_upper = (grid.velocity[grid.upper_face_yz_position_of_x(p_top, dim).index] + grid.velocity[grid.upper_face_yz_position_of_x(p_bottom, dim).index]) * 0.5;
		double u_avg_lower = (grid.velocity[grid.lower_face_yz_position_of_x(p_top, dim).index] + grid.velocity[grid.lower_face_yz_position_of_x(p_bottom, dim).index]) * 0.5;
		bool top1_bigger = (grid.get_cell_position(top1)[0] > grid.get_cell_position(top1)[0]);
		bool bottom1_bigger = (grid.get_cell_position(bottom1)[0] > grid.get_cell_position(bottom2)[0]);
		double u_avg;
		int dir;
		double A;
		if (upper) {
			dir = 1;
			u_avg = u_avg_upper;
		}
		else {
			dir = -1;
			u_avg = u_avg_lower;
		}
		if ( !top1_bigger && !bottom1_bigger && !((top1==top2) && (bottom1==bottom2))) {
			double a = (top1.pi[0] - p_bottom.pi[0] + 1) / (double)(top1.pi[0] - bottom1.pi[0] + 1);
			double t = (grid.get_cell_position(top1)[0] - (grid.get_cell_position(p_bottom)[0] + grid.get_cell_position(p_top)[0])/2.0) / (grid.get_cell_position(top1)[0] - grid.get_cell_position(bottom1)[0]);
			A = dx*dx * (p_top.pi[0] - p_bottom.pi[0] + 1);
			if (utility::is_nan(A * dir * a * u_avg * (1-t) / dx) || utility::is_nan( A * dir * a * u_avg * t / dx)) {
				console::dump("1.A=%e, u_avg=%e, a=%e, t=%e, dir=%d\n",A,u_avg,a,t,dir);
			}
			rhs->add(p[1], A * dir * a * u_avg * (1-t) / dx);
			rhs->add(p[0], A * dir * a * u_avg * t / dx);
		}
		else if (top1_bigger && bottom1_bigger) {
			double a = (p_top.pi[0] - bottom1.pi[0] + 1) / (double)(top1.pi[0] - bottom1.pi[0] + 1);
			double t = (grid.get_cell_position(top1)[0] - (grid.get_cell_position(p_bottom)[0] + grid.get_cell_position(p_top)[0])/2.0)/(grid.get_cell_position(top1)[0] - grid.get_cell_position(bottom1)[0]);
			A = dx*dx * (p_top.pi[0] - p_bottom.pi[0] + 1);
			if (utility::is_nan(A * dir * a * u_avg * (1-t) / dx) || utility::is_nan( A * dir * a * u_avg * t / dx)) {
				console::dump("2.A=%e, u_avg=%e, a=%e, t=%e, dir=%d\n",A,u_avg,a,t,dir);
			}
			rhs->add(p[1], A * dir * a * u_avg * (1-t)/ dx);
			rhs->add(p[0], A * dir * a * u_avg * t / dx);
		}
		else {
			A = dx*dx * (p_top.pi[0] - p_bottom.pi[0] + 1);
			if (utility::is_nan(A * dir * u_avg / 2.0 / dx) || utility::is_nan( A * dir * u_avg / 2.0 / dx)) {
				console::dump("3.A=%e, u_avg=%e, dir=%d\n",A,u_avg,dir);
			}
			rhs->add(p[1], A * dir * u_avg / 2.0 / dx);
			rhs->add(p[0], A * dir * u_avg / 2.0 / dx);
		}
	}
	else {
		double u_avg = uf.velocity_sum_x[uf.root(top1.index)] / uf.number_of_cells[uf.root(top1.index)];
		double u_top = grid.velocity[grid.layers[top1.depth]->active_faces[0](top1.pi+vec3i(1,0,0))];
		double u_bottom = grid.velocity[grid.layers[bottom1.depth]->active_faces[0](bottom1.pi+vec3i(0,0,0))];
		if (top1 == bottom1) {
			if (utility::is_nan(dx*dx * u_top / dx) || utility::is_nan(dx*dx * (- u_bottom) / dx)) {
				console::dump("4.u_top=%e, u_bottom=%e\n", u_top, u_bottom);
			}
			rhs->add(p[1], dx*dx * u_top / dx);
			rhs->add(p[0], dx*dx * (- u_bottom) / dx);
		}
		else {
			if (utility::is_nan(dx*dx * (u_top - u_avg) / dx) || utility::is_nan(dx*dx * (u_avg - u_bottom) / dx)) {
				console::dump("5.u_top=%e, u_bottom=%e, u_avg=%e\n", u_top, u_bottom, u_avg);
			}
			rhs->add(p[1], dx*dx * (u_top - u_avg) / dx);
			rhs->add(p[0], dx*dx * (u_avg - u_bottom) / dx);
		}
	}
}
void macoctreeproject3::project_merged_cell( grid3 &grid, double dt, bool use_Eigen, int step_count,
					  size_t region_count,
					  const std::vector<uint_type> &regions,
					  const std::vector<Real> &current_volumes,
					  const std::vector<Real> &target_volumes,
					  std::vector<Real> &y_list,
					  UnionFind &uf,
					  int all_cell_count,
					  int merged_cell_count,
					  std::vector<cell_id3_and_is_merged> &cell_ids_included_merged_cells,
					  std::vector<int> &merged_cell_id,
					  std::function<vec3d(const vec3d &p)> solid_velocity,
					  std::vector<Real> *pressure_vector ) {
	//
	scoped_timer timer(this);
	timer.tick(); console::dump( ">>> Projection step...\n" );
	//
	double rhs_correct (0.0);
	std::vector<Real> rhs_corrects(region_count);
	//
	if( m_param.volume_correction ) {
		const double dt0 = std::max(m_param.minimal_volume_correct_timestep,dt);
		if( region_count ) {
			timer.tick(); console::dump( "Computing regional volume corrector..." );
			for( unsigned n=0; n<region_count; ++n ) {
				double x = (current_volumes[n]-target_volumes[n])/target_volumes[n];
				double y = y_list[n] + x*dt0; y_list[n] = y;
				double kp = -log(1.0-m_param.volume_recover_ratio)/(m_param.np*dt0);
				double ki = kp*kp/16.0;
				rhs_corrects[n] = -(kp*x+ki*y)/(x+1.0);

				assert( ! utility::is_nan(rhs_corrects[n]));
				if( rhs_corrects[n] ) rhs_correct = 1.0;
			}
			console::dump( "Done. Took %s.\n", timer.stock("regional_volume_compute").c_str());
		} else {
			//
			timer.tick(); console::dump( "Computing global volume corrector..." );
			double current_volume = grid.get_volume_density();//ここの値がlevelsetを参照しているためNanになっている
			console::dump( "current volume = %e", current_volume);
			if( ! m_initial_volume ) m_initial_volume = current_volume;
			double x = (current_volume-m_initial_volume)/m_initial_volume;
			double y = m_y_prev + x*dt0; m_y_prev = y;
			double kp = -log(1.0-m_param.volume_recover_ratio)/(m_param.np*dt0);
			double ki = kp*kp/16.0;
			rhs_correct = -(kp*x+ki*y)/(x+1.0);
				// assert( ! utility::is_nan(x));
				// assert( ! utility::is_nan(y));
				// assert( ! utility::is_nan(kp));
				// assert( ! utility::is_nan(ki));
			// console::dump("x=%e, y=%e, kp=%e, ki=%e\n",x,y,kp,ki);
			assert( ! utility::is_nan(rhs_correct));
			//
			const double volume_ratio = current_volume / m_initial_volume;
			console::write("current_volume",current_volume);
			console::write("volume_ratio",volume_ratio);
			console::dump( "Done. Took %s. Volume ratio=%.2e\n", volume_ratio, timer.stock("volume_compute").c_str());
		}
	}
	//
	auto is_surface_cell = [&]( const cell_id3 &cell_id ) {
		//
		if( m_param.volume_correct_skip_surfaces ) {
			bool touching_air (false);
			grid.iterate_cell_neighbors(cell_id,[&]( char dim, const cell_id3 &cell_n_id ) {
				if( grid.levelset[cell_n_id.index] > 0.0 ) touching_air = true;
			});
			return touching_air;
		} else {
			return false;
		}
	};
	//
	assert( m_matrix.assembled );
	//
	//こっからRHSの作成
	timer.tick(); console::dump( "Allocating vectors..." );
	auto rhs = m_factory->allocate_vector(grid.valid_cell_count);
	auto result = m_factory->allocate_vector(grid.valid_cell_count);
	console::dump("vector size: %zu\n", grid.valid_cell_count);
	console::dump( "Done. Took %s.\n", timer.stock("allocate_vector").c_str());
	//
	console::dump("vector size: %zu\n", grid.valid_cell_count);
	timer.tick(); console::dump( "Computing divergence..." );
	//

	//rhsの新バージョン
	grid.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
		uint_type point = grid.cell_map[cell_id.index]-1;
		if (grid.cell_map[cell_id.index]) {
			for (char dim : DIMS3) {
				std::vector<uint_type> p(4); //圧力点のmatrix上のidを格納する配列
				cell_id3 neighbor_cell_id;
				neighbor_cell_id.depth = cell_id.depth;
				neighbor_cell_id.pi = cell_id.pi + vec3i(dim==0, dim == 1,dim == 2);
				if (!grid.layers[neighbor_cell_id.depth]->active_cells.safe_active(neighbor_cell_id.pi)) {
					continue;
				}
				neighbor_cell_id.index = grid.layers[neighbor_cell_id.depth]->active_cells(neighbor_cell_id.pi);
				// p1はbottom1、p2はtop1、p3はbottom2、p4はtop2
				// topはbottom+merged_cell_countで表現
				p[0] = grid.cell_map[uf.bottom[uf.root(cell_id.index)].index] -1;
				p[1] = grid.cell_map[uf.top[uf.root(cell_id.index)].index] -1;
				p[2] = grid.cell_map[uf.bottom[uf.root(neighbor_cell_id.index)].index] -1;
				p[3] = grid.cell_map[uf.top[uf.root(neighbor_cell_id.index)].index] -1;

				bool is_uniform_this = p[0] == p[1];
				bool is_uniform_neighbor = p[2] == p[3];
				// ★ root を事前計算してキャッシュ
				uint_type root_cell = uf.root(cell_id.index);

				uint_type root_neighbor = uf.root(neighbor_cell_id.index);

				// ★ top/bottom を事前取得
				cell_id3 top1 = uf.top[root_cell];
				cell_id3 bottom1 = uf.bottom[root_cell];
				cell_id3 top2 = uf.top[root_neighbor];
				cell_id3 bottom2 = uf.bottom[root_neighbor];
				double dx = grid.get_cell_dx(cell_id);
				//bool upper = true;
				bool is_top1_bigger = (grid.get_cell_position(cell_id)[0] > grid.get_cell_position(top2)[0]);
				bool is_bottom1_bigger = (grid.get_cell_position(bottom1)[0] > grid.get_cell_position(bottom2)[0]);
				//以下はassemble_RHSの中身。関数から呼び出すと遅かったので直接書くことで関数呼び出しオーバーヘッドを削減（？）

				if (dim == 1 || dim == 2) {//x方向以外
					//p[1]は今回top1, p[0]はbottom1に対応
					//LHSと現状逆になっていることに注意

					cell_id3 p_top = (grid.get_cell_position(top1)[0] < grid.get_cell_position(top2)[0]) ? top1 : top2;
					cell_id3 p_bottom = (grid.get_cell_position(bottom1)[0] > grid.get_cell_position(bottom2)[0]) ? bottom1 : bottom2;

					//int n = (p_top.pi - p_bottom.pi)[0];
					if (cell_id.pi[0] == top2.pi[0] || cell_id.pi[0] == bottom2.pi[0]) { //x方向以外の均一セル 対称性◯

						if (grid.face_map[grid.layers[neighbor_cell_id.depth]->active_faces[dim](neighbor_cell_id.pi)])
						rhs->add(point, dx*dx* (grid.velocity[grid.layers[neighbor_cell_id.depth]->active_faces[dim](neighbor_cell_id.pi)] ));

						if (grid.face_map[grid.layers[cell_id.depth]->active_faces[dim](cell_id.pi)])
						rhs->add(point, dx*dx*(- grid.velocity[grid.layers[cell_id.depth]->active_faces[dim](cell_id.pi)]));

					}
					else if ( (is_top1_bigger && is_bottom1_bigger) || (!is_top1_bigger && !is_bottom1_bigger) ) { //mergedセルのうち両方ともtop2/bottom2より大きいor小さい場合
						if (is_top1_bigger) {

							if (grid.face_map[grid.layers[neighbor_cell_id.depth]->active_faces[dim](neighbor_cell_id.pi)])
							rhs->add(point, dx*dx* (grid.velocity[grid.layers[neighbor_cell_id.depth]->active_faces[dim](top2.pi)]));

							if (grid.face_map[grid.layers[cell_id.depth]->active_faces[dim](cell_id.pi)])
							rhs->add(point, dx*dx*( - grid.velocity[grid.layers[cell_id.depth]->active_faces[dim](cell_id.pi)]));

						} else if (!is_top1_bigger) {

							if (grid.face_map[grid.layers[neighbor_cell_id.depth]->active_faces[dim](neighbor_cell_id.pi)])
							rhs->add(point, dx*dx* (grid.velocity[grid.layers[neighbor_cell_id.depth]->active_faces[dim](cell_id.pi)]));

							if (grid.face_map[grid.layers[cell_id.depth]->active_faces[dim](bottom1.pi)])
							rhs->add(point, dx*dx*( - grid.velocity[grid.layers[cell_id.depth]->active_faces[dim](bottom1.pi)]));

						}
					}
				}
				else {//dim==0(x方向)
					double u_upper = grid.velocity[grid.layers[top1.depth]->active_faces[dim](top1.pi+vec3i(1,0,0))];
					double u_lower = grid.velocity[grid.layers[bottom1.depth]->active_faces[dim](bottom1.pi)];
					double u_avg = (u_upper + u_lower) * 0.5;

					if (top1.index==bottom1.index) {//均一セルの場合
						if (grid.face_map[grid.layers[neighbor_cell_id.depth]->active_faces[dim](neighbor_cell_id.pi+vec3i(1,0,0))]) 
						rhs->add(point, dx*dx* ( grid.velocity[grid.layers[neighbor_cell_id.depth]->active_faces[dim](neighbor_cell_id.pi)]));
						if (grid.face_map[grid.layers[cell_id.depth]->active_faces[dim](cell_id.pi)])
						rhs->add(point, dx*dx* (- grid.velocity[grid.layers[cell_id.depth]->active_faces[dim](cell_id.pi)]));
					}
					else {//mergedセルの場合
						if (cell_id.index == top1.index) {
							if (grid.face_map[grid.layers[top1.depth]->active_faces[dim](top1.pi+vec3i(1,0,0))]) 
							rhs->add(point, dx*dx* ( u_upper));
							rhs->add(point, dx*dx* (- u_avg));
						}
						else if (cell_id.index == bottom1.index) {
							if (grid.face_map[grid.layers[bottom1.depth]->active_faces[dim](bottom1.pi)]) 
							rhs->add(point, dx*dx* ( u_avg));
							rhs->add(point, dx*dx* (- u_lower));
						}
					}
				}
			}
		}
	});
		if (step_count == 40) {
			char filename_rhs[256];
			sprintf(filename_rhs, "/home/takebe/shiokaze/project_merged_cell/rhs_step_40.txt");
			std::ofstream file_rhs(filename_rhs);
			file_rhs << "Row RHS_Value\n";
            rhs->const_for_each([&]( size_t row, double value ) {
                file_rhs << row << " " << value << "\n";
            });
			file_rhs.close();
		}
	// rhs設計の旧バージョン
	// grid.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
	// 	for (int dir = -1; dir <= 1; dir += 2) {
	// 		for (char dim : DIMS3) {
	// 			if (should_Lhs_calc(grid, cell_id, uf, dim) && grid.cell_map[uf.top[uf.root(cell_id.index)].index]) {
	// 				std::vector<uint_type> p(4); //圧力点のmatrix上のidを格納する配列
	// 				cell_id3 neighbor_cell_id;
	// 				neighbor_cell_id.depth = cell_id.depth;
	// 				neighbor_cell_id.pi = cell_id.pi + dir * vec3i(dim==0, dim == 1,dim == 2);
	// 				if (!grid.layers[neighbor_cell_id.depth]->active_cells.safe_active(neighbor_cell_id.pi)) {
	// 					continue;
	// 				}
	// 				neighbor_cell_id.index = grid.layers[neighbor_cell_id.depth]->active_cells(neighbor_cell_id.pi);
	// 				// p1はtop1、p2はbottom1、p3はtop2、p4はbottom2
	// 				// topはbottom+merged_cell_countで表現
	// 				p[0] = grid.cell_map[uf.bottom[uf.root(cell_id.index)].index] -1;
	// 				p[1] = grid.cell_map[uf.top[uf.root(cell_id.index)].index] -1;
	// 				p[2] = grid.cell_map[uf.bottom[uf.root(neighbor_cell_id.index)].index] -1;
	// 				p[3] = grid.cell_map[uf.top[uf.root(neighbor_cell_id.index)].index] -1;

	// 				// ★ root を事前計算してキャッシュ
	// 				uint_type root_cell = uf.root(cell_id.index);

	// 				uint_type root_neighbor = uf.root(neighbor_cell_id.index);

	// 				// ★ top/bottom を事前取得
	// 				cell_id3 top1 = uf.top[root_cell];
	// 				cell_id3 bottom1 = uf.bottom[root_cell];
	// 				cell_id3 top2 = uf.top[root_neighbor];
	// 				cell_id3 bottom2 = uf.bottom[root_neighbor];
	// 				double dx = grid.get_cell_dx(cell_id);
	// 				bool upper = (dir == 1);
	// 				//以下はassemble_RHSの中身。関数から呼び出すと遅かったので直接書くことで関数呼び出しオーバーヘッドを削減（？）
	// 				if (dim == 1 || dim == 2) {
	// 					//p[1]は今回top1, p[0]はbottom1に対応
	// 					//LHSと現状逆になっていることに注意
	// 					cell_id3 p_top = (grid.get_cell_position(top1)[0] < grid.get_cell_position(top2)[0]) ? top1 : top2;
	// 					cell_id3 p_bottom = (grid.get_cell_position(bottom1)[0] > grid.get_cell_position(bottom2)[0]) ? bottom1 : bottom2;
	// 					int n = (p_top.pi - p_bottom.pi)[0];
	// 					double u_avg_upper = (grid.velocity[grid.upper_face_yz_position_of_x(p_top, dim).index] + grid.velocity[grid.upper_face_yz_position_of_x(p_bottom, dim).index]) * 0.5;
	// 					double u_avg_lower = (grid.velocity[grid.lower_face_yz_position_of_x(p_top, dim).index] + grid.velocity[grid.lower_face_yz_position_of_x(p_bottom, dim).index]) * 0.5;
	// 					bool top1_bigger = (grid.get_cell_position(top1)[0] > grid.get_cell_position(top1)[0]);
	// 					bool bottom1_bigger = (grid.get_cell_position(bottom1)[0] > grid.get_cell_position(bottom2)[0]);
	// 					double u_avg;
	// 					int dir;
	// 					double A;
	// 					if (upper) {
	// 						dir = 1;
	// 						u_avg = u_avg_upper;
	// 					}
	// 					else {
	// 						dir = -1;
	// 						u_avg = u_avg_lower;
	// 					}
	// 					if ( !top1_bigger && !bottom1_bigger && !((top1==top2) && (bottom1==bottom2))) {
	// 						double a = (top1.pi[0] - p_bottom.pi[0] + 1) / (double)(top1.pi[0] - bottom1.pi[0] + 1);
	// 						double t = (grid.get_cell_position(top1)[0] - (grid.get_cell_position(p_bottom)[0] + grid.get_cell_position(p_top)[0])/2.0) / (grid.get_cell_position(top1)[0] - grid.get_cell_position(bottom1)[0]);
	// 						A = dx*dx * (p_top.pi[0] - p_bottom.pi[0] + 1);
	// 						if (!utility::is_nan(u_avg)) {
	// 							if (use_Eigen) {
	// 								b(p[1]) += A * dir * a * u_avg * (1-t) / dx;
	// 								b(p[0]) += A * dir * a * u_avg * t / dx;
	// 							}
	// 							else {
	// 								rhs->add(p[1], A * dir * a * u_avg * (1-t) / dx);
	// 								rhs->add(p[0], A * dir * a * u_avg * t / dx);
	// 							}

	// 						}
	// 					}
	// 					else if (top1_bigger && bottom1_bigger) {
	// 						double a = (p_top.pi[0] - bottom1.pi[0] + 1) / (double)(top1.pi[0] - bottom1.pi[0] + 1);
	// 						double t = (grid.get_cell_position(top1)[0] - (grid.get_cell_position(p_bottom)[0] + grid.get_cell_position(p_top)[0])/2.0)/(grid.get_cell_position(top1)[0] - grid.get_cell_position(bottom1)[0]);
	// 						A = dx*dx * (p_top.pi[0] - p_bottom.pi[0] + 1);
	// 						if (!utility::is_nan(u_avg)) {
								
	// 							if (use_Eigen) {
	// 								b(p[1]) += A * dir * a * u_avg * (1-t) / dx;
	// 								b(p[0]) += A * dir * a * u_avg * t / dx;
	// 							}
	// 							else {
	// 								rhs->add(p[1], A * dir * a * u_avg * (1-t)/ dx);
	// 								rhs->add(p[0], A * dir * a * u_avg * t / dx);
	// 							}
	// 						}
	// 					}
	// 					else {
	// 						A = dx*dx * (p_top.pi[0] - p_bottom.pi[0] + 1);
	// 						if (!utility::is_nan(u_avg)) {
	// 							if (use_Eigen) {
	// 								b(p[1]) += A * dir * u_avg / 2.0 / dx;
	// 								b(p[0]) += A * dir * u_avg / 2.0 / dx;
	// 							}
	// 							else {
	// 								rhs->add(p[1], A * dir * u_avg / 2.0 / dx);
	// 								rhs->add(p[0], A * dir * u_avg / 2.0 / dx);
	// 							}
	// 						}

	// 					}
	// 				}
	// 				else {
	// 					double u_avg = uf.velocity_sum_x[uf.root(top1.index)] / uf.number_of_cells[uf.root(top1.index)];
	// 					double u_top = grid.velocity[grid.layers[top1.depth]->active_faces[0](top1.pi+vec3i(1,0,0))];
	// 					if (grid.layers[bottom1.depth]->active_faces[dim].safe_active(bottom1.pi)){
							
	// 					}
	// 					double u_bottom = grid.velocity[grid.layers[bottom1.depth]->active_faces[0](bottom1.pi+vec3i(0,0,0))];
						
	// 					if (top1 == bottom1) {
	// 						if (!utility::is_nan(u_top) && !utility::is_nan(u_bottom)) {
	// 							if (use_Eigen) {
	// 								b(p[1]) += dx*dx * u_top / dx;
	// 								b(p[0]) += dx*dx * (- u_bottom) / dx;
	// 							}
	// 							else {
	// 								// rhs->add(p[1], dx*dx * u_top / dx);
	// 								// rhs->add(p[0], dx*dx * (- u_bottom) / dx);
	// 							}
	// 							}
	// 					}
	// 					else {
	// 						if (!utility::is_nan(u_top) && !utility::is_nan(u_avg) && !utility::is_nan(u_bottom)) {
	// 							if (use_Eigen) {
	// 								b(p[1]) += dx*dx * (u_top - u_avg) / dx;
	// 								b(p[0]) += dx*dx * (u_avg - u_bottom) / dx;
	// 							}
	// 							else {
	// 								// rhs->add(p[1], dx*dx * (u_top - u_avg) / dx);
	// 								// rhs->add(p[0], dx*dx * (u_avg - u_bottom) / dx);
	// 							}
	// 						}
	// 					}
	// 				}
	// 							}
	// 						}
	// 					}
	// });
	
	grid.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
		if( grid.cell_map[cell_id.index] ) {
			uint_type row = grid.cell_map[cell_id.index]-1;
			if( rhs_correct ) {
				bool solid_cell (false);
				for( char dim : DIMS3 ) {
					grid.iterate_face_neighbors(cell_id,dim,[&]( const face_id3 &face_id ){
						if( grid.area[face_id.index] < 1.0 ) solid_cell = true;
					});
					if( solid_cell ) break;
				}
				if( rhs_correct && ! solid_cell && ! is_surface_cell(cell_id)) {
					const double dx = grid.get_cell_dx(cell_id);
					if( region_count ) {
						rhs->add(row,(dx*dx*dx)*rhs_corrects[regions[cell_id.index]-1]);
					} else {
						console::dump("rhs_correct=%e\n",rhs_correct);
						if (use_Eigen) {
							b(row) += (dx*dx*dx)*rhs_correct;
						}
						else {
							rhs->add(row,(dx*dx*dx)*rhs_correct);
						}
					}
				}
			}
			
			if( m_param.fix_divergence ) {
				double err (0.0);
				if( solid_velocity ) {
					grid.get_unmofidied_divergence(cell_id,[&]( const face_id3 &face_id, double value ) {
						err += value*grid.sample_solid_face_velocity(face_id,solid_velocity);
					});
				}
				bool boundary_flag (false);
				for( int dim : DIMS3 ) {
					boundary_flag = cell_id.pi[dim] == 0 || cell_id.pi[dim] == grid.layers[cell_id.depth]->shape[dim]-1;
					if( boundary_flag ) break;
				}
				if( ! boundary_flag ) {
					if (use_Eigen) {
						b(row) += -err;
					}
					else {
						rhs->add(row,-err);
					}
				}
			}
			//
			const vec3d p = grid.get_cell_position(cell_id);
			const double dx = grid.get_cell_dx(cell_id);
			//
			// grid.get_divergence_qc(cell_id,[&]( const face_id3 &face_id, double value0, double value1 ) {
			// 	bool add_divergence (true);
			// 	if( add_divergence ) {
			// 		if( solid_velocity && value1 ) {
			// 			if( cell_id.pi[0] != grid.layers[cell_id.depth]->shape[0]-1 ) {
			// 				if (use_Eigen) {
			// 					b(row) += value1*grid.sample_solid_face_velocity(face_id,solid_velocity);
			// 				}
			// 				else {
			// 					rhs->add(row,value1*grid.sample_solid_face_velocity(face_id,solid_velocity));
			// 				}
			// 			}
			// 		}
			// 		if( value0 ) {
			// 			if (use_Eigen) {
			// 				b(row) += value0*grid.velocity[face_id.index];
			// 			}
			// 			else {
			// 				rhs->add(row,value0*grid.velocity[face_id.index]);
			// 			}
			// 		}
			// 	}
			// });
		}
	});
	console::dump( "Done. Took %s.\n", timer.stock("divergence_compute").c_str());
	//
	auto compute_vector_kind = [&]( std::vector<unsigned char> &result ) {
		unsigned total_kinds (0);
		result.resize(grid.valid_cell_count);
		grid.serial_iterate_active_cells([&]( const cell_id3 &cell_id ) {
			if( grid.cell_map[cell_id.index] ) {
				uint_type row = grid.cell_map[cell_id.index]-1;
				result[row] = cell_id.depth;
				total_kinds = std::max(total_kinds,(unsigned)result[row]);
			}
		});
		return total_kinds;
	};
	//
	timer.tick(); console::dump( "Counting and forming the vector norm kind..." );
	std::vector<unsigned char> vector_kind; unsigned total_kinds = compute_vector_kind(vector_kind);
	m_solver->register_vector_norm_kind(vector_kind);
	console::dump( "Done. Total=%u. Took %s.\n", total_kinds, timer.stock("count_vector_norm_kind").c_str());
	//
	timer.tick(); console::dump( "Solving the linear system..." );
	rhs->const_for_each([]( size_t row, double value ){
		assert( ! utility::is_nan(value));
	});
	//
	auto rhs_dummy = m_factory->allocate_vector(grid.valid_cell_count);
	double sum = 0.0;
	for (size_t n=0; n<grid.valid_cell_count; ++n ) {
		sum += rhs->at(n);
	}
	console::dump("RHS sum=%e\n",sum);
	auto status = m_solver->solve(m_matrix.Lhs.get(),rhs.get(),result.get());
	//bool is_solved = solveWithICCG(A,b,x);
	//
	// if (!is_solved) {
	// 	console::dump("Linear system not solved!\n");
	// }

	// auto build_residual_string = [&]( const std::vector<double> &rhs ) {
	// 	std::string str;
	// 	if( rhs.size() <= 4 ) {
	// 		for( const auto &e : rhs ) {
	// 			str += " " + console::format_str("%.2e",e);
	// 		}
	// 	} else {
	// 		double max_resid (0.0);
	// 		double min_resid (-1.0);
	// 		unsigned slot_max (0), slot_min(0);
	// 		for( int n=0; n<rhs.size(); ++n ) {
	// 			const auto &e = rhs[n];
	// 			if( e > max_resid ) {
	// 				max_resid = e;
	// 				slot_max = n;
	// 			}
	// 			if( min_resid < 0.0 ) min_resid = max_resid;
	// 			else if( e < min_resid ) {
	// 				min_resid = e;
	// 				slot_min = n;
	// 			}
	// 		}
	// 		str = console::format_str("(num=%u, min[%u]=%.2e, max[%u]=%.2e)",status.vector_reresid.size(),min_resid,slot_min,max_resid,slot_max);
	// 	}
	// 	return str;
	// };
	// for( int n=0; n<status.vector_reresid.size(); ++n ) {
	// 	console::write(console::format_str("reresid%d",n),status.vector_reresid[n]);
	// 	console::write(console::format_str("absresid%d",n),status.vector_absresid[n]);
	// }
	console::dump( "Done. Reresid=%.2e. Iterations = %u. Took %s.\n", status.reresid, status.count, timer.stock("solve_linear_system").c_str());
	console::write("num_iteration_count",status.count);
	//
	// if( ! status.vector_reresid.empty()) {
	// 	std::string reresid_str = build_residual_string(status.vector_reresid);
	// 	std::string absresid_str = build_residual_string(status.vector_absresid);
	// 	console::dump( "Reresid=%s\n", reresid_str.c_str());
	// 	console::dump( "Absresid=%s\n", absresid_str.c_str());
	// }
	
	timer.tick(); console::dump( "Updating velocity..." );
	for( uint_type n=0; n<grid.face_count; ++n ) {
		if( ! grid.face_map[n] ) grid.velocity[n] = 0.0;
	}
	// grid.iterate_active_faces([&]( const face_id3 &face_id, int tid ) {
	// 	grid.get_scaled_gradient_qc(face_id,[&]( const cell_id3 &cell_id, double value, const grid3::gradient_info3_cloud &info ) {
	// 		if( grid.face_map[face_id.index] && grid.cell_map[cell_id.index] && value ) {
	// 			grid.velocity[face_id.index] -= value * result->at(grid.cell_map[cell_id.index]-1);
	// 		}
	// 	});
	// });
	// //
	// if( pressure_vector ) {
	// 	pressure_vector->resize(grid.cell_count);
	// 	for( uint_type n=0; n<grid.cell_count; ++n ) {
	// 		if( grid.cell_map[n] ) {
	// 			(*pressure_vector)[n] = result->at(grid.cell_map[n]-1);
	// 		} else {
	// 			(*pressure_vector)[n] = 0.0;
	// 		}
	// 	}
	// }
	//
	console::dump( "Done. Took %s.\n", timer.stock("update_velocity").c_str());
	//
	grid.clear_map();
	console::dump( "<<< Done. Took %s.\n", timer.stock("project").c_str());
}


//