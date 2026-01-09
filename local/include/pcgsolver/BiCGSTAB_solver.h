#ifndef BICGSTAB_SOLVER_H
#define BICGSTAB_SOLVER_H

#include <cmath>
#include <functional>
#include <numeric>
#include <limits>
#include "sparse_matrix.h"
#include "blas_wrapper.h"
#include "pcg_solver.h"

// BiCGSTAB (BiCG Stabilized) solver with Modified Incomplete Cholesky preconditioner
template <class T>
struct BiCGSTABSolver
{
    BiCGSTABSolver(void)
    {
        set_solver_parameters(1e-5, 100, 0.97, 0.25);
    }

    void set_solver_parameters(T tolerance_factor_, int max_iterations_, T modified_incomplete_cholesky_parameter_=0.97, T min_diagonal_ratio_=0.25)
    {
        tolerance_factor=tolerance_factor_;
        if(tolerance_factor<1e-30) tolerance_factor=1e-30;
        max_iterations=max_iterations_;
        modified_incomplete_cholesky_parameter=modified_incomplete_cholesky_parameter_;
        min_diagonal_ratio=min_diagonal_ratio_;
    }

    bool solve(const SparseMatrix<T> &matrix, const std::vector<T> &rhs, std::vector<T> &result, T &residual_out, int &iterations_out, double abs_tol=0.0 ) 
    {
        printf( "Starting BiCGSTAB solver...\n" );
        unsigned int n=matrix.n;
        if(r.size()!=n){ r.resize(n); r0.resize(n); p.resize(n); v.resize(n); s.resize(n); t.resize(n); y.resize(n); z.resize(n); }
        
        zero(result);
        r=rhs;
        r0=r;  // r0* = initial residual (arbitrary)
        residual_out=BLAS::abs_max(r);
        T residual0 = residual_out;
        
        if(residual_out==0) {
            iterations_out=0;
            return true;
        }
        
        double tol= abs_tol ? abs_tol : tolerance_factor*residual_out;
        
        form_preconditioner(matrix);
        apply_preconditioner(r, y);
        p=y;
        
        fixed_matrix.construct_from_matrix(matrix);
        T rho=BLAS::dot(r0, r);
        T alpha, omega;
        int iteration;
        
        for(iteration=0; iteration<max_iterations; ++iteration){
            multiply(fixed_matrix, p, v);
            apply_preconditioner(v, z);
            
            T pv_dot=BLAS::dot(r0, z);
            if(std::abs(pv_dot)<1e-30) {
                iterations_out=iteration;
                residual_out /= residual0;
                printf( "BiCGSTAB did not converge because of small pv_dot. Iteration %d: residual = %e\n", iteration, residual_out );
                return false;
            }
            
            alpha=rho/pv_dot;
            BLAS::add_scaled(-alpha, z, r);  // s=r-alpha*v
            s=r;
            
            residual_out=BLAS::abs_max(s);
            if(residual_out<=tol) {
                BLAS::add_scaled(alpha, p, result);
                iterations_out=iteration+1;
                residual_out /= residual0;
                printf( " converged. Iteration %d: residual = %e alpha = %e rho = %e\n", iteration, residual_out, alpha, rho );
                return true;
            }
            
            multiply(fixed_matrix, s, v);
            //apply_preconditioner(v, z);
            
            T st_dot=BLAS::dot(z, s);
            T tt_dot=BLAS::dot(z, z);
            
            if(std::abs(tt_dot)<1e-30) {
                iterations_out=iteration;
                residual_out /= residual0;
                return false;
            }
            
            omega=st_dot/tt_dot;
            
            BLAS::add_scaled(alpha, p, result);
            BLAS::add_scaled(omega, s, result);
            BLAS::add_scaled(-omega, z, r);
            
            residual_out=BLAS::abs_max(r);
            if(residual_out<=tol) {
                iterations_out=iteration+1;
                residual_out /= residual0;
                return true;
            }
            
            if(iteration % 5000 == 0)
                printf("  BiCGSTAB Iteration %d: residual = %e alpha = %e omega = %e\n", iteration, residual_out, alpha, omega);
            
            T rho_new=BLAS::dot(r0, r);
            T beta=(rho_new/rho)*(alpha/omega);
            rho=rho_new;
            
            // p = y + beta*(p - omega*v)
            BLAS::add_scaled(-omega, z, p);
            BLAS::add_scaled(beta, p, y);
            p=y;
        }
        
        iterations_out=iteration;
        residual_out /= residual0;
        printf( "BiCGSTAB did not converge. Iteration %d: residual = %e\n", iteration, residual_out );
        return false;
    }

    protected:

    SparseColumnLowerFactor<T> ic_factor;
    std::vector<T> r, r0, p, v, s, t, y, z;
    FixedSparseMatrix<T> fixed_matrix;

    T tolerance_factor;
    int max_iterations;
    T modified_incomplete_cholesky_parameter;
    T min_diagonal_ratio;

    void form_preconditioner(const SparseMatrix<T>& matrix)
    {
        factor_modified_incomplete_cholesky0(matrix, ic_factor);
    }

    void apply_preconditioner(const std::vector<T> &x, std::vector<T> &result)
    {
        solve_lower(ic_factor, x, result);
        solve_lower_transpose_in_place(ic_factor, result);
    }
};

#endif