#include "BiCGSTAB_solver_CPU.hpp"
#include "utils.hpp"

#include <math.h>
#include <cstdio>

std::optional<BiCGSTABSolverCpu> BiCGSTABSolverCpu::create(int n)
{
    BiCGSTABSolverCpu instance;

    try{
        instance.x.resize(n);
        instance.r.resize(n);
        instance.r_hat.resize(n);
        instance.p.resize(n);
        instance.v.resize(n);
        instance.h.resize(n);
        instance.s.resize(n);
        instance.t.resize(n);
    }
    catch(std::exception& e){
        return std::nullopt;
    }

    return instance;
}

std::optional<std::shared_ptr<BiCGSTABSolverCpu>> BiCGSTABSolverCpu::create_shared(int n)
{
    std::shared_ptr<BiCGSTABSolverCpu> instance = std::shared_ptr<BiCGSTABSolverCpu>(new BiCGSTABSolverCpu());
    try{
        instance->x.resize(n);
        instance->r.resize(n);
        instance->r_hat.resize(n);
        instance->p.resize(n);
        instance->v.resize(n);
        instance->h.resize(n);
        instance->s.resize(n);
        instance->t.resize(n);
    }
    catch(std::exception& e){
        return std::nullopt;
    }

    return instance;
}

double BiCGSTABSolverCpu::solve(CSRMatrixCPU &A, std::vector<double> &b, int iter, double prec)
{
    double prec_sq;
    /* Start with initial guess (0, 0 ..., 0) */
    std::fill(x.begin(), x.end(), 0.0f);
    /* Set initial residum vector */
    copy_a_to_b(b, r);
    /* Find r_hat such that <r_hat, t> != 0 */
    set_one_on_first_non_zero(r, r_hat);
    /* Set rho0 */
    rho = vec_dot_vec(r, r_hat);
    /* Set p = r */
    std::memcpy(p.data(), r.data(), sizeof(double) * p.size());
    /* Main solution loop */
    for(int i = 0; i < iter; i++){
        vec_CSRmat_prod(A, p, v);
        alpha = rho / vec_dot_vec(r_hat, v);
        vec_scale(alpha, p);
        /* At this point p is scaled by alpha */
        vec_plus_vec(x, p, h);
        vec_scale(-alpha, v);
        /* At this point v is scaled by -alpha */
        vec_plus_vec(r, v, s);
        /* Check if s meets expected precisiton */
        prec_sq = vec_dot_vec(s, s);
        if(prec_sq < prec * prec){
            std::memcpy(x.data(), h.data(), sizeof(double) * x.size());
            return std::sqrt(prec_sq);
        }
        vec_CSRmat_prod(A, s, t);
        omega = vec_dot_vec(t, s) / vec_dot_vec(t, t);
        vec_scale(omega, s);
        /* At this point s is scaled by omega*/
        vec_plus_vec(h, s, x);
        vec_scale(-omega, t);
        vec_plus_vec(s, t, r);
        /* Check if residum small enough */
        prec_sq = vec_dot_vec(r, r);
        if(prec_sq < prec * prec){
            /* x has already solution of requested precision */
            return std::sqrt(prec_sq);
        }
        prev_rho = rho;
        rho = vec_dot_vec(r, r_hat);
        /* We can skip alpha since we already scaled p */
        beta = (rho / prev_rho) * omega;
        /* We scaled p by alpha previously */
        vec_scale(beta, p);
        vec_scale(omega / alpha, v);
        vec_plus_vec(p, v, p);
        vec_plus_vec(p, r, p);
    }

    return std::sqrt(prec_sq);
}   

