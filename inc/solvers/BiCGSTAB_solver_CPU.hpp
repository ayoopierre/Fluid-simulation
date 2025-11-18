#ifndef MY_BICGSTAB_SOLVER_CPU
#define MY_BICGSTAB_SOLVER_CPU

#include <memory>
#include <vector>
#include <optional>

#include "csr_matrix_cpu.hpp"

class BiCGSTABSolverCpu{
    public:
        /*
        Create instance of solver in CPU context.
        For performance all used objects will be pre-allocated for given
        size of system to be solved. Since allocation of underlying resources
        might fail if system is very large, we will controll if the object
        is actually instantiated.
        */
        static std::optional<BiCGSTABSolverCpu> create(int n);
        static std::optional<std::shared_ptr<BiCGSTABSolverCpu>> create_shared(int n);

        /* Solve linear system using, BiCGSTAB method, returns norm2 of residum vector */
        double solve(CSRMatrixCPU& A, std::vector<double>& b, int iter, double prec);

        inline BiCGSTABSolverCpu(BiCGSTABSolverCpu& other){
            x = other.x;
            r = other.r;
            r_hat = other.r_hat;
            p = other.p;
            v = other.v;
            h = other.h;
            s = other.s;
            t = other.t;
        }

        inline BiCGSTABSolverCpu& operator=(BiCGSTABSolverCpu& other){
            x = other.x;
            r = other.r;
            r_hat = other.r_hat;
            p = other.p;
            v = other.v;
            h = other.h;
            s = other.s;
            t = other.t;

            return *this;
        }

        inline BiCGSTABSolverCpu(BiCGSTABSolverCpu&& other){
            x = std::move(other.x);
            r = std::move(other.r);
            r_hat = std::move(other.r_hat);
            p = std::move(other.p);
            v = std::move(other.v);
            h = std::move(other.h);
            s = std::move(other.s);
            t = std::move(other.t);
        }

        inline BiCGSTABSolverCpu& operator=(BiCGSTABSolverCpu&& other){
            x = std::move(other.x);
            r = std::move(other.r);
            r_hat = std::move(other.r_hat);
            p = std::move(other.p);
            v = std::move(other.v);
            h = std::move(other.h);
            s = std::move(other.s);
            t = std::move(other.t);

            return *this;
        }

        std::vector<double> x;
        std::vector<double> r;
        std::vector<double> r_hat;
        std::vector<double> p;
        std::vector<double> v;
        std::vector<double> h;
        std::vector<double> s;
        std::vector<double> t;

        double alpha;
        double beta;
        double prev_rho;
        double rho;
        double omega;
    private:
        BiCGSTABSolverCpu(){};
};

#endif