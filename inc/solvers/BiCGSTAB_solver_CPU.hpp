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

        /* We need custom assigment to avoid copies */
        BiCGSTABSolverCpu(BiCGSTABSolverCpu& other);
        BiCGSTABSolverCpu& operator=(BiCGSTABSolverCpu& other);

        BiCGSTABSolverCpu(BiCGSTABSolverCpu&& other);
        BiCGSTABSolverCpu& operator=(BiCGSTABSolverCpu&& other);

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