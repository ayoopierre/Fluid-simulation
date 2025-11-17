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
        static std::shared_ptr<BiCGSTABSolverCpu> create_shared(int n);

        std::vector<double> solve(CSRMatrixCPU& A, std::vector<double>& b);

        std::vector<double> x;
};

#endif