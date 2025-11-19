#include <vector>
#include <optional>
#include <cstdio>

#include "csr_matrix_cpu.hpp"
#include "BiCGSTAB_solver_CPU.hpp"

int main(void)
{
    std::vector<double> val = {12.0, 1.5, 3.0, 2.5, 10.0, 6.0};
    std::vector<int> col = {0, 3, 1, 0, 2, 3};
    std::vector<int> row_idx = {0, 3, 1, 0, 2, 3};

    std::optional<CSRMatrixCPU> opt_A = CSRMatrixCPU::create(0, 0);
    if (!opt_A.has_value())
    {
        std::printf("Failed to create CSR matrix\n");
        return -1;
    }

    CSRMatrixCPU A = (opt_A.value());
    // CSRMatrixCPU A = std::move(opt_A.value());
    A.rows = 4;
    A.nnz = col.size();
    A.col = std::move(col);
    A.val = std::move(val);
    A.row_ptr = std::move(row_idx);
    // std::swap(A.val, val);

    std::vector<double> b = {1.0, 2.0, 3.0, 4.0};

    std::optional<BiCGSTABSolverCpu> opt_solver = BiCGSTABSolverCpu::create(4);

    if (!opt_solver.has_value())
    {
        std::printf("Failed to create solver\n");
        return -1;
    }

    BiCGSTABSolverCpu solver = std::move(opt_solver.value());

    double prec = solver.solve(A, b, 10, 0.0000005);

    for (int i = 0; i < 4; i++)
    {
        std::printf("%lf ", solver.x[i]);
    }
    std::printf("\n");
    std::printf("Prec: %lf\n", prec);

    return 0;
}