#ifndef MY_FLUID_SIMULATION_BACKEND_CPU
#define MY_FLUID_SIMULATION_BACKEND_CPU

#include <memory>
#include <array>
#include <cstdio>
#include <cmath>
#include <sstream>

#include "fluid_simulation_backend.hpp"
#include "csr_matrix_cpu.hpp"
#include "BiCGSTAB_solver_CPU.hpp"
#include "stencil_cpu.hpp"
#include "utils.hpp"

#define FLUID_SIM_LHS1_NNZ 11
#define FLUID_SIM_LHS2_NNZ 11
#define FLUID_SIM_LHS3_NNZ 11
#define FLUID_SIM_MATRIX_NNZ 33
#define FLUID_SIM_EQ_TYPES 3

#define U_OFFSET 0
#define V_OFFSET (width * height)
#define RHO_OFFSET (2 * width * height)

#define AT(i, j) (((j) * width) + (i))

#define SKIP_IF_NOT_IN_BOUNDS(i, j)                                      \
    {                                                                    \
        if ((i) < 0 || (i) > (width - 1) || (j) < 0 || (j) > height - 1) \
            continue;                                                    \
    }

class FluidSimulationBackendCPU : public FluidSimulationBackend
{
public:
    FluidSimulationBackendCPU(int x_resolution, int y_resolution);

    void step()
    {


        for (int j = 0; j < height; j++)
        {
            for (int i = 0; i < width; i++)
            {
                double u = std::rand() / (double)RAND_MAX;
                u -= 0.5;

                solver->x[AT(i, j)] = (i != 0 && i != width - 1 && j != 0 && j != height - 1) ? u : 0.0;

                double v = std::rand() / (double)RAND_MAX;
                v -= 0.5;

                solver->x[V_OFFSET + AT(i, j)] = (i != 0 && i != width - 1 && j != 0 && j != height - 1) ? v : 0.0;

                double density = std::rand() / (double)RAND_MAX;
                density += 0.2;

                solver->x[RHO_OFFSET + AT(i, j)] = (i != 0 && i != width - 1 && j != 0 && j != height - 1) ? density : 0.0;
            }
        }

        // writeStateToCSV("pre.csv", solver->x, width, height, 0);
        // writeStateToCSV("rho_pre.csv", solver->x, width, height, RHO_OFFSET);

        build_CSR_matrix();
        writeCSRtoCSV("matrix.csv", A->val, A->col, A->row_ptr, FLUID_SIM_EQ_TYPES * width * height,
            FLUID_SIM_EQ_TYPES * width * height);

        // exit(1);

        // writeCSRtoCSV("temp.csv", A->val, A->col, A->row_ptr, 3 * width * height, 3 * width * height);
        double prec = solver->solve(*A.get(), RHS, 30, 1e-9);

        // writeStateToCSV("post.csv", solver->x, width, height, 0);
        // writeStateToCSV("rho_post.csv", solver->x, width, height, RHO_OFFSET);
#ifdef VERBOSE
        std::printf("Managed to get %lf precision\n", prec);
#endif

        for (int i = 0; i < 1000; i++)
        {
            build_CSR_matrix();
            prec = solver->solve(*A.get(), RHS, 30, 1e-9);

            std::ostringstream name_matrix;
            std::ostringstream name_u;
            std::ostringstream name_v;
            std::ostringstream name_rho;
            name_v << "v" << i << ".csv";
            name_u << "u" << i << ".csv"; 
            name_rho << "rho" << i << ".csv";
            name_matrix << "matrix" << i << ".csv";

            writeCSRtoCSV(name_matrix.str(), A->val, A->col, A->row_ptr, FLUID_SIM_EQ_TYPES * width * height, FLUID_SIM_EQ_TYPES * width * height);
            writeStateToCSV(name_v.str(), solver->x, width, height, V_OFFSET);
            writeStateToCSV(name_u.str(), solver->x, width, height, 0);
            writeStateToCSV(name_rho.str(), solver->x, width, height, RHO_OFFSET);
#ifdef VERBOSE
            std::printf("%d Managed to get %5.14lf precision\n", i, prec);
#endif
        }
    };
    void temp() { build_CSR_matrix(); };

protected:
    void init_stencils();
    void init_stencils_2();
    void init_test_stencils();

    void apply_user_input();
    void build_CSR_matrix();
    void apply_wall_conditions() {};
    void run_BiCSTAB();
    void update_pressure() {};
    void write_heatmap() {};

    std::shared_ptr<BiCGSTABSolverCpu> solver;

    std::shared_ptr<CSRMatrixCPU> A;
    std::vector<double> RHS;

    std::vector<bool> is_wall;

    double dt = 1e-4;
    // double D = 1.27e-3; /* Diffusion rate */
    double mu = 1e-3;    /* Viscosities */
    double lambda = -2.0 / 3.0 * mu;
    #define DX 1e-3
    #define DY 1e-3

    double dx = DX;
    double dy = DY;
    double dx_sq = DX * DX;
    double dy_sq = DY * DY;
    double cs_sq = 1481.0; /* Squared speed of sound */

    /* Base discretized Navier-Stokes stencil */
    std::vector<StencilCpu> u_type_stencil;
    std::vector<StencilCpu> v_type_stencil;
    std::vector<StencilCpu> rho_type_stencil;

    /* Wall stencils */
    std::vector<StencilCpu> wall_u_type_stencil;
    std::vector<StencilCpu> wall_v_type_stencil;
    std::vector<StencilCpu> wall_rho_type_stencil;
};

#endif