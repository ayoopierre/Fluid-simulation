#ifndef MY_FLUID_SIMULATION_BACKEND_CPU
#define MY_FLUID_SIMULATION_BACKEND_CPU

#include <memory>
#include <array>

#include "fluid_simulation_backend.hpp"
#include "csr_matrix_cpu.hpp"
#include "BiCGSTAB_solver_CPU.hpp"

#define FLUID_SIM_LHS1_NNZ 11
#define FLUID_SIM_LHS2_NNZ 11
#define FLUID_SIM_LHS3_NNZ 9
#define FLUID_SIM_MATRIX_NNZ 31
#define FLUID_SIM_EQ_TYPES 3

#define U_OFFSET 0
#define V_OFFSET (width * height)
#define RHO_OFFSET (2 * width * height)

#define AT(i, j) (j * width) + i

class FluidSimulationBackendCPU : FluidSimulationBackend
{
public:
    FluidSimulationBackendCPU(int x_resolution, int y_resolution);

    void step(double dt);

protected:
    void init_stencils();

    void apply_user_input();
    void build_CSR_matrix();
    void build_CSR_matrix_2();
    void apply_wall_conditions();
    void run_BiCSTAB();
    void update_pressure();
    void write_heatmap();

    std::shared_ptr<BiCGSTABSolverCpu> solver;

    std::shared_ptr<CSRMatrixCPU> A;
    std::vector<double> RHS;
    
    std::vector<bool> is_wall;

    double dt;
    double D;  /* Diffusion rate */
    double mu; /* Viscosities */
    double lambda;
    double dx;
    double dy;
    double dx_sq;
    double dy_sq;
    double cs_sq; /* Squared speed of sound */

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