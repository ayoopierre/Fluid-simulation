#ifndef MY_FLUID_SIMULATION_BACKEND_CPU
#define MY_FLUID_SIMULATION_BACKEND_CPU

#include <memory>

#include "fluid_simulation_backend.hpp"
#include "csr_matrix_cpu.hpp"

#define FLUID_SIM_LHS1_NNZ 11
#define FLUID_SIM_LHS2_NNZ 11
#define FLUID_SIM_LHS3_NNZ 9
#define FLUID_SIM_MATRIX_NNZ 31
#define FLUID_SIM_EQ_TYPES 3

#define U_OFFSET 0
#define V_OFFSET (width * height)
#define RHO_OFFSET (2 * width * height)

#define AT(i, j) (j * width) + i

class FluidSimulationBackendCPU : FluidSimulationBackend{
    public:
        FluidSimulationBackendCPU(int x_resolution, int y_resolution);
        
        void step(double dt);
    protected:
        void apply_user_input();
        void build_CSR_matrix();
        void run_BiCSTAB();
        void update_pressure();
        void write_heatmap();

        std::shared_ptr<CSRMatrixCPU> A;
        std::vector<double> RHS;
};

#endif