#include "fluid_simulation_backend_cpu.hpp"

#include <optional>
#include <cassert>

FluidSimulationBackendCPU::FluidSimulationBackendCPU(int x_resolution, int y_resolution)
{
    width = x_resolution;
    height = y_resolution;

    u.reserve(width * height);
    v.reserve(width * height);
    rho.reserve(width * height);

    user_u.reserve(width * height);
    user_v.reserve(width * height);
    user_rho.reserve(width * height);

    /*
    This is method specific. From linear system of equations
    that will be solved in this simulation we preciesly know
    how many non-zero entries there is in the matrix. For
    LHS1 and LHS2 equations there will be 11 non-zero entries
    per row, and for LHS3 there will be 9 non-zero entries per
    row and there is width x height equations of each type.
    Therfore we have 31 * width * height non-zero entries in
    whole matrix. 
    */
    std::optional a = CSRMatrixCPU::create_shared(
        FLUID_SIM_EQ_TYPES * width * height,
        FLUID_SIM_MATRIX_NNZ * width * height
    );
    assert(a != std::nullopt);
    A = a.value();

    RHS.resize(FLUID_SIM_EQ_TYPES * width * height);
}

void FluidSimulationBackendCPU::apply_user_input()
{
    for(int i = 0; i < width * height; i++){
        u[i] += user_u[i];
        v[i] += user_v[i];
        rho[i] += user_rho[i];
    }
}

void FluidSimulationBackendCPU::build_CSR_matrix()
{
    int offset = 0;
    /* LHS1 equations */
    for(int j = 0; j < height; j++){
        for(int i = 0; i < width; i++){
            A->row_ptr[offset + i * j] = FLUID_SIM_LHS1_NNZ * i * j;
            /* Cols of U components */
            A->col[offset + i * j + 0];
            A->val[offset + i * j + 0];

            A->col[offset + i * j + 1];
            A->val[offset + i * j + 1];

            A->col[offset + i * j + 2];
            A->val[offset + i * j + 2];

            A->col[offset + i * j + 3];
            A->val[offset + i * j + 3];

            A->col[offset + i * j + 4];
            A->val[offset + i * j + 4];
            /* Cols of V components */

            /* Cols of RHO components */
        }
    }

    offset += width * height;
    /* LHS2 equations */
    for(int j = 0; j < height; j++){
        for(int i = 0; i < width; i++){
            A->row_ptr[offset + i * j] = FLUID_SIM_LHS2_NNZ * i * j;
            /* Cols of U components */

            /* Cols of V components */

            /* Cols of RHO components */
        }
    }

    offset += width * height;
    /* LHS3 equations */
    for(int j = 0; j < height; j++){
        for(int i = 0; i < width; i++){
            A->row_ptr[offset + i * j] = FLUID_SIM_LHS3_NNZ * i * j;
            /* Cols of U components */

            /* Cols of V components */

            /* Cols of RHO components */
        }
    }
}
