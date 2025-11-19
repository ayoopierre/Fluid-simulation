#include "fluid_simulation_backend_cpu.hpp"

#include <optional>
#include <cassert>

#include <stencil_cpu.hpp>

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
        FLUID_SIM_MATRIX_NNZ * width * height);
    assert(a != std::nullopt);
    A = std::move(a.value());

    std::optional s = BiCGSTABSolverCpu::create_shared(
        FLUID_SIM_EQ_TYPES * width * height);
    assert(s != std::nullopt);
    solver = std::move(s.value());

    RHS.resize(FLUID_SIM_EQ_TYPES * width * height);
}

void FluidSimulationBackendCPU::apply_user_input()
{
    for (int i = 0; i < width * height; i++)
    {
        u[i] += user_u[i];
        v[i] += user_v[i];
        rho[i] += user_rho[i];
    }
}

void FluidSimulationBackendCPU::build_CSR_matrix()
{
    double dt;
    // clang-format off
    StencilCpu u_type_stencils[] = {
        {.offset = U_OFFSET, .di = 0, .dj = 0, 
            .get_val = [&](int i, int j)
            {
                return 1 / dt + (mu * (2 / (dx_sq)) + 2 / (dy_sq) + (mu + lambda) * 2 / dx_sq) / solver->x[RHO_OFFSET + AT(i, j)];
            }
        },
        {.offset = U_OFFSET, .di = -1, .dj = 0, 
            .get_val = [&](int i, int j)
            {
                return -(mu / dx_sq + (mu + lambda) / dx_sq) / solver->x[RHO_OFFSET + AT(i, j)] - solver->x[U_OFFSET + AT(i,j)] / (2 * dx);
            }
        },
        {.offset = U_OFFSET, .di = 1, .dj = 0, 
            .get_val = [&](int i, int j)
            {
                return -(mu / dx_sq + (mu + lambda) / dx_sq) / solver->x[RHO_OFFSET + AT(i, j)] + solver->x[U_OFFSET + AT(i,j)] / (2 * dx);
            }
        },
        {.offset = U_OFFSET, .di = 0, .dj = -1, 
            .get_val = [&](int i, int j)
            {
                return -(mu / dy_sq) / solver->x[RHO_OFFSET + AT(i, j)] - solver->x[V_OFFSET + AT(i,j)] / (2 * dy);
            }
        },
        {.offset = U_OFFSET, .di = 0, .dj = 1, 
            .get_val = [&](int i, int j)
            {
                return -(mu / dy_sq) / solver->x[RHO_OFFSET + AT(i, j)] + solver->x[V_OFFSET + AT(i,j)] / (2 * dy);
            }
        },
        {.offset = V_OFFSET, .di = 1, .dj = 1, 
            .get_val = [&](int i, int j)
            {
                return -1 / (4 * solver->x[RHO_OFFSET + AT(i, j)] * dx * dy);
            }
        },
        {.offset = V_OFFSET, .di = -1, .dj = -1, 
            .get_val = [&](int i, int j)
            {
                return -1 / (4 * solver->x[RHO_OFFSET + AT(i, j)] * dx * dy);
            }
        },
        {.offset = V_OFFSET, .di = -1, .dj = 1, 
            .get_val = [&](int i, int j)
            {
                return 1 / (4 * solver->x[RHO_OFFSET + AT(i, j)] * dx * dy);
            }
        },
        {.offset = V_OFFSET, .di = 1, .dj = -1, 
            .get_val = [&](int i, int j)
            {
                return 1 / (4 * solver->x[RHO_OFFSET + AT(i, j)] * dx * dy);
            }
        },
        {.offset = RHO_OFFSET, .di = 1, .dj = 0, 
            .get_val = [&](int i, int j)
            {
                return cs_sq / (2 * solver->x[RHO_OFFSET + AT(i, j)] * dx);
            }
        },
        {.offset = RHO_OFFSET, .di = -1, .dj = 0, 
            .get_val = [&](int i, int j)
            {
                return -cs_sq / (2 * solver->x[RHO_OFFSET + AT(i, j)] * dx);
            }
        }
    };
    // clang-format on

    int offset = 0;
    /* LHS1 equations */
    for (int j = 0; j < height; j++)
    {
        for (int i = 0; i < width; i++)
        {
            int idx = FLUID_SIM_LHS1_NNZ * (j * width + i);
            A->row_ptr[offset + j * width + i] = idx;
            /* Set row using u-type stencil */
            int p = 0;
            for(StencilCpu& s : u_type_stencils){
                A->col[idx + p] = s.offset + (j + s.dj) * width + (i + s.di);
                p++;
            }
        }
    }

    offset += width * height;
    /* LHS2 equations */
    for (int j = 0; j < height; j++)
    {
        for (int i = 0; i < width; i++)
        {
            A->row_ptr[offset + i * j] = FLUID_SIM_LHS2_NNZ * i * j;
            /* Cols of U components */

            /* Cols of V components */

            /* Cols of RHO components */
        }
    }

    offset += width * height;
    /* LHS3 equations */
    for (int j = 0; j < height; j++)
    {
        for (int i = 0; i < width; i++)
        {
            A->row_ptr[offset + i * j] = FLUID_SIM_LHS3_NNZ * i * j;
            /* Cols of U components */

            /* Cols of V components */

            /* Cols of RHO components */
        }
    }
}
