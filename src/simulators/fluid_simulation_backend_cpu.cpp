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

    is_wall.resize(width * height);
    RHS.resize(FLUID_SIM_EQ_TYPES * width * height);
    init_stencils();
}

void FluidSimulationBackendCPU::init_stencils()
{
    // clang-format off
    /* Create u-type stencil for LHS1 equations */
    u_type_stencil = {
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

    /* Create v-type stencil for LHS2 equations */
    v_type_stencil = {
        {.offset = V_OFFSET, .di = 0, .dj = 0, 
            .get_val = [&](int i, int j)
            {
                return 1 / dt + (mu * (-2 / (dx_sq)) - 2 / (dy_sq) - (mu + lambda) * 2 / dx_sq) / solver->x[RHO_OFFSET + AT(i, j)];
            }
        },
        {.offset = V_OFFSET, .di = -1, .dj = 0, 
            .get_val = [&](int i, int j)
            {
                return -(mu / (dx_sq)) / solver->x[RHO_OFFSET + AT(i, j)] + solver->x[U_OFFSET + AT(i,j)] / (2 * dx);
            }
        },
        {.offset = V_OFFSET, .di = 1, .dj = 0, 
            .get_val = [&](int i, int j)
            {
                return -(mu / (dx_sq)) / solver->x[RHO_OFFSET + AT(i, j)] - solver->x[U_OFFSET + AT(i,j)] / (2 * dx);
            }
        },
        {.offset = V_OFFSET, .di = 0, .dj = -1, 
            .get_val = [&](int i, int j)
            {
                return -((2 * mu + lambda) / dy_sq) / solver->x[RHO_OFFSET + AT(i, j)] + solver->x[V_OFFSET + AT(i,j)] / (2 * dy);
            }
        },
        {.offset = V_OFFSET, .di = 0, .dj = 1, 
            .get_val = [&](int i, int j)
            {
                return -((2 * mu + lambda) / dy_sq) / solver->x[RHO_OFFSET + AT(i, j)] - solver->x[V_OFFSET + AT(i,j)] / (2 * dy);
            }
        },
        {.offset = U_OFFSET, .di = 1, .dj = 1, 
            .get_val = [&](int i, int j)
            {
                return -1 / ( 4 * solver->x[RHO_OFFSET + AT(i, j)] * dx * dy);
            }
        },
        {.offset = U_OFFSET, .di = -1, .dj = -1, 
            .get_val = [&](int i, int j)
            {
                return -1 / ( 4 * solver->x[RHO_OFFSET + AT(i, j)] * dx * dy);
            }
        },
        {.offset = U_OFFSET, .di = -1, .dj = 1, 
            .get_val = [&](int i, int j)
            {
                return 1 / ( 4 * solver->x[RHO_OFFSET + AT(i, j)] * dx * dy);
            }
        },
        {.offset = U_OFFSET, .di = 1, .dj = -1, 
            .get_val = [&](int i, int j)
            {
                return -1 / ( 4 * solver->x[RHO_OFFSET + AT(i, j)] * dx * dy);
            }
        },
        {.offset = RHO_OFFSET, .di = 0, .dj = 1, 
            .get_val = [&](int i, int j)
            {
                return cs_sq / (2 * solver->x[RHO_OFFSET + AT(i, j)] * dy);
            }
        },
        {.offset = RHO_OFFSET, .di = 1, .dj = 0, 
            .get_val = [&](int i, int j)
            {
                return -cs_sq / (2 * solver->x[RHO_OFFSET + AT(i, j)] * dy);
            }
        }
    };

    rho_type_stencil = {
        {.offset = RHO_OFFSET, .di = 0, .dj = 0, 
            .get_val = [&](int i, int j)
            {
                return 1 / dt + D *(2 / dx_sq + 2 / dy_sq);
            }
        },
        {.offset = RHO_OFFSET, .di = -1, .dj = 0, 
            .get_val = [&](int i, int j)
            {
                return -D / dx_sq;
            }
        },
        {.offset = RHO_OFFSET, .di = 1, .dj = 0, 
            .get_val = [&](int i, int j)
            {
                return -D / dx_sq;
            }
        },
        {.offset = RHO_OFFSET, .di = 0, .dj = -1, 
            .get_val = [&](int i, int j)
            {
                return -D / dy_sq;
            }
        },
        {.offset = RHO_OFFSET, .di = 1, .dj = 1, 
            .get_val = [&](int i, int j)
            {
                return -D / dy_sq;
            }
        }
    };

    /* Wall stencils */
    wall_u_type_stencil =   {
        {.offset = U_OFFSET, .di = 0, .dj = 0, 
            .get_val = [&](int i, int j)
            {
                return 1.0;
            }
        }
    };

    wall_v_type_stencil =  {
        {.offset = V_OFFSET, .di = 0, .dj = 0, 
            .get_val = [&](int i, int j)
            {
                return 1.0;
            }
        }
    };

    /* Rho wall stencil needs 4 neighbours, and calc non-walls */
    wall_rho_type_stencil =  {
        {.offset = RHO_OFFSET, .di = 0, .dj = 0, 
            .get_val = [&](int i, int j)
            {
                int fluid_neighbors = 0;
                if(!is_wall[AT(i + 1, j)]) fluid_neighbors++;
                if(!is_wall[AT(i - 1, j)]) fluid_neighbors++;
                if(!is_wall[AT(i, j + 1)]) fluid_neighbors++;
                if(!is_wall[AT(i, j - 1)]) fluid_neighbors++;
                return fluid_neighbors;
            }
        },
        {.offset = RHO_OFFSET, .di = 1, .dj = 0, 
            .get_val = [&](int i, int j)
            {
                return (is_wall[AT(i + 1, j)]) ? 0.0f : 1.0f;
            }
        },
        {.offset = RHO_OFFSET, .di = -1, .dj = 0, 
            .get_val = [&](int i, int j)
            {
                return (is_wall[AT(i - 1, j)]) ? 0.0f : 1.0f;
            }
        },
        {.offset = RHO_OFFSET, .di = 0, .dj = 1, 
            .get_val = [&](int i, int j)
            {
                return (is_wall[AT(i, j + 1)]) ? 0.0f : 1.0f;
            }
        },
        {.offset = RHO_OFFSET, .di = 0, .dj = -1, 
            .get_val = [&](int i, int j)
            {
                return (is_wall[AT(i, j - 1)]) ? 0.0f : 1.0f;
            }
        }
    };
    // clang-format on
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
    /* !!! TODO !!!
        There is absolutley no need to loop this many times,
    since each stencil works on private part of matrix we can apply
    all types of stencils at the same time. We also iterate over all
    cells, so we rither fill out fluid cell or solid wall cell, so logic
    is in the same place. SIMPLY FILL OUT CELL BY CELL EITHER WALL OR FLUID
    */
    int offset = 0;
    /* LHS1 equations */
    for (int j = 0; j < height; j++)
    {
        for (int i = 0; i < width; i++)
        {
            int idx = FLUID_SIM_LHS1_NNZ * (j * width + i);
            A->row_ptr[offset + j * width + i] = idx;
            /* Set row using u-type stencil */
            for (int p = 0; p < FLUID_SIM_LHS1_NNZ; p++)
            {
                A->col[idx + p] = u_type_stencil[p].offset + (j + u_type_stencil[p].dj) * width + (i + u_type_stencil[p].di);
                A->val[idx + p] = *u_type_stencil[p].get_val(i, j);
            }
        }
    }

    offset += width * height;
    /* LHS2 equations */
    for (int j = 0; j < height; j++)
    {
        for (int i = 0; i < width; i++)
        {
            int idx = FLUID_SIM_LHS2_NNZ * (j * width + i);
            A->row_ptr[offset + j * width + i] = idx;
            /* Set row using u-type stencil */
            for (int p = 0; p < FLUID_SIM_LHS2_NNZ; p++)
            {
                A->col[idx + p] = v_type_stencil[p].offset + (j + v_type_stencil[p].dj) * width + (i + v_type_stencil[p].di);
                A->val[idx + p] = v_type_stencil[p].get_val(i, j);
            }
        }
    }

    offset += width * height;
    /* LHS3 equations */
    for (int j = 0; j < height; j++)
    {
        for (int i = 0; i < width; i++)
        {
            int idx = FLUID_SIM_LHS3_NNZ * (j * width + i);
            A->row_ptr[offset + j * width + i] = idx;
            /* Set row using u-type stencil */
            for (int p = 0; p < FLUID_SIM_LHS3_NNZ; p++)
            {
                A->col[idx + p] = rho_type_stencil[p].offset + (j + rho_type_stencil[p].dj) * width + (i + rho_type_stencil[p].di);
                A->val[idx + p] = rho_type_stencil[p].get_val(i, j);
            }
        }
    }

    /* Apply wall constriants
    (since we are sure that wall stencils will not enforce
    any non-zeros that could be used previously we still have
    keep estimates on NNZs and their layout in CSR matrix)
    */
    /*
    For each wall cell:
    1. Zero out row for given row
    2. Apply appropriate stencil
    3. Set appropriate RHS
    */
    for (int j = 0; j < height; j++)
    {
        for (int i = 0; i < width; i++)
        {
            if (!is_wall[AT(i, j)])
                continue;
            std::fill(A->val[U_OFFSET + AT(i, j)], A->val[U_OFFSET + AT(i, j) + FLUID_SIM_LHS1_NNZ], 0.0);
            std::fill(A->val[V_OFFSET + AT(i, j)], A->val[U_OFFSET + AT(i, j) + FLUID_SIM_LHS2_NNZ], 0.0);
            std::fill(A->val[RHO_OFFSET + AT(i, j)], A->val[U_OFFSET + AT(i, j) + FLUID_SIM_LHS3_NNZ], 0.0);

            int p = 0;
            for (StencilCpu &s : wall_u_type_stencil)
            {
                A->col[U_OFFSET + p] = s.offset + (j + s.dj) * width + (i + s.di);
                A->val[U_OFFSET + p] = s.get_val(i, j);
                p++;
            }

            p = 0;
            for (StencilCpu &s : wall_v_type_stencil)
            {
                A->col[V_OFFSET + p] = s.offset + (j + s.dj) * width + (i + s.di);
                A->val[V_OFFSET + p] = s.get_val(i, j);
                p++;
            }

            p = 0;
            for (StencilCpu &s : wall_rho_type_stencil)
            {
                A->col[RHO_OFFSET + p] = s.offset + (j + s.dj) * width + (i + s.di);
                A->val[RHO_OFFSET + p] = s.get_val(i, j);
                p++;
            }
        }
    }
}

void FluidSimulationBackendCPU::build_CSR_matrix_2()
{
    for (int j = 0; j < height; j++)
    {
        for (int i = 0; i < width; i++)
        {
            int off1 = width * height * FLUID_SIM_LHS1_NNZ;
            int off2 = 2 * off1;

            if (is_wall[AT(i, j)])
            {
                std::fill(
                    A->val[AT(i, j) * FLUID_SIM_LHS1_NNZ],
                    A->val[(AT(i, j) + 1) * FLUID_SIM_LHS1_NNZ], 0.0);
                std::fill(
                    A->val[off1 + AT(i, j) * FLUID_SIM_LHS2_NNZ],
                    A->val[off1 + (AT(i, j) + 1) * FLUID_SIM_LHS2_NNZ], 0.0);
                std::fill(
                    A->val[off2 + AT(i, j) * FLUID_SIM_LHS1_NNZ],
                    A->val[off2 + (AT(i, j) + 1) * FLUID_SIM_LHS1_NNZ], 0.0);

                std::fill(
                    A->col[AT(i, j) * FLUID_SIM_LHS1_NNZ],
                    A->col[(AT(i, j) + 1) * FLUID_SIM_LHS1_NNZ], 0);
                std::fill(
                    A->col[off1 + AT(i, j) * FLUID_SIM_LHS2_NNZ],
                    A->col[off1 + (AT(i, j) + 1) * FLUID_SIM_LHS2_NNZ], 0);
                std::fill(
                    A->col[off2 + AT(i, j) * FLUID_SIM_LHS1_NNZ],
                    A->col[off2 + (AT(i, j) + 1) * FLUID_SIM_LHS1_NNZ], 0);

                int p = 0;
                for (StencilCpu &s : wall_u_type_stencil)
                {
                    A->col[AT(i, j) * FLUID_SIM_LHS1_NNZ + p] = s.offset + (j + s.dj) * width + (i + s.di);
                    A->val[AT(i, j) * FLUID_SIM_LHS1_NNZ + p] = s.get_val(i, j);
                    p++;
                }

                p = 0;
                for (StencilCpu &s : wall_v_type_stencil)
                {
                    int off = width * height * FLUID_SIM_LHS1_NNZ;
                    A->col[off1 + AT(i, j) * FLUID_SIM_LHS2_NNZ + p] = s.offset + (j + s.dj) * width + (i + s.di);
                    A->val[off1 + AT(i, j) * FLUID_SIM_LHS2_NNZ + p] = s.get_val(i, j);
                    p++;
                }

                p = 0;
                for (StencilCpu &s : wall_rho_type_stencil)
                {
                    A->col[off2 + AT(i, j) * FLUID_SIM_LHS3_NNZ + p] = s.offset + (j + s.dj) * width + (i + s.di);
                    A->val[off2 + AT(i, j) * FLUID_SIM_LHS3_NNZ + p] = s.get_val(i, j);
                    p++;
                }
            }
            else
            {
                /* u-type equation : LHS1 */

                /* v_type equation : LHS2 */

                /* rho_type equation :LHS3 */
                
            }
        }
    }
}

void FluidSimulationBackendCPU::run_BiCSTAB()
{
    solver->solve(*A.get(), RHS, 50, 1e-10);
}
