#include "fluid_simulation_backend_cpu.hpp"

#include <optional>
#include <cassert>
#include <cstdio>
#include <cmath>

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

#ifdef VERBOSE
    std::printf("Managed to resize needed vectors\n");
#endif
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
    if (a.has_value())
    {
        A = std::move(a.value());
#ifdef VERBOSE
        std::printf("Allocated CSR matrix\n");
#endif
    }

    std::optional s = BiCGSTABSolverCpu::create_shared(
        FLUID_SIM_EQ_TYPES * width * height);
    if (s.has_value())
    {
        solver = std::move(s.value());
#ifdef VERBOSE
        std::printf("Allocated BiCGSTAB solver\n");
#endif
    }

    is_wall.resize(width * height);

    for (int i = 0; i < width; i++)
    {
        is_wall[i] = true;
        is_wall[(height - 1) * width + i] = true;
    }
#ifdef VERBOSE
    std::printf("Set horizontal walls\n");
#endif

    for (int j = 0; j < height; j++)
    {
        is_wall[j * width] = true;
        is_wall[j * width + width - 1] = true;
    }
#ifdef VERBOSE
    std::printf("Set vertical walls\n");
#endif

#ifdef VERBOSE
    std::printf("Allocated wall bitmap\n");
#endif

    RHS.resize(FLUID_SIM_EQ_TYPES * width * height);

#ifdef VERBOSE
    std::printf("Allocated RHS buffer\n");
#endif

    init_stencils();

#ifdef VERBOSE
    std::printf("Created stencils\n");
#endif
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
                double u_tn = solver->x[U_OFFSET + AT(i, j)];
                double v_tn = solver->x[V_OFFSET + AT(i, j)];

                double max_u_tn = std::max(u_tn, 0.0);
                double min_u_tn = std::min(u_tn, 0.0);

                double max_v_tn = std::max(v_tn, 0.0);
                double min_v_tn = std::min(v_tn, 0.0);

                return 1 / dt + (max_u_tn - min_u_tn) / dx + (max_v_tn - min_v_tn) / dy;
            }
        },
        {.offset = RHO_OFFSET, .di = -1, .dj = 0, 
            .get_val = [&](int i, int j)
            {
                double u_tn = solver->x[U_OFFSET + AT(i, j)];
                return -std::max(u_tn, 0.0) / dx;
            }
        },
        {.offset = RHO_OFFSET, .di = 1, .dj = 0, 
            .get_val = [&](int i, int j)
            {
                double u_tn = solver->x[U_OFFSET + AT(i, j)];
                return std::min(u_tn, 0.0) / dx;
            }
        },
        {.offset = RHO_OFFSET, .di = 0, .dj = -1, 
            .get_val = [&](int i, int j)
            {
                double v_tn = solver->x[V_OFFSET + AT(i, j)];
                return -std::max(v_tn, 0.0) / dy;
            }
        },
        {.offset = RHO_OFFSET, .di = 0, .dj = 1, 
            .get_val = [&](int i, int j)
            {
                double v_tn = solver->x[V_OFFSET + AT(i, j)];
                return std::min(v_tn, 0.0) / dy;
            }
        },
        {.offset = U_OFFSET, .di = -1, .dj = 0, 
            .get_val = [&](int i, int j)
            {
                return -solver->x[RHO_OFFSET + AT(i, j)] / (2 * dx);
            }
        },
        {.offset = U_OFFSET, .di = 1, .dj = 0, 
            .get_val = [&](int i, int j)
            {
                return solver->x[RHO_OFFSET + AT(i, j)] / (2 * dx);
            }
        },
        {.offset = V_OFFSET, .di = 0, .dj = -1, 
            .get_val = [&](int i, int j)
            {
                return -solver->x[RHO_OFFSET + AT(i, j)] / (2 * dy);
            }
        },
        {.offset = V_OFFSET, .di = 0, .dj = 1, 
            .get_val = [&](int i, int j)
            {
                return solver->x[RHO_OFFSET + AT(i, j)] / (2 * dy);
            }
        },
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
                if(0 < AT(i + 1, j) && AT(i + 1, j) < width * height && !is_wall[AT(i + 1, j)]) fluid_neighbors++;
                if(0 < AT(i - 1, j) && AT(i - 1, j) < width * height && !is_wall[AT(i - 1, j)]) fluid_neighbors++;
                if(0 < AT(i, j + 1) && AT(i, j + 1) < width * height && !is_wall[AT(i, j + 1)]) fluid_neighbors++;
                if(0 < AT(i, j - 1) && AT(i, j - 1) < width * height && !is_wall[AT(i, j - 1)]) fluid_neighbors++;
                return fluid_neighbors != 0 ? (double)fluid_neighbors : 1.0;
            }
        },
        {.offset = RHO_OFFSET, .di = 1, .dj = 0, 
            .get_val = [&](int i, int j)
            {
                return ( is_wall[AT(i + 1, j)]) ? 0.0 : -1.0;
            }
        },
        {.offset = RHO_OFFSET, .di = -1, .dj = 0, 
            .get_val = [&](int i, int j)
            {
                return (is_wall[AT(i - 1, j)]) ? 0.0 : -1.0;
            }
        },
        {.offset = RHO_OFFSET, .di = 0, .dj = 1, 
            .get_val = [&](int i, int j)
            {
                return (is_wall[AT(i, j + 1)]) ? 0.0 : -1.0;
            }
        },
        {.offset = RHO_OFFSET, .di = 0, .dj = -1, 
            .get_val = [&](int i, int j)
            {
                return (is_wall[AT(i, j - 1)]) ? 0.0 : -1.0;
            }
        }
    };
    // clang-format on

    auto lbd = [this](StencilCpu &s1, StencilCpu &s2)
    {
        return s1.offset + AT(s1.di, s1.dj) < s2.offset + AT(s2.di, s2.dj);
    };

    std::sort(u_type_stencil.begin(), u_type_stencil.end(), lbd);
    std::sort(v_type_stencil.begin(), v_type_stencil.end(), lbd);
    std::sort(rho_type_stencil.begin(), rho_type_stencil.end(), lbd);

    std::sort(wall_u_type_stencil.begin(), wall_u_type_stencil.end(), lbd);
    std::sort(wall_v_type_stencil.begin(), wall_v_type_stencil.end(), lbd);
    std::sort(wall_rho_type_stencil.begin(), wall_rho_type_stencil.end(), lbd);
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
    A->row_ptr[A->rows] = A->nnz;

    for (int j = 0; j < height; j++)
    {
        for (int i = 0; i < width; i++)
        {
            int off1 = width * height * FLUID_SIM_LHS1_NNZ;
            int off2 = off1 + width * height * FLUID_SIM_LHS2_NNZ;
            int p = 0;

            /* Set row pointers */
            A->row_ptr[U_OFFSET + AT(i, j)] = AT(i, j) * FLUID_SIM_LHS1_NNZ;

            A->row_ptr[V_OFFSET + AT(i, j)] = V_OFFSET * FLUID_SIM_LHS1_NNZ +
                                              AT(i, j) * FLUID_SIM_LHS2_NNZ;

            /* We know LHS1 and LHS2 has same NNZ */
            A->row_ptr[RHO_OFFSET + AT(i, j)] = RHO_OFFSET * FLUID_SIM_LHS1_NNZ +
                                                AT(i, j) * FLUID_SIM_LHS3_NNZ;

            // std::printf("i: %5d, j: %5d, U: %5d, V: %5d, RHO: %5d\n", i, j,
            //             A->row_ptr[U_OFFSET + AT(i, j)],
            //             A->row_ptr[V_OFFSET + AT(i, j)],
            //             A->row_ptr[RHO_OFFSET + AT(i, j)]);

            if (is_wall[AT(i, j)])
            {
                std::fill(
                    A->val.begin() + AT(i, j) * FLUID_SIM_LHS1_NNZ,
                    A->val.begin() + (AT(i, j) + 1) * FLUID_SIM_LHS1_NNZ, 0.0);
                std::fill(
                    A->val.begin() + off1 + AT(i, j) * FLUID_SIM_LHS2_NNZ,
                    A->val.begin() + off1 + (AT(i, j) + 1) * FLUID_SIM_LHS2_NNZ, 0.0);
                std::fill(
                    A->val.begin() + off2 + AT(i, j) * FLUID_SIM_LHS3_NNZ,
                    A->val.begin() + off2 + (AT(i, j) + 1) * FLUID_SIM_LHS3_NNZ, 0.0);

                std::fill(
                    A->col.begin() + AT(i, j) * FLUID_SIM_LHS1_NNZ,
                    A->col.begin() + (AT(i, j) + 1) * FLUID_SIM_LHS1_NNZ, 0);
                std::fill(
                    A->col.begin() + off1 + AT(i, j) * FLUID_SIM_LHS2_NNZ,
                    A->col.begin() + off1 + (AT(i, j) + 1) * FLUID_SIM_LHS2_NNZ, 0);
                std::fill(
                    A->col.begin() + off2 + AT(i, j) * FLUID_SIM_LHS3_NNZ,
                    A->col.begin() + off2 + (AT(i, j) + 1) * FLUID_SIM_LHS3_NNZ, 0);

                p = 0;
                for (StencilCpu &s : wall_u_type_stencil)
                {
                    SKIP_IF_NOT_IN_BOUNDS(i + s.di, j + s.dj);
                    A->col[A->row_ptr[U_OFFSET + AT(i, j)] + p] = s.offset + (j + s.dj) * width + (i + s.di);
                    A->val[A->row_ptr[U_OFFSET + AT(i, j)] + p] = s.get_val(i, j);
                    p++;
                }
                /* We can also set RHS1 */
                RHS[U_OFFSET + AT(i, j)] = 0;

                p = 0;
                for (StencilCpu &s : wall_v_type_stencil)
                {
                    SKIP_IF_NOT_IN_BOUNDS(i + s.di, j + s.dj);
                    A->col[A->row_ptr[V_OFFSET + AT(i, j)] + p] = s.offset + (j + s.dj) * width + (i + s.di);
                    A->val[A->row_ptr[V_OFFSET + AT(i, j)] + p] = s.get_val(i, j);
                    p++;
                }
                /* We can also set RHS2 */
                RHS[V_OFFSET + AT(i, j)] = 0;

                p = 0;
                for (StencilCpu &s : wall_rho_type_stencil)
                {
                    SKIP_IF_NOT_IN_BOUNDS(i + s.di, j + s.dj);
                    A->col[A->row_ptr[RHO_OFFSET + AT(i, j)] + p] = s.offset + (j + s.dj) * width + (i + s.di);
                    A->val[A->row_ptr[RHO_OFFSET + AT(i, j)] + p] = s.get_val(i, j);
                    p++;
                }
                /* We can also set RHS3 */
                RHS[RHO_OFFSET + AT(i, j)] = 0;
            }
            else
            {
                /* No need to fill with zeros, since we will alaredy fill out all NNZ fields for given row */
                /* u-type equation : LHS1 */
                p = 0;
                for (StencilCpu &s : u_type_stencil)
                {
                    SKIP_IF_NOT_IN_BOUNDS(i + s.di, j + s.dj);
                    A->col[A->row_ptr[U_OFFSET + AT(i, j)] + p] = s.offset + (j + s.dj) * width + (i + s.di);
                    A->val[A->row_ptr[U_OFFSET + AT(i, j)] + p] = s.get_val(i, j);
                    p++;
                }
                /* Set RHS1 */
                RHS[U_OFFSET + AT(i, j)] = solver->x[U_OFFSET + AT(i, j)] / dt;

                /* v_type equation : LHS2 */
                p = 0;
                for (StencilCpu &s : v_type_stencil)
                {
                    SKIP_IF_NOT_IN_BOUNDS(i + s.di, j + s.dj);
                    A->col[A->row_ptr[V_OFFSET + AT(i, j)] + p] = s.offset + (j + s.dj) * width + (i + s.di);
                    A->val[A->row_ptr[V_OFFSET + AT(i, j)] + p] = s.get_val(i, j);
                    p++;
                }
                /* Set RHS2 */
                RHS[V_OFFSET + AT(i, j)] = solver->x[V_OFFSET + AT(i, j)] / dt;

                /* rho_type equation : LHS3 */
                p = 0;
                for (StencilCpu &s : v_type_stencil)
                {
                    SKIP_IF_NOT_IN_BOUNDS(i + s.di, j + s.dj);
                    A->col[A->row_ptr[RHO_OFFSET + AT(i, j)] + p] = s.offset + (j + s.dj) * width + (i + s.di);
                    A->val[A->row_ptr[RHO_OFFSET + AT(i, j)] + p] = s.get_val(i, j);
                    p++;
                }
                /* Set RHS3 */
                RHS[RHO_OFFSET + AT(i, j)] = solver->x[RHO_OFFSET + AT(i, j)] / dt;
            }
        }
    }
}

void FluidSimulationBackendCPU::run_BiCSTAB()
{
    solver->solve(*A.get(), RHS, 50, 1e-10);
}
