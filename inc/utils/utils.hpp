#ifndef MY_UTILS
#define MY_UTILS

#include <cstring>
#include <vector>
#include <algorithm>
#include <cstdio>

#include "csr_matrix_cpu.hpp"

inline void vec_plus_vec(std::vector<double> &a,
                         std::vector<double> &b,
                         std::vector<double> &c)
{
    for (int i = 0; i < a.size(); i++)
        c[i] = a[i] + b[i];
}

inline void vec_scale(double scale, std::vector<double> &a)
{
    for (int i = 0; i < a.size(); i++)
        a[i] *= scale;
}

inline double vec_dot_vec(
    std::vector<double> &a,
    std::vector<double> &b)
{
    double dot = 0.0f;

    for (int i = 0; i < a.size(); i++)
        dot += a[i] * b[i];

    return dot;
}

inline void copy_a_to_b(std::vector<double> &a, std::vector<double> &b)
{
    std::memcpy(b.data(), a.data(), sizeof(double) * a.size());
}

inline int set_one_on_first_non_zero(std::vector<double> &r,
                                     std::vector<double> &r_hat)
{
    std::fill(r_hat.begin(), r_hat.end(), 0.0f);
    for (int i = 0; i < r.size(); i++)
    {
        if (r[i] != 0)
        {
            r_hat[i] = 1.0f;
            return i;
        }
    }
    return -1;
}

inline void vec_CSRmat_prod(CSRMatrixCPU &A,
                            std::vector<double> &x,
                            std::vector<double> &b)
{
    std::fill(b.begin(), b.end(), 0.0f);
    for (int row_index = 0; row_index < A.rows; row_index++)
    {
        for (int col_index = A.row_ptr[row_index]; col_index < A.row_ptr[row_index + 1]; col_index++)
        {
            b[row_index] += A.val[col_index] * x[A.col[col_index]];
        }
    }
}

#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

inline bool writeCSRtoCSV(
    const std::string& filename,
    const std::vector<double>& values,
    const std::vector<int>& col_idx,
    const std::vector<int>& row_ptr,
    int rows,
    int cols)
{
    std::ofstream out(filename);
    if (!out.is_open())
        return false;

    // For each row
    for (int r = 0; r < rows; ++r)
    {
        int start = row_ptr[r];
        int end   = row_ptr[r + 1];
        int nz    = start;  // pointer to next non-zero in this row

        // For each column
        for (int c = 0; c < cols; ++c)
        {
            auto b_it = col_idx.begin() + start;
            auto e_it =  col_idx.begin() + end;

            auto it = std::find(b_it, e_it, c);
            if(it != e_it){
                auto v_it = values.begin() + std::distance(col_idx.begin(), it);
                out << (*v_it);
            }
            else{
                out << 0.0;
            }
            if(c + 1 < cols)
                out << ",";
        }
        out << "\n";
    }

    return true;
}

#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

inline bool writeStateToCSV(
    const std::string& filename,
    const std::vector<double>& values,
    int rows,
    int cols, int offset)
{
    std::ofstream out(filename);
    if (!out.is_open())
        return false;

    // For each row
    for (int r = 0; r < rows; ++r)
    {
        // For each column
        for (int c = 0; c < cols; ++c)
        {
            out << values[offset + r * rows + c];
            if(c + 1 < cols)
                out << ",";
        }
        out << "\n";
    }

    return true;
}


#endif