#include "csr_matrix_cpu.hpp"

std::optional<CSRMatrixCPU> CSRMatrixCPU::create(int rows, int nnz)
{
    CSRMatrixCPU instance;
    try
    {
        instance.row_ptr.resize(rows + 1);
        instance.col.resize(nnz);
        instance.val.resize(nnz);

        instance.rows = rows + 1;
        instance.nnz = nnz;
    }
    catch (const std::exception &e)
    {
        /*
        If reallocation for vector fails std::bad_alloc is thrown,
        catch an do not return instance of an object.
        */
        return std::nullopt;
    }
    return std::move(instance);
}

std::optional<std::shared_ptr<CSRMatrixCPU>> CSRMatrixCPU::create_shared(int rows, int nnz)
{
    std::shared_ptr<CSRMatrixCPU> instance;
    try
    {
        instance->row_ptr.resize(rows + 1);
        instance->col.resize(nnz);
        instance->val.resize(nnz);

        instance->rows = rows;
        instance->nnz = nnz;
    }
    catch (const std::exception &e)
    {
        /*
        If reallocation for vector fails std::bad_alloc is thrown,
        catch an do not return instance of an object.
        */
        return std::nullopt;
    }
    return instance;
}

CSRMatrixCPU::CSRMatrixCPU(CSRMatrixCPU &other)
{
    nnz = other.nnz;
    rows = other.rows;
    col = other.col;
    val = other.val;
    row_ptr = other.row_ptr;
}

CSRMatrixCPU &CSRMatrixCPU::operator=(CSRMatrixCPU &other)
{
    nnz = other.nnz;
    rows = other.rows;
    col = other.col;
    val = other.val;
    row_ptr = other.row_ptr;

    return *this;
}

CSRMatrixCPU::CSRMatrixCPU(CSRMatrixCPU &&other)
{
    nnz = other.nnz;
    rows = other.rows;
    col = std::move(other.col);
    val = std::move(other.val);
    row_ptr = std::move(other.row_ptr);
}

CSRMatrixCPU &CSRMatrixCPU::operator=(CSRMatrixCPU &&other)
{
    nnz = other.nnz;
    rows = other.rows;
    col = std::move(other.col);
    val = std::move(other.val);
    row_ptr = std::move(other.row_ptr);

    return *this;
}

double CSRMatrixCPU::get_val(int i, int j)
{
    if (i < 0 || i >= rows)
        return 0.0;

    int start = row_ptr[i];
    int end = row_ptr[i + 1];

    if (end - start < 16)
    {
        for (int k = start; k < end; k++)
            if (col[k] == j)
                return val[k];
        return 0.0;
    }

    while (start < end)
    {
        int mid = start + (end - start) / 2;

        if (col[mid] < j)
        {
            start = mid + 1;
        }
        else
        {
            end = mid;
        }
    }

    return (start < row_ptr[i + 1] && col[start] == j) ? val[start] : 0.0;
}
