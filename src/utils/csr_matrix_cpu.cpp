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
    catch(const std::exception& e)
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
    catch(const std::exception& e)
    {
        /*
        If reallocation for vector fails std::bad_alloc is thrown,
        catch an do not return instance of an object.
        */
        return std::nullopt;
    }
    return instance;
}
