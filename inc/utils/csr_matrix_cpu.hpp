#ifndef MY_CSR_MATRIX_CPU
#define MY_CSR_MATRIX_CPU

#include <memory>
#include <vector>
#include <optional>

class CSRMatrixCPU{
    public:
        /*
        For note why there is no public constructor check BiCGSTAB solver note
        */
        static std::optional<CSRMatrixCPU> create(int rows, int nnz);
        static std::optional<std::shared_ptr<CSRMatrixCPU>> create_shared(int rows, int nnz);

        /* We need custom assigment to avoid copies */
        CSRMatrixCPU(CSRMatrixCPU& other);
        CSRMatrixCPU& operator=(CSRMatrixCPU& other);

        CSRMatrixCPU(CSRMatrixCPU&& other);
        CSRMatrixCPU& operator=(CSRMatrixCPU&& other);

        double get_val(int i, int j);
        
        int nnz;
        int rows;

        std::vector<int> row_ptr;
        std::vector<int> col;
        std::vector<double> val;
    private:
        CSRMatrixCPU() {};
};

#endif