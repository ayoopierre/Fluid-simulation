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

        inline double get_val(int i, int j){
            if(i < 0 || i >= rows) return 0.0;

            int start = row_ptr[i];
            int end = row_ptr[i + 1];

            if(end - start < 16){
                for(int k = start; k < end; k++)
                    if(col[k] == j)
                        return val[k];
                return 0.0;
            }

            while(start < end){
                int mid = start + (end - start) / 2;

                if(col[mid] < j){
                    start = mid + 1;
                }
                else{
                    end = mid;
                }
            }

            return (start < row_ptr[i + 1] && col[start] == j) ? val[start] : 0.0;
        }

        int nnz;
        int rows;

        std::vector<int> row_ptr;
        std::vector<int> col;
        std::vector<double> val;
    private:
        CSRMatrixCPU();
};

#endif