#ifndef MY_FLUID_SIMULATION_BACKEND
#define MY_FLUID_SIMULATION_BACKEND

#include <vector>

class FluidSimulationBackend{
    public:
        virtual void step(double dt) = 0;

    protected:
        virtual void apply_user_input() = 0;
        virtual void build_CSR_matrix() = 0;
        virtual void run_BiCSTAB() = 0;
        virtual void update_pressure() = 0;
        virtual void write_heatmap() = 0;

        int width;
        int height;

        std::vector<double> u;
        std::vector<double> v;
        std::vector<double> rho;

        std::vector<double> user_u;
        std::vector<double> user_v;
        std::vector<double> user_rho;
};

#endif