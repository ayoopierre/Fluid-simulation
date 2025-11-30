# Fluid simulation code
This is monolithic Navier-Stokes equation sovler

# Buliding
g++ -o main.exe src/main.cpp src/solvers/* src/simulators/* src/utils/* -Iinc/solvers -Iinc/utils -Iinc/simulators