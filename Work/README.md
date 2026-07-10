# Brusselator Reaction-Diffusion System Solver

This project implements a finite difference method for solving the two-dimensional Brusselator reaction-diffusion system. The solver is written in C++ using the Eigen library and supports both explicit Runge-Kutta and implicit Euler time integration methods. A Bash script is also provided for automated parameter variation.

## How to build your code

Ensure you have CMake (>=3.12) and the Eigen3 library installed.

1. Create a `build` directory:
   ```bash
   mkdir build && cd build
   ```
2. Configure the project with CMake:
   ```bash
   cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=.. ..
   ```
3. Compile and install:
   ```bash
   make -j$(nproc)
   make install
   ```
   The binaries (`brusselator` and `test_solver`) will be installed in the `bin/` directory.

## How to run a single simulation

You can run a single simulation using the compiled `brusselator` executable. It accepts arguments from the command line.
The output VTK files will be stored in an `output/` directory created in your current working directory.

**Usage:**
```bash
./bin/brusselator <nx> <ny> <duration> <dt> <write_interval> <Ca> <Cb> [method] [reactions]
```
- `nx`, `ny`: Number of grid points in x and y directions
- `duration`: Simulated time duration (s)
- `dt`: Time step (s)
- `write_interval`: Simulation time between writing output files (s)
- `Ca`, `Cb`: Parameters of the Brusselator ODE System
- `method` (Optional): `rk4` (default) or `implicit`
- `reactions` (Optional): `on` (default) or `off`

**Example:**
```bash
./bin/brusselator 40 40 30.0 0.01 1.0 1.0 3.0 rk4 on
```

## How to use your parameter bash script

The `scripts/run_variation.sh` script automates parameter variations. It generates a `start_run.sh` script for each run, executes the simulation, and organizes the output files.

**Predefined Variations:**
- `var1`: Uses Explicit RK4, Ca=1.0, varies Cb in {1.5, 2.0, 2.5, 3.0}
- `var2`: Uses Implicit Euler, Ca=1.0, Cb=3.0, varies dt in {0.01, 0.02, 0.03, 0.04, 0.05}

**Run Predefined:**
```bash
./scripts/run_variation.sh var1
./scripts/run_variation.sh var2
```

**Run Custom Variation:**
You can vary any specific parameter across multiple values:
```bash
./scripts/run_variation.sh <parameter_name> "<values>" <variation_name>
```
*Example:* `./scripts/run_variation.sh Cb "1.5 2.0 2.5 3.0" my_custom_var`

## How the results will be organized

When using the parameter variation script, results will be saved in the `results/` directory, nested by the variation name and the specific parameter value. Each run directory will contain the generated `start_run.sh`, the VTK output files (`*.vtk`), and the simulation console output (`simulation.log`).

Example structure:
```
results/
└── var1/
    ├── run_1.5/
    │   ├── start_run.sh
    │   ├── output_0001.vtk
    │   ├── output_0002.vtk
    │   └── simulation.log
    └── run_2.0/
        ├── start_run.sh
        ├── output_0001.vtk
        ├── output_0002.vtk
        └── simulation.log
```
