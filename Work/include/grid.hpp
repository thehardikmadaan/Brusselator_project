#ifndef GRID_HPP
#define GRID_HPP

#include <string>

// Structure to store all grid, physical, and solver parameters
struct GridInfo {
    int nx = 40;                 // Number of grid points in x-direction
    int ny = 40;                 // Number of grid points in y-direction
    double dx = 0.005;            // Grid spacing in x-direction [m]
    double dy = 0.005;            // Grid spacing in y-direction [m]
    double duration = 30.0;       // Total simulation time [s]
    double dt = 0.01;             // Time step size [s]
    double write_interval = 1.0;  // Time interval between writing output files [s]
    double Ca = 1.0;              // Brusselator parameter Ca
    double Cb = 3.0;              // Brusselator parameter Cb
    double D1 = 1e-5;             // Diffusion coefficient for species 1 (activator) [m^2/s]
    double D2 = 1e-6;             // Diffusion coefficient for species 2 (inhibitor) [m^2/s]
    std::string method = "rk4";   // Time integration method ("rk4" or "implicit")
    bool reactions = true;        // Toggle for chemical reactions (on/off)
};

#endif // GRID_HPP
