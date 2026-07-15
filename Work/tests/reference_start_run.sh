#!/bin/bash
# This script runs a single simulation of the 2D Brusselator reaction-diffusion system.
# Positional Arguments:
#   1. nx             - Number of grid points in x-direction (e.g., 40)
#   2. ny             - Number of grid points in y-direction (e.g., 40)
#   3. duration       - Total simulation time in seconds (e.g., 5.0)
#   4. dt             - Time step size in seconds (e.g., 0.01)
#   5. write_interval - Time interval between writing output files in seconds (e.g., 0.1)
#   6. Ca             - Brusselator parameter Ca (e.g., 1.0)
#   7. Cb             - Brusselator parameter Cb (e.g., 3.0)
#   8. method         - Time integration method (rk4 or implicit) (e.g., implicit)
#   9. reactions      - Reaction toggle (on or off) (e.g., on)

./bin/brusselator 40 40 5.0 0.01 0.1 1.0 3.0 implicit on
