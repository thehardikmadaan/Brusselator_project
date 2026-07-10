#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Predefined Variation 1 (var1)
# Explicit RK4, Ca=1.0, vary Cb in {1.5, 2.0, 2.5, 3.0}, grid 40x40, duration 30s, dt=0.01, write_interval=1.0
run_var1() {
    VAR_NAME="var1"
    VALUES=(1.5 2.0 2.5 3.0)
    echo "Running Parameter Variation 1 (var1)..."
    
    for cb in "${VALUES[@]}"; do
        RUN_DIR="results/${VAR_NAME}/run_${cb}"
        mkdir -p "${RUN_DIR}"
        
        # Define simulation arguments
        NX=40
        NY=40
        DUR=30.0
        DT=0.01
        W_INT=1.0
        CA=1.0
        CB=$cb
        METHOD="rk4"
        REACT="on"
        
        # Create start_run.sh with mandatory explanatory comments in the root dir
        cat <<EOF > "start_run.sh"
#!/bin/bash
# This script runs a single simulation of the 2D Brusselator reaction-diffusion system.
# Positional Arguments:
#   1. nx             - Number of grid points in x-direction (e.g., $NX)
#   2. ny             - Number of grid points in y-direction (e.g., $NY)
#   3. duration       - Total simulation time in seconds (e.g., $DUR)
#   4. dt             - Time step size in seconds (e.g., $DT)
#   5. write_interval - Time interval between writing output files in seconds (e.g., $W_INT)
#   6. Ca             - Brusselator parameter Ca (e.g., $CA)
#   7. Cb             - Brusselator parameter Cb (e.g., $CB)
#   8. method         - Time integration method (rk4 or implicit) (e.g., $METHOD)
#   9. reactions      - Reaction toggle (on or off) (e.g., $REACT)

./bin/brusselator $NX $NY $DUR $DT $W_INT $CA $CB $METHOD $REACT
EOF
        chmod +x "start_run.sh"
        
        echo "  Running simulation for Cb = ${cb}..."
        # Clean output directory before running
        rm -rf output && mkdir -p output
        
        # Run the start script
        ./start_run.sh > "simulation.log" 2>&1
        
        # Move generated files to the run directory
        mv start_run.sh "${RUN_DIR}/"
        mv simulation.log "${RUN_DIR}/"
        mv output/*.vts "${RUN_DIR}/"
        rm -rf output
    done
    echo "Parameter Variation 1 completed. Results saved in results/var1/"
}

# Predefined Variation 2 (var2)
# Implicit Euler, Ca=1.0, Cb=3.0, vary dt in {0.01, 0.02, 0.03, 0.04, 0.05}, grid 40x40, duration 5s, write_interval=0.1
run_var2() {
    VAR_NAME="var2"
    VALUES=(0.01 0.02 0.03 0.04 0.05)
    echo "Running Parameter Variation 2 (var2)..."
    
    for dt in "${VALUES[@]}"; do
        RUN_DIR="results/${VAR_NAME}/run_${dt}"
        mkdir -p "${RUN_DIR}"
        
        # Define simulation arguments
        NX=40
        NY=40
        DUR=5.0
        DT=$dt
        W_INT=0.1
        CA=1.0
        CB=3.0
        METHOD="implicit"
        REACT="on"
        
        # Create start_run.sh with mandatory explanatory comments in the root dir
        cat <<EOF > "start_run.sh"
#!/bin/bash
# This script runs a single simulation of the 2D Brusselator reaction-diffusion system.
# Positional Arguments:
#   1. nx             - Number of grid points in x-direction (e.g., $NX)
#   2. ny             - Number of grid points in y-direction (e.g., $NY)
#   3. duration       - Total simulation time in seconds (e.g., $DUR)
#   4. dt             - Time step size in seconds (e.g., $DT)
#   5. write_interval - Time interval between writing output files in seconds (e.g., $W_INT)
#   6. Ca             - Brusselator parameter Ca (e.g., $CA)
#   7. Cb             - Brusselator parameter Cb (e.g., $CB)
#   8. method         - Time integration method (rk4 or implicit) (e.g., $METHOD)
#   9. reactions      - Reaction toggle (on or off) (e.g., $REACT)

./bin/brusselator $NX $NY $DUR $DT $W_INT $CA $CB $METHOD $REACT
EOF
        chmod +x "start_run.sh"
        
        echo "  Running simulation for dt = ${dt}..."
        # Clean output directory before running
        rm -rf output && mkdir -p output
        
        # Run the start script
        ./start_run.sh > "simulation.log" 2>&1
        
        # Move generated files to the run directory
        mv start_run.sh "${RUN_DIR}/"
        mv simulation.log "${RUN_DIR}/"
        mv output/*.vts "${RUN_DIR}/"
        rm -rf output
    done
    echo "Parameter Variation 2 completed. Results saved in results/var2/"
}

# Custom Parameter Variation
run_custom() {
    PARAM=$1
    VALUES=($2)
    VAR_NAME=$3
    
    echo "Running Custom Parameter Variation '${VAR_NAME}'..."
    echo "Varying parameter '${PARAM}' over values: ${VALUES[*]}"
    
    # Defaults
    NX=40
    NY=40
    DUR=10.0
    DT=0.01
    W_INT=1.0
    CA=1.0
    CB=3.0
    METHOD="rk4"
    REACT="on"
    
    for val in "${VALUES[@]}"; do
        RUN_DIR="results/${VAR_NAME}/run_${val}"
        mkdir -p "${RUN_DIR}"
        
        # Apply custom parameter override
        if [ "$PARAM" = "nx" ]; then NX=$val
        elif [ "$PARAM" = "ny" ]; then NY=$val
        elif [ "$PARAM" = "duration" ]; then DUR=$val
        elif [ "$PARAM" = "dt" ]; then DT=$val
        elif [ "$PARAM" = "write_interval" ]; then W_INT=$val
        elif [ "$PARAM" = "Ca" ]; then CA=$val
        elif [ "$PARAM" = "Cb" ]; then CB=$val
        elif [ "$PARAM" = "method" ]; then METHOD=$val
        elif [ "$PARAM" = "reactions" ]; then REACT=$val
        else
            echo "Error: Unknown parameter '$PARAM'"
            exit 1
        fi
        
        # Create start_run.sh
        cat <<EOF > "start_run.sh"
#!/bin/bash
# This script runs a single simulation of the 2D Brusselator reaction-diffusion system.
# Positional Arguments:
#   1. nx             - Number of grid points in x-direction
#   2. ny             - Number of grid points in y-direction
#   3. duration       - Total simulation time in seconds
#   4. dt             - Time step size in seconds
#   5. write_interval - Time interval between writing output files in seconds
#   6. Ca             - Brusselator parameter Ca
#   7. Cb             - Brusselator parameter Cb
#   8. method         - Time integration method (rk4 or implicit)
#   9. reactions      - Reaction toggle (on or off)

./bin/brusselator $NX $NY $DUR $DT $W_INT $CA $CB $METHOD $REACT
EOF
        chmod +x "start_run.sh"
        
        echo "  Running simulation for ${PARAM} = ${val}..."
        rm -rf output && mkdir -p output
        
        ./start_run.sh > "simulation.log" 2>&1
        
        mv start_run.sh "${RUN_DIR}/"
        mv simulation.log "${RUN_DIR}/"
        mv output/*.vts "${RUN_DIR}/"
        rm -rf output
    done
    echo "Custom Parameter Variation completed. Results saved in results/${VAR_NAME}/"
}

# Parse Command Line Arguments
if [ "$#" -eq 1 ]; then
    if [ "$1" = "var1" ]; then
        run_var1
    elif [ "$1" = "var2" ]; then
        run_var2
    else
        echo "Error: Predefined mode must be 'var1' or 'var2'"
        echo "Usage: $0 [var1|var2] or $0 <parameter_name> \"<values>\" <variation_name>"
        exit 1
    fi
elif [ "$#" -eq 3 ]; then
    run_custom "$1" "$2" "$3"
else
    echo "Usage:"
    echo "  Predefined variations: $0 [var1|var2]"
    echo "  Custom variations:     $0 <parameter_name> \"<values>\" <variation_name>"
    exit 1
fi
