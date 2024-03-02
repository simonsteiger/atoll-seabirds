#!/bin/bash

# Check if the project directory path and additional arguments are provided
if [ $# -lt 3 ]; then
    echo "Usage: $0 <project_directory> <arg1> <arg2>"
    exit 1
fi

# Change directory to the provided directory path
cd "$1" || { echo "Error: Directory not found"; exit 1; }

# Define the paths to the scripts
julia_path="$1/julia/reproduce.jl"
renv_path="$1/renv/activate.R"
r_path="$1/R/steibl_et_al_2024_atoll_seabird_analysis_R-script.R"

# Execute the julia pipeline
julia --project "$julia_path" "$2" "$3"

# Activate renv
Rscript "$renv_path"

# Create figures
Rscript "$r_path"
