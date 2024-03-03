#!/bin/bash

# Check if the project directory path and additional arguments are provided
if [ $# -lt 4 ]; then
    echo "Usage: $0 <project_directory> <sample?> <loocv?> <sensitivity?>"
    exit 1
fi

# Change directory to the provided directory path
cd "$1" || { echo "Error: Directory not found"; exit 1; }

# Define the paths to the scripts
julia_path="$1/julia/reproduce.jl"
rprofile_path="$1/.Rprofile"
r_path="$1/R/steibl_et_al_2024_atoll_seabird_analysis_R-script.R"

# Execute the julia pipeline
julia --project --threads 4 "$julia_path" "$2" "$3" "$4"

# Activate Rproject
Rscript "$rprofile_path"

# Create figures
Rscript "$r_path"
