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

# Note on ARGS to julia script
# ----------------------------
# ARGS allows you to make slight modifications on how the scripts are run
# ARGS[1]: if `true`, load previously created chains
# This can save a lot of time if want to recreate summaries without refitting the models
# ARGS[2]: if `true`, run loo_cv
# Running loo_cv takes the longest out of all downstream analyses, and skipping it further speeds up recreating summaries
# => The fastest setup is `ARGS = [true, false]`, but it assumes that you have saved chains!
