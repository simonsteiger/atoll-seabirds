# NetCDF data frames are too large to share on GitHub
# We are only uploading the specifications we used to download them
# and the scripts we used to clean them.

import Pkg

Pkg.instantiate()

const ROOT = dirname(Base.active_project())

# Set ARGS to true if you have saved chains from a previous run which you want to load
ARGS = false

include("$ROOT/scripts/pca.jl")
include("$ROOT/scripts/presence.jl")
include("$ROOT/scripts/count.jl")

# TODO decide if Rmd file should be included here
# using RCall
# 
# R"""
# library(here)
# source(here("R", "ggtestplots.Rmd"))
# """
