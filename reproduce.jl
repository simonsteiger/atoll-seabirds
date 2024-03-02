# NetCDF data frames are too large to share on GitHub
# We are only uploading the specifications we used to download them
# and the scripts we used to clean them.

import Pkg

Pkg.instantiate()

const ROOT = dirname(Base.active_project())

# ARGS allows you to make slight modifications on how the scripts are run
# ARGS[1]: if `true`, load previously created chains
# This can save a lot of time if want to recreate summaries without refitting the models
# ARGS[2]: if `true`, run loo_cv
# Running loo_cv takes the longest out of all downstream analyses, and skipping it further speeds up recreating summaries
# => The fastest setup is `ARGS = [true, false]`, but it assumes that you have saved chains!
ARGS = [true, false]

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
