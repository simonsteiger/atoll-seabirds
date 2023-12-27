using RCall

R"""
library(here)
renv::restore()
source(here("R/write.R"))
"""

const ROOT = dirname(Base.active_project())

include("$ROOT/scripts/pca.jl")
include("$ROOT/scripts/presence.jl")
include("$ROOT/scripts/count.jl")
include("$ROOT/scripts/globalestimates.jl")
include("$ROOT/scripts/plotprep.jl")

# ggtestplots.R to be added
