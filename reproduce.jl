# NetCDF data frames are too large to share on GitHub
# We are only uploading the specifications we used to download them
# and the scripts we used to clean them.

const ROOT = dirname(Base.active_project())

include("$ROOT/scripts/pca.jl")
include("$ROOT/scripts/presence.jl")
include("$ROOT/scripts/count.jl")
include("$ROOT/scripts/globalestimates.jl")
include("$ROOT/scripts/plotprep.jl")
