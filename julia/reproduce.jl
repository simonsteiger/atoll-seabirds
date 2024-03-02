# NetCDF data frames are too large to share on GitHub
# We are only uploading the specifications we used to download them
# and the scripts we used to clean them.

import Pkg

Pkg.instantiate()

# If you are sampling from the posterior locally, you can adjust the number of samples and chains in each script.
# We recommend running at least 10 000 samples on four chains.

include(joinpath(dirname(Base.active_project()), "julia", "src", "constants.jl"))
include(joinpath(ROOT, "julia", "scripts", "pca.jl"))
include(joinpath(ROOT, "julia", "scripts", "presence.jl"))
include(joinpath(ROOT, "julia", "scripts", "count.jl"))
