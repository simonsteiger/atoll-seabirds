# This script is part of the project associated with
# Article: Atolls are globally significant hubs for tropical seabirds
# Authors: Steibl S, Steiger S, Wegmann AS, Holmes ND, Young, HS, Carr P, Russell JC 
# Last edited: 2024-03-10

import Pkg

@info "Instantiating project"

Pkg.instantiate()

# If you are sampling from the posterior locally, you can adjust the number of samples and chains in each script.
# We recommend running at least 10 000 samples on four chains.

include(joinpath(dirname(Base.active_project()), "julia", "src", "constants.jl"))
include(joinpath(ROOT, "julia", "scripts", "pca.jl"))
include(joinpath(ROOT, "julia", "scripts", "presence.jl"))
include(joinpath(ROOT, "julia", "scripts", "count.jl"))
