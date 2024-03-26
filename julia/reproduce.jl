# This script is part of the project associated with
# Article: Atolls are globally significant sites for tropical seabirds
# Authors: Steibl S, Steiger S, Wegmann AS, Holmes ND, Young, HS, Carr P, Russell JC 
# Last edited: 2024-03-24

import Pkg

@info "Instantiating project"

Pkg.instantiate()

# If you are sampling from the posterior locally, you can adjust the number of samples and chains in each script.
# We recommend running at least 10 000 samples on four chains.

load = ARGS[1] == "false"
run_loocv = ARGS[2] == "true"
run_sensitivity = ARGS[3] == "true"

include(joinpath(dirname(Base.active_project()), "julia", "src", "constants.jl"));

active_presence_priorsettings = run_sensitivity ? PRESENCE_PRIORSETTINGS : ["default"]
active_count_priorsettings = run_sensitivity ? COUNT_PRIORSETTINGS : ["default"]
priorsetting = ""

include(joinpath(ROOT, "julia", "scripts", "pca.jl"));

for current_priorsetting in active_presence_priorsettings
    global priorsetting = current_priorsetting
    include(joinpath(ROOT, "julia", "scripts", "presence.jl"));
end

for current_priorsetting in active_count_priorsettings
    global priorsetting = current_priorsetting
    include(joinpath(ROOT, "julia", "scripts", "count.jl"));
    include(joinpath(ROOT, "julia", "scripts", "summaries.jl"));
end
