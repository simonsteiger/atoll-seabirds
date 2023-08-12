using CSV, DataFrames, Chain

pop = @chain begin
   CSV.read("data/atoll_seabird_populations_29Jul.csv", DataFrame)
   DataFrames.transform(_, All() .=> ByRow(x -> ismissing(x) ? 0 : x) => identity)
   subset(_, All() .=> ByRow(x -> x != -1))
end
