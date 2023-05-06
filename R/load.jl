using NetCDF, Statistics, BenchmarkTools

var = ["nppv", "chl", "phyc"]

function readcp!(var)
  x = Dict()
  [x[v] = NetCDF.open(string("data/copernicus_", v, ".nc"), v) for v in var];
  [x[v] = abs.(x[v]) for v in var];
  [x[v] = mean(x[v], dims = 4) for v in var];
end

@benchmark $readcp!($var)