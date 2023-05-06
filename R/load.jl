using NetCDF # read
using Statistics, DataFrames, Compat, Pipe # wrangling

# First element is file path
# Second entry is loaded variable
vars = [
  ["nppv", "nppv"],
  ["chl", "chl"],
  ["phyc", "phyc"],
  #=
  name of sea surface temperature variable?
  ["temp_1", "???"],
  ["temp_2", "???"],
  =#
  ["nppv", "longitude"],
  ["nppv", "latitude"],
]

# Create a dictionary to write to
d = Dict()

# Read nc files
[d[v[2]] = NetCDF.open(string("data/copernicus_", v[1], ".nc"), v[2]) for v in vars];

# Calculate absolute values of all cells
[d[v[2]] = abs.(d[v[2]]) for v in vars];

# Calculate mean of array slices containing vars[2]
[d[v[2]] = mean(d[v[2]], dims=4) for v in vars];

# Drop singleton dimensions of arrays
[d[v[2]] = dropdims(d[v[2]], dims=n) for v in vars[1:3], n in [4, 3]]

iterlat = @pipe collect(Base.Iterators.repeated(d["latitude"], length(d["longitude"]))) |>
                collect(Base.Iterators.flatten(_))

test = @pipe DataFrame(d["nppv"], :auto) |>
             insertcols(_, 1, :longitude => d["longitude"]) |>
             stack(_, Not(:longitude)) |>
             insertcols(_, 2, :latitude => iterlat)

# Pipeline for higher-resolution velo data will be different
# DataFrame has different dims!
