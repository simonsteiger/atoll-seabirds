using NetCDF # read
using Statistics, DataFrames, Compat, Chain # wrangling

# First element is key in dictionary
# Second entry is csv file name
# Third entry is variable name in nc file

vars = [
  ["nppv", "nppv", "nppv"],
  ["chl", "chl", "chl"],
  ["phyc", "phyc", "phyc"],
  ["temp_1", "temp_1", "thetao"],
  ["temp_2", "temp_2", "thetao"],
  ["velo_1", "velo_1", "vo"],
  ["velo_2", "velo_2", "vo"],
  ["velo_3", "velo_3", "vo"],
  ["velo_4", "velo_4", "vo"],
  ["longitude", "nppv", "longitude"],
  ["latitude", "nppv", "latitude"],
  ["longitude_velo", "velo_1", "longitude"],
  ["latitude_velo", "velo_1", "latitude"],
  ["longitude_temp", "temp_1", "longitude"],
  ["latitude_temp", "temp_1", "latitude"],
]

# Create a dictionary to write to
dct = Dict()

# Read nc files
[dct[v[1]] = NetCDF.open(string("data/copernicus_", v[2], ".nc"), v[3]) for v in vars];

# Calculate absolute values of all cells
[dct[v[1]] = abs.(dct[v[1]]) for v in vars];

# Calculate mean of array slices containing vars[2]
[dct[v[1]] = mean(dct[v[1]], dims=4) for v in vars];

# Concatenate temp arrays on time axis
dct["temp"] = @chain Base.cat(dct["temp_1"], dct["temp_2"], dims=3) begin
  mean(_, dims=3)
end

# Concatenate velo arrays on time axis
dct["velo"] = @chain Base.cat(dct["velo_1"], dct["velo_2"], dct["velo_3"], dct["velo_4"], dims=3) begin
  mean(_, dims=3)
end

# Remove no longer necessary arrays from dictionary
[delete!(dct, k) for k in ["temp_1", "temp_2", "velo_1", "velo_2", "velo_3", "velo_4"]]

# Final names for matrices
nms = ["nppv", "chl", "phyc", "temp", "velo"]

# Drop singleton dimensions of arrays
[dct[n] = dropdims(dct[n], dims=d) for n in nms, d in [4, 3]]

dctdf = Dict{String,DataFrame}()

function expandlat(lat, long)
  out = @chain collect(Base.Iterators.repeated(lat, length(long))) begin
    collect(Base.Iterators.flatten(_))
  end

  return out
end

function writecp(dct, var, long, lat)
  iterlat = expandlat(lat, long)

  out = @chain DataFrame(dct[var], :auto) begin
    insertcols(_, 1, :longitude => long)
    stack(_, Not(:longitude))
    insertcols(_, 2, :latitude => iterlat)
    select(_, Not(:variable))
    rename(_, :value => var)
  end

  return out
end

[dctdf[n] = writecp(dct, n, dct["longitude"], dct["latitude"]) for n in ["nppv", "chl", "phyc"]]

[dctdf[n] = writecp(dct, n, dct[string("longitude_", n)], dct[string("latitude_", n)]) for n in ["velo", "temp"]]

cp_data =
  let v = collect(dctdf)
    @chain [v[i][2] for i in eachindex(v)] begin
      reduce((x, y) -> innerjoin(x, y, on=[:latitude, :longitude]), _)
    end
  end

function filter_copernicus(data, lat, long)
  out = @chain data begin
    subset(_, (:latitude, :longitude) => ByRow((x, y) -> x >= lat - 1 & x <= lat + 1 & y >= long - 1 & y <= long + 1))
    # also mean per variable
  end

  return out
end

