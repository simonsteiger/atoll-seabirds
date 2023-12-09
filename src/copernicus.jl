module Copernicus

using NetCDF, CSV
using Statistics, DataFrames, Compat, Chain, StatsFuns

const KEYS = [
  "nppv", 
  "chl",
  "phyc",
  "temp_" .* string.(1:2)...,
  "velo_"  .* string.(1:4)...,
  "longitude",
  "latitude",
  [x .* y for x in ["longitude", "latitude"], y in ["_velo", "_temp"]]...,
]

const FILENAMES = [
  "nppv",
  "chl",
  "phyc",
  "temp_" .* string.(1:2)...,
  "velo_"  .* string.(1:4)...,
  "nppv", "nppv",
  "velo_1", "velo_1",
  "temp_1", "temp_1",
]

const VARS = [
  "nppv",
  "chl",
  "phyc",
  "thetao",
  "thetao",
  "vo",
  "vo",
  "vo",
  "vo",
  "longitude",
  "latitude",
  "longitude",
  "latitude",
  "longitude",
  "latitude",
]

# Create a dictionary to write to
dct = Dict{String,Any}(Pair.(KEYS, NetCDF.open.("data/copernicus_" .* FILENAMES .* ".nc", VARS)))

# Calculate absolute values of all cells
[dct[k] = abs.(dct[k]) for k in KEYS[1:9]]; # only nppv, chl, phyc

# Calculate mean of array slices
[dct[k] = mean(dct[k], dims=4) for k in KEYS[1:9]]

# Temperature: Multiply by scale factor and add offset
[dct[k] = dct[k] .* 0.0007324442267417908 .+ 21 for k in ["temp_1", "temp_2"]]

# Velocity: Multiply by scale factor
[dct[k] = dct[k] .* 0.0006103701889514923 for k in ["velo_1", "velo_2", "velo_3", "velo_4"]]


for k in keys(dct)
  dct[k] = allowmissing(Float64.(dct[k]))
end

[map!(x -> isinf(x) ? missing : x, dct[k], dct[k]) for k in KEYS[1:9]];

# Concatenate temp arrays on time axis
dct["temp"] = @chain Base.cat(dct["temp_1"], dct["temp_2"], dims=3) begin
  mean(_, dims=3)
  map(x -> x == 45.0 ? missing : x, _) # 45 is the new missing after rescale
end

# Concatenate velo arrays on time axis
dct["velo"] = @chain Base.cat(dct["velo_1"], dct["velo_2"], dct["velo_3"], dct["velo_4"], dims=3) begin
  mean(_, dims=3)
  map(x -> x == 20.0 ? missing : x, _) # 20 is the new missing after rescale
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

# Join all dfs stored in dict to single df
cp_df = @chain begin
  get.(Ref(dctdf), ["nppv", "chl"], missing)
  reduce((x, y) -> leftjoin(x, y, on=[:latitude, :longitude]), _)
end

# problem here might be (is!) that the longitude and latitude data are set to abs

# Calculate outcome means for ± 1 long/lat around each atoll
function extract(df, env_lat, env_long; tol=1)
  cp_ll = [:latitude, :longitude]

  out = @chain df begin
    subset(_, cp_ll => ByRow((x, y) -> ≈(x, env_lat; atol=tol) && ≈(y, env_long; atol=tol)))
    transform(_, Not(cp_ll) .=> ByRow(x -> isinf(x) ? missing : x) .=> identity) # Inf represents land, recode to missing
    combine(_, Not(cp_ll) .=> (x -> mean(skipmissing(x))) .=> identity)
  end

  return out
end

envs = @chain begin
  CSV.read("data/seabird_atolls_envs_02May.csv", DataFrame)
  transform(:long => (x -> ifelse.(x .< 0, x .+ 360, x)) => identity)
end

cp_summary = reduce(vcat, [extract(cp_df, i...) for i in eachrow(envs[!, [:lat, :long]])])
cp_summary.atoll = envs.atoll

# Join by atoll

end
