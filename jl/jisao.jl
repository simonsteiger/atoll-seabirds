module JISAO

using NetCDF
using Statistics

vars = [
    ["lat", "latitude"],
    ["lon", "longitude"],
    ["data", "precip_anom"]
]

dct = Dict()

# Open NC
[dct[v[2]] = NetCDF.open("data/preci_anom_cor.nc", v[1]) for v in vars]

# Why are some numbers in dct["precip_anom"][:, :, 2] so huge? Inf = 32767? Probably!

ifelse.(dct["precip_anom"][:, :, 2] .== 32767, missing, dct["precip_anom"][:, :, 2])

# Calculate absolute values of all cells
[dct[v[2]] = abs.(dct[v[2]]) for v in vars];

[dct[v[2]] = mean(dct[v[2]], dims=4) for v in vars];

end