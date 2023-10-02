# split data function
# function split_data(df, species; p=0.50)
#     speciesdf = @chain df begin
#         subset(_, :species => ByRow(x -> x .== species))
#     end
#     #shuffled = shuffleobs(speciesdf)
#     train, test = splitobs(speciesdf; at=p, shuffle=true) .|> getobs
#     # Below code upsamples PC1, but we want to upsample 1-6 simultaneously. Concatenate vectors to array - how?
#     # return trainset, testset = oversample(speciesdf.PC1, speciesdf.presence))
# end

# Dicts for trainset, testset
# trainset = Dict()
# testset = Dict()

# Probably fails because some species can't be split without filter conditions
# [(trainset[s], testset[s]) = split_data(pop, s) for s in unique(pop.species)];

# for f in numerics, k in keys(trainset)
#     μ, σ = rescale!(trainset[k][!, f])
#     rescale!(testset[k][!, f], μ, σ)
# end

# Throw out atolls with mostly 0 counts

# Dicts for train, test, train_label, test_label
# train = Dict{String,Matrix}()
# test = Dict{String,Matrix}()
# train_label = Dict{String,Vector}()
# test_label = Dict{String,Vector}()

# # random intercept only
# @model function random_intercept(atoll, presence)
#     a = length(unique(atoll))
#     ᾱ ~ Normal()
#     τ ~ Exponential(1)
#     α ~ filldist(Normal(ᾱ, τ), a)
# 
#     for i in eachindex(presence)
#         v = logistic(α[atoll[i]])
#         presence[i] ~ Bernoulli(v)
#     end
# end;
# 
# m1 = random_intercept(all_atoll, all_presence)
# #chain1 = sample(m1, HMC(0.05, 10), 15_000, discard_initial=5000)
# 
# # species and atoll intercept
# @model function species_intercept(atoll, species, presence)
#     a = length(unique(atoll))
#     s = length(unique(species))
#     ᾱ ~ Normal()
#     τ ~ Exponential(1)
#     α ~ filldist(Normal(ᾱ, τ), a)
#     β ~ filldist(Normal(0, 1), s)
# 
#     for i in eachindex(presence)
#         v = logistic(α[atoll[i]] + β[species[i]])
#         presence[i] ~ Bernoulli(v)
#     end
# end;
# 
# m2 = species_intercept(all_atoll, all_species, all_presence)
# #chain2 = sample(m2, HMC(0.05, 10), 8000, discard_initial=2000)
# 
# @model function pc_species_intercept(atoll, species, pc, presence)
#     # Number of groups per predictor
#     a = length(unique(atoll))
#     s = length(unique(species))
# 
#     # Priors for atolls
#     ᾱ ~ Normal()
#     τ ~ Exponential(1)
#     α ~ filldist(Normal(ᾱ, τ), a)
#     # Priors for species effect
#     β ~ filldist(Normal(0, 1), s)
#     # Priors for PC
#     θ1 ~ filldist(Normal(0, 1), s)
#     θ2 ~ filldist(Normal(0, 1), s)
#     θ3 ~ filldist(Normal(0, 1), s)
#     θ4 ~ filldist(Normal(0, 1), s)
#     θ5 ~ filldist(Normal(0, 1), s)
#     θ6 ~ filldist(Normal(0, 1), s)
# 
#     for i in eachindex(presence)
#         v = logistic(
#             α[atoll[i]] +
#             β[species[i]] +
#             θ1[species[i]] * pc[i, 1] +
#             θ2[species[i]] * pc[i, 2] +
#             θ3[species[i]] * pc[i, 3] +
#             θ4[species[i]] * pc[i, 4] +
#             θ5[species[i]] * pc[i, 5] +
#             θ6[species[i]] * pc[i, 6]
#         )
#         presence[i] ~ Bernoulli(v)
#     end
# end;

# 
# 
# @model function interaction_model(atoll, species, pc, nesting, presence)
#     # Number of groups per predictor
#     a = length(unique(atoll))
#     s = length(unique(species))
#     n = length(unique(nesting))
# 
#     # Priors for atolls
#     ᾱ ~ Normal()
#     τ ~ Exponential(1)
#     α ~ filldist(Normal(ᾱ, τ), a)
#     # Priors for species effect
#     β ~ filldist(Normal(0, 1), s)
#     # Priors for species PC
#     θ1 ~ filldist(Normal(0, 1), s)
#     θ2 ~ filldist(Normal(0, 1), s)
#     θ3 ~ filldist(Normal(0, 1), s)
#     θ4 ~ filldist(Normal(0, 1), s)
#     θ5 ~ filldist(Normal(0, 1), s)
#     θ6 ~ filldist(Normal(0, 1), s)
#     # Priors for nesting PC
#     γ1 ~ filldist(Normal(0, 1), n)
#     γ2 ~ filldist(Normal(0, 1), n)
#     γ3 ~ filldist(Normal(0, 1), n)
#     γ4 ~ filldist(Normal(0, 1), n)
#     γ5 ~ filldist(Normal(0, 1), n)
#     γ6 ~ filldist(Normal(0, 1), n)
# 
# 
#     for i in eachindex(presence)
#         v = logistic(
#             α[atoll[i]] +
#             β[species[i]] +
#             θ1[species[i]] * pc[i, 1] +
#             θ2[species[i]] * pc[i, 2] +
#             θ3[species[i]] * pc[i, 3] +
#             θ4[species[i]] * pc[i, 4] +
#             θ5[species[i]] * pc[i, 5] +
#             θ6[species[i]] * pc[i, 6]
#         )
#         presence[i] ~ Bernoulli(v)
#     end
# end;
# 
# m3 = pc_species_intercept(all_atoll, all_species, all_pc, all_presence)
#chain3 = sample(m3, HMC(0.05, 10), 8000, discard_initial=2000)

# @model function pc_species_intercept(atoll, species, pc, presence)
#     # Number of groups per predictor
#     a = length(unique(atoll))
#     s = length(unique(species))
# 
#     # Priors for atolls
#     # ᾱ ~ Normal()
#     # τ ~ Exponential(1)
#     # α ~ filldist(Normal(ᾱ, τ), a)
#     # Priors for species effect
#     β̄ ~ Normal()
#     κ ~ Exponential(1)
#     β ~ filldist(Normal(β̄, κ), a, s)
#     # Priors for PC
#     θ1 ~ filldist(Normal(0, 1), s)
#     θ2 ~ filldist(Normal(0, 1), s)
#     θ3 ~ filldist(Normal(0, 1), s)
#     θ4 ~ filldist(Normal(0, 1), s)
#     θ5 ~ filldist(Normal(0, 1), s)
#     θ6 ~ filldist(Normal(0, 1), s)
# 
#     for i in eachindex(presence)
#         v = logistic(
#             # α[atoll[i]] +
#             β[atoll[i], species[i]] +
#             θ1[species[i]] * pc[i, 1] +
#             θ2[species[i]] * pc[i, 2] +
#             θ3[species[i]] * pc[i, 3] +
#             θ4[species[i]] * pc[i, 4] +
#             θ5[species[i]] * pc[i, 5] +
#             θ6[species[i]] * pc[i, 6]
#         )
#         presence[i] ~ Bernoulli(v)
#     end
# end;

# m4 = pc_species_intercept(all_atoll, all_species, all_pc, all_presence);
# chain4 = sample(m4, HMC(0.01, 10), 20_000, discard_initial=5000)

# using Serialization

# serialize("chain4.jls", chain4)

# Fit smaller sample
# sdf = @chain Preprocess.envscores begin
#     subset(_, :species => ByRow(x -> x in ("Anous_stolidus", "Phaethon_lepturus", "Onychoprion_fuscatus")), skipmissing=true)
#     # subset(_, :atoll => ByRow(x -> x in ("Tetiaroa", "Palmyra", "Cocos", "Midway", "Wake", "Laysan")), skipmissing=true)
# end
# 
# using StatsBase
# 
# sub_atoll = denserank(sdf.atoll) .|> Int64
# sub_species = denserank(sdf.species) .|> Int64
# sub_pc = Matrix(sdf[:, 1:6])
# sub_presence = sdf.presence .|> Float64
# 
# m4 = pc_species_intercept(sub_atoll, sub_species, sub_pc, sub_presence)
# chain4 = sample(m4, HMC(0.01, 10), 20_000, discard_initial=2000)
# 
# function printchain(chain)
#     x = describe(chain)[1]
#     show(DataFrame(x), allrows=true)
# end
# 
# printchain(chain4)

# df_anst = @chain CSV.read("data/envscores.csv", DataFrame) begin
#     subset(_, :species => x -> x .== "Anous_stolidus")
#     select(_, Not([:cond, :region, :Column1]))
# end