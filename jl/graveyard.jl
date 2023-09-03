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