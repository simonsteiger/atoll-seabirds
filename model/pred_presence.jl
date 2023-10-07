using Serialization, Turing

chain = deserialize("model/chains/2023-10-03_chain_nesting_noexclude.jls")

α = get_params(chain).α
β = get_params(chain).β
# Spatial parameters
etasq = get_params(chain).etasq
rhosq = get_params(chain).rhosq
γ = get_params(chain).γ
# PC 1
ω11 = get_params(chain).ω11
ω12 = get_params(chain).ω12
ω13 = get_params(chain).ω13
ω1 = [ω11, ω12, ω13]
# PC 2
ω21 = get_params(chain).ω21
ω22 = get_params(chain).ω22
ω23 = get_params(chain).ω23
ω2 = [ω21, ω22, ω23]
# PC 3
ω31 = get_params(chain).ω31
ω32 = get_params(chain).ω32
ω33 = get_params(chain).ω33
ω3 = [ω31, ω32, ω33]
# PC 4
ω41 = get_params(chain).ω41
ω42 = get_params(chain).ω42
ω43 = get_params(chain).ω43
ω4 = [ω41, ω42, ω43]
# PC 5
ω51 = get_params(chain).ω51
ω52 = get_params(chain).ω52
ω53 = get_params(chain).ω53
ω5 = [ω51, ω52, ω53]
# PC 6
ω61 = get_params(chain).ω61
ω62 = get_params(chain).ω62
ω63 = get_params(chain).ω63
ω6 = [ω61, ω62, ω63]

function prediction(atoll, species, nesting, α, β, γ, ω1, ω2, ω3, ω4, ω5, ω6)
    for i in eachindex(atoll)
        v = logistic(
            α[atoll[i]] +
            β[species[i]] +
            γ[atoll[i]] + # not sure if this enters the model like that
            ω1[nesting[i]][species[i]] * pc[i, 1] +
            ω2[nesting[i]][species[i]] * pc[i, 2] +
            ω3[nesting[i]][species[i]] * pc[i, 3] +
            ω4[nesting[i]][species[i]] * pc[i, 4] +
            ω5[nesting[i]][species[i]] * pc[i, 5] +
            ω6[nesting[i]][species[i]] * pc[i, 6]
        )
        presence[i] = Bernoulli(v)
    end

    return presence
end