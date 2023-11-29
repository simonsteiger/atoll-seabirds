module ParamPlots

export plotparams

using Turing, Chain, StatsPlots

function plotparams(chain::Chains, param::String, n::Int64; lab=nothing)

    # add regex to automatically figure out the n parameters matching param

    p = fill(plot(), n)

    for i in 1:n
        
        sm = @chain DataFrame(group(chain, "$param$i")) begin
            stack(_)
            groupby(_, :variable)
            combine(_, :value => (x -> (mean=mean(logistic.(x)), std=std(logistic.(x)))) => AsTable)
        end

        p[i] = scatter(sm.mean, eachindex(sm.mean), xerror=sm.std, title="$param$i", color=i, label=false)
        
        if !isnothing(lab)
            yt = i in (1, 4) ? lab : fill(" ", length(lab))
            yticks!([0.5:1:length(lab)-0.5;], yt) 
            yflip!()
        end

    end

    out = plot([i for i in p]..., size=(1000, 500), layout=(1,n))

    return out
end

#nestingtypes = String.(unique(Preprocess.envs_known.nestingtype))
#species_names = replace.(sort(unique(Preprocess.envs_known.species)), r"[a-z]+_" => " ")
#
#ps = [begin plotparams(chain, "θ$i", 3, lab=nothing); xlims!(0, 1) end for i in 1:6]
#
#df_chain = @chain DataFrame(chain) begin
#    select(_, r"ω̄")
#    stack(_)
#    groupby(_, :variable)
#    combine(_, :value => (x -> (std=std(logistic.(x)), mean=mean(logistic.(x)))) .=> AsTable)
#end
#
#groups = chop.(df_chain.variable, tail=1)
#scatter(df_chain.mean, df_chain.variable, xerror=df_chain.std, group=groups)
#yticks!((1:18).-0.5, df_chain.variable)
#xlims!(0, 1)
#yflip!()

end
