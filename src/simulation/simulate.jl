function simulate(parms::BDParameters, Nₑ::Float64; N_max::Union{Float64, Int64}=Inf, S_max::Union{Float64, Int64}=Inf, iter=1000)::Union{Outbreak, Nothing}
    proc = simulate(parms, N_max=N_max, S_max=S_max)
    if n_sampled(proc) < S_max
        count = 0
        while n_sampled(proc) < S_max && count <= iter
            count += 1
            proc = simulate(parms, N_max=N_max, S_max=S_max)
        end
        if n_sampled(proc) < S_max
            println("Sampled target not reached within $(iter) iterations")
            return
        end
    end
    tree = simulate(proc, Nₑ)
    return Outbreak(proc, PhyloTree(tree, Nₑ))
end