
abstract type BDParameters end


@with_kw mutable struct CRBDParameters <: BDParameters
    N₀::Int64 = 1   # initial population size
    λ::Float64      # birth rate
    μ::Float64      # death rate
    ψ::Float64      # extinct / ancestral sampling rate
    ρ₀::Float64     # extant sampling rate
    r::Float64      # removal probability (upon sampling)
    t_max::Float64  # maximum simulation time
end


@with_kw mutable struct MTBDParameters <: BDParameters
    n_types::Int64          # number of distinct subgroups / types
    N₀::Vector{Int64}       # initial population distribution of each type a
    λ::Matrix{Float64}      # birth rate of type a -> b
    μ::Vector{Float64}      # death rate
    γ::Matrix{Float64}      # mutation rate of type a -> b
    ψ::Vector{Float64}      # extinct / ancestral sampling rate
    ρ₀::Vector{Float64}     # extant sample rate
    r::Vector{Float64}      # removal probability (upon sampling)
    t_max::Float64  # maximum simulation time
end


function Base.convert(::Type{MTBDParameters}, parms::CRBDParameters)::MTBDParameters
    @unpack N₀, λ, μ, ψ, ρ₀, r, t_max = parms
    return MTBDParameters(n_types=1, N₀=[N₀], λ=diagm([λ]), μ=[μ], γ=diagm([0.]), ψ=[ψ], ρ₀=[ρ₀], r=[r], t_max=t_max)
end


function Base.copy(parms::T) where T <: BDParameters
    return T(parms)
end





mutable struct Host
	id::Int64
    type::Int64
    infector::Union{Int64, Nothing}
    infector_type::Union{Int64, Nothing}
	root::Union{Float64, Nothing}
	leaf_times::Vector{Float64}
	leaf_ids::Vector{Int64}
end






struct BDProcess
    parms::BDParameters
    traj::Matrix{Float64}
    linelist::DataFrame
end


struct ProcessSummary
    δbar::Vector{Float64}
    λbar::Matrix{Float64}
    ψbar::Vector{Float64}
    Ribar::Matrix{Float64}
    Rbar::Float64
    frequency::Vector{Float64}
    sampled::Vector{Int64}
end


function Base.show(io::IO, summ::ProcessSummary)
    println(io, "Sampled          : ", summ.sampled)
    println(io, "Removal rate     : ", round.(summ.δbar, digits=3))
    println(io, "Birth rate       : ", round.(summ.λbar, digits=3))
    println(io, "Sample rate      : ", round.(summ.ψbar, digits=3))
    println(io, "Sample prop.     : ", round.(summ.ψbar ./ summ.δbar, digits=3))
    println(io, "R (type)         : ", round.(summ.Ribar, digits=3))
    println(io, "R (total)        : ", round.(summ.Rbar, digits=3))
    print(io, "Frequency        : ", round.(summ.frequency, digits=3))
end


function Base.show(io::IO, proc::BDProcess)
    println(io, "Initial pop.     : ", proc.parms.N₀)
    println(io, "Final pop. (cum) : ", proc.traj[2:end, end], " (", nrow(proc.linelist), ")")
    print(io, "Sampled          : ", n_sampled(proc))
end


mutable struct PhyloTree
    Nₑ::Float64
    tree::DataFrame
    height::Float64
    origin::Float64
    tip_labels::Vector{String}
    nwk::Node
end


function PhyloTree(tree::DataFrame, Nₑ::Float64)::PhyloTree
    tree = binarize(tree)
    update_branch_lengths!(tree)
    tip_labels = get_tip_labels(tree)
    nwk = newick(tree)
    height = tree[1, :t] - tree[end-1, :t]
    origin = tree[1, :t] - tree[end, :t]
    return PhyloTree(Nₑ, tree, height, origin, tip_labels, nwk)
end


function Base.show(io::IO, tree::PhyloTree)
    println(io, "Tips   : ", tree.tip_labels)
    println(io, "Height : ", round(tree.height, digits=3))
    println(io, "Origin : ", round(tree.origin, digits=3))
end

mutable struct Outbreak
    proc::BDProcess
    tree::PhyloTree
    summ::ProcessSummary    
end


function Outbreak(proc::BDProcess, tree::PhyloTree)::Outbreak
    summ = summarize(proc)
    return Outbreak(proc, tree, summ)
end


function Base.show(io::IO, out::Outbreak)
    println(io, "Transmission tree ")
    println(io, out.proc)
    println(io, "")
    println(io, "Process summary ")
    println(io, out.summ)
    println(io, "")
    println(io, "Phylogenetic tree ")
    println(io, out.tree)
end
