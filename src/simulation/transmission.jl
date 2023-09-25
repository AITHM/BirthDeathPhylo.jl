# Single-type constant rate birth death model


function simulate(N₀::Int64, 
                  λ::Float64, 
                  μ::Float64, 
                  ψ::Float64,
                  ρ₀::Float64,
                  r::Float64;
                  N_max::Union{Float64, Int64} = Inf,
                  S_max::Union{Int64, Float64} = Inf,
                  t_max::Union{Int64, Float64} = Inf)::BDProcess

    parms = CRBDParameters(N₀, λ, μ, ψ, ρ₀, r, t_max)

    return simulate(parms, N_max=N_max, S_max=S_max)
end


function simulate(parms::CRBDParameters; N_max::Union{Float64, Int64}=Inf, S_max::Union{Float64, Int64}=Inf)::BDProcess

    @unpack N₀, λ, μ, ψ, ρ₀, r, t_max = parms

    # Make sure some stopping criteria is specified
    N_max = isinf(min(N_max, S_max, t_max)) ? 1e3 : N_max + 0.

    linelist = DataFrame(event = fill(1, N₀),
                         child_id = collect(1:N₀),
                         child_type = fill(1, N₀),
                         parent_id = fill(0, N₀),
                         parent_type = fill(1, N₀),
                         t_birth = fill(0., N₀),
                         t_death = fill(Inf, N₀),
                         t_sam = fill(-1., N₀))

    # Convert population size to Float64
    N = N₀ + 0.

    # Initialize trajectory
    t = 0.
    t_out = [t]
    N_out = [N]

    # Pre-calculate event rates
    total_event_rate = λ + μ + ψ
    birth_threshold = λ / total_event_rate
    death_threshold = (λ + μ) / total_event_rate

    # Create a cohort of individuals
    active = collect(1:N₀)
    cum_inc = N₀

    # Track cumulative sampled
    S = 0.

    # Main simulation loop
    while (N > 0.) && (cum_inc < N_max) && (S < S_max)
        rnd = rand()
        t -= log(rnd) / (total_event_rate * N)
        t > t_max && break
        if rnd ≤ birth_threshold                                       # Birth event
            cum_inc += 1                                                # Increment cumulative incidence
            parent = sample(active)                                     # Sample random parent from current cohort
            push!(linelist, (1, cum_inc, 1, parent, 1, t, Inf, -1.))    # Add event to linelist
            push!(active, cum_inc)                                      # Update list of active individuals
            N += 1.                                                     # Update total population size

        elseif rnd ≤ death_threshold                                   # Death event
            deceased_idx = sample(1:length(active))                     # Sample index of deceased individual
            deceased = active[deceased_idx]                             # Retrieve id of deceased individual
            linelist[deceased, :t_death] = t                            # Set the time of death for deceased individual
            deleteat!(active, deceased_idx)                             # Update list of active individuals
            N -= 1.                                                     # Update total population size

        else                                                            # Sampling event
            sampled_idx = sample(1:length(active))                      # Sample index of sampled individual
            sampled = active[sampled_idx]                               # Retrieve id of sampled individual
            if linelist[sampled, :t_sam] < 0.                           # Increment cumulative sampled count (if not already sampled)
                S += 1.
            end
            linelist[sampled, :t_sam] = t                               # Update time of sampling for sampled individual
            if rand() < r                                               # Removal upon sampling
                linelist[sampled, :t_death] = t                         # Update time of death / removal time
                deleteat!(active, sampled_idx)                          # Update list of active individuals
                N -= 1.                                                 # Update total population size
            end
        end
        push!(t_out, t)                                                 # Update trajectory
        push!(N_out, N)
    end

    # Present day sampling
    if t ≥ t_max && N > 0.
        for i in active                                        
            linelist[i, :t_death] = t_max                               # Truncate survival time for remaining active individuals
            if rand() <  ρ₀                                             # Sample with probability ρ₀
                linelist[i, :t_sam] = t_max
            end
        end
    end

    # height = maximum(linelist.t_sam)

    return BDProcess(
                parms,
                vcat(t_out', N_out'),
                linelist#,
                # height
    )
end



# Multi-type birth death model


function simulate(N₀::Vector{Int64},                   # Initial population distribution (across types)
                  λ::Matrix{Float64},                  # Birth rate: a -> b
                  μ::Vector{Float64},                  # Death rate
                  γ::Matrix{Float64},                  # Mutation rate: a -> b
                  ψ::Vector{Float64},                  # Sampling rate
                  ρ₀::Vector{Float64},                 # Extant sampling probability
                  r::Vector{Float64};                  # Removal probability
                  N_max::Union{Float64, Int64} = Inf,  # Maximum population size
                  S_max::Union{Int64, Float64} = Inf,  # Maximum sample size
                  t_max::Union{Int64, Float64} = Inf)::BDProcess

    n_types = size(λ)[1]

    parms = MTBDParameters(n_types, N₀, λ, μ, γ, ψ, ρ₀, r, t_max)

    return simulate(parms, N_max=N_max, S_max=S_max)
end


function simulate(parms::MTBDParameters; N_max::Union{Float64, Int64}=Inf, S_max::Union{Float64, Int64}=Inf)::BDProcess

    @unpack n_types, N₀, λ, μ, γ, ψ, ρ₀, r, t_max = parms

    # Check all inputs are equal size / length
    @assert size(λ)[1] == size(λ)[2] == length(μ) == size(γ)[1] == size(γ)[2] == length(ψ) == length(ρ₀) == length(r)

    # Check all inputs are physical
    @assert all(λ .≥ 0.) && all(μ .≥ 0.) && all(γ .≥ 0.) && all(ψ .≥ 0.) && all(ρ₀ .≥ 0.) && all(ρ₀ .≤ 1.) &&  all(r .≥ 0.) && all(r .≤ 1.)

    # Determine the number of types
    n_types = length(μ)

    # Specify a maximum population size
    N_max = isinf(min(N_max, S_max, t_max)) ? 1e3 : N_max + 0.

    # Calculate total birth and mutation rates for each type
    Λ = sum(λ, dims=2)
    Γ = sum(γ, dims=2)

    # Throughout the simulation keep track of two key vectors
    # 1. [N]ᵢ - number of individuals of type i
    # 2. [active]ⱼ - ids of active individuals j ∈ 1, ..., ∑_(i=1)^(n_type) Nᵢ

    # Generate the initial distribution of N and active given inputs
    N_tot = 0
    active = [Vector{Int64}() for _ in 1:n_types]
    type = Vector{Int64}()
    for t in eachindex(N₀)
        for i in 1:N₀[t]
            N_tot += 1
            push!(active[t], N_tot)
            push!(type, t)
        end
    end

    # Initialize linelist of all events
    # - each birth and mutation event appears in a new row
    # - sampling and death update existing entries
    linelist = DataFrame(event = fill(1, N_tot),
                         child_id = collect(1:N_tot),
                         child_type = type,
                         parent_id = fill(0, N_tot),
                         parent_type = fill(0, N_tot),
                         t_birth = fill(0., N_tot),
                         t_death = fill(Inf, N_tot),
                         t_sam = fill(-1., N_tot))

    # Initialize the simulation
    N = reshape(N₀ .+ 0., (1, n_types))
    t = 0.
    t_out = [t]
    N_out = VectorOfArray([N[:]])

    # Record the cumulative event count
    row_num = N_tot

    # Track cumulative sampled
    S = 0.

    while (N_tot > 0) && (row_num < N_max) && (S < S_max)

        birth_rate = (N * Λ)[1]     # Cumulate birth rate (across all individuals and types)
        death_rate = (N * μ)[1]     # Cumulate death rate (across all individuals and types)
        mutation_rate = (N * Γ)[1]  # Cumulate mutation rate (across all individuals and types)
        sample_rate = (N * ψ)[1]    # Cumulate sample rate (across all individuals and types)

        # Compute the total event rate
        total_event_rate = birth_rate + death_rate + mutation_rate + sample_rate

        # Sample a random value
        rnd = rand()
        t -= log(rnd) / total_event_rate            # Update time
        t > t_max && break                          # Compare time with stopping criteria

        # Birth event
        if rnd <= birth_rate / total_event_rate
            row_num += 1
            parent_type = wsample(1:n_types, (N' .* Λ)[:])              # Sample type of parent
            parent_id = sample(active[parent_type])                     # Sample id of parent
            child_type = wsample(1:n_types, @view λ[parent_type, :])    # Sample type of new active individual
            push!(linelist, (1, row_num, child_type, 
                             parent_id, parent_type, t, Inf, -1.))      # Update linelist
            push!(active[child_type], row_num)                          # Update active
            N[child_type] += 1.                                         # Update population sizes
            N_tot += 1                                                  # Update total popuation size

        # Death event
        elseif rnd <= (birth_rate + death_rate) / total_event_rate
            deceased_type = wsample(1:n_types, (N' .* μ)[:])            # Sample type of deceased
            deceased_idx = sample(1:length(active[deceased_type]))      # Sample index of deceased individual
            deceased = active[deceased_type][deceased_idx]              # Determine id of deceased individual
            linelist[deceased, :t_death] = t                            # Update time of death in linelist
            deleteat!(active[deceased_type], deceased_idx)              # Update list of active individuals
            N[deceased_type] -= 1.                                      # Update population sizes
            N_tot -= 1

        # Mutation event
        elseif rnd <= (birth_rate + death_rate + mutation_rate) / total_event_rate
            row_num += 1                                        
            old_type = wsample(1:n_types, (N' .* Γ)[:])                 # Sample type of mutating individual
            mutant_idx = sample(1:length(active[old_type]))             # Sample index of mutating individual 
            mutant_id = active[old_type][mutant_idx]                    # Retrieve id of mutating individual
            new_type = wsample(1:n_types, @view γ[old_type, :])         # Sample new mutated type
            linelist[mutant_id, :t_death] = t                           # Record death of old individual from linelist
            push!(linelist, (0, row_num, new_type, 
                             mutant_id, old_type, t, Inf, -1.))         # Update linelist with new mutant
            deleteat!(active[old_type], mutant_idx)                     # Update list of active individuals
            push!(active[new_type], row_num)                            # Update list of active individuals
            N[old_type] -= 1.                                           # Update population sizes
            N[new_type] += 1.                                           

        # Sampling event
        else
            sampled_type = wsample(1:n_types, (N' .* ψ)[:])             # Sample type of sampled individual
            sampled_idx = sample(1:length(active[sampled_type]))        # Sample index of sampled individual
            sampled = active[sampled_type][sampled_idx]                 # Retrieve id of sampled individual
            if linelist[sampled, :t_sam] < 0.                           # Check if sampled individual has been sampled previously
                S += 1.                                                 # If not, increment cumulative sample count
            end
            linelist[sampled, :t_sam] = t                               # Update linelist with new sample time

            # Removal
            if rand() < r[sampled_type]                                 # Check if sampled individal is removed
                linelist[sampled, :t_death] = t                         # Update removal time of sampled individual
                deleteat!(active[sampled_type], sampled_idx)            # Update list of active individuals
                N[sampled_type] -= 1.                                   # Update population sizes
                N_tot -= 1                                              # Update total population size
            end
        end

        # Update outputs
        push!(t_out, t)
        push!(N_out, N[:])
    end

    # Present day sampling
    if t >= t_max && sum(N) > 0.                 # If any active individuals remain
        for t in 1:n_types
            for i in active[t]
                linelist[i, :t_death] = t_max    # Truncate their time of death
                if rand() <  ρ₀[t]                  # If they are sampled
                    linelist[i, :t_sam] = t_max  # Update their sampling time
                end
            end
        end
    end


    return BDProcess(
                parms,
                vcat(t_out', convert(Array, N_out)),
                linelist
    )
end
