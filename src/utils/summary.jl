
function summarize(parms::CRBDParameters)
    removal_rate = parms.μ + parms.r * parms.ψ
    println("Lifespan   : ", 1. / removal_rate)
    println("Growth rate: ", parms.λ - removal_rate)
    println("R₀         : ", parms.λ / removal_rate)
    print("Sample rate: ", parms.ψ / removal_rate)
end


# TODO: Update to incorporate mutation (currently not included)
function summarize(parms::MTBDParameters)
    removal_rate = parms.μ .+ parms.r .* parms.ψ
    println("Lifespan   : ", 1. ./ removal_rate)
    println("Growth rate: ", parms.λ .- removal_rate)
    Rᵢ = parms.λ ./ removal_rate
    println("Rᵢ         : ", Rᵢ)
    Λ = parms.λ[1, 1] - parms.λ[2,2] - removal_rate[1] + removal_rate[2]
    c = sqrt(Λ^2 + 4 * parms.λ[1,2] * parms.λ[2,1])
    f₁ = (c + Λ) / (c + Λ + 2. * parms.λ[1,2])
    println("R₀         : ", f₁ * (Rᵢ[1,1] + Rᵢ[1, 2]) + (1. - f₁) * (Rᵢ[2,1] + Rᵢ[2, 2]))
    println("Sample rate: ", parms.ψ ./ removal_rate)
    print("Frequencies :", [f₁, 1. - f₁])
end


"""
    n_sampled(linelist::DataFrame)

    Compute the number of sampled individuals in `linelist`.
"""
function n_sampled(linelist::DataFrame)
    sampled = 0
    @eachrow! linelist begin
        if :t_sam > 0.
            sampled += 1
        end
    end
    return sampled
end

@forward BDProcess.linelist n_sampled


"""
    n_deceased(linelist::DataFrame)

    Compute the number of deceased individuals in `linelist`.
"""
function n_deceased(linelist::DataFrame)
    deceased = 0
    @eachrow! linelist begin
        if :t_death > 0.
            deceased += 1
        end
    end
    return deceased
end


@forward BDProcess.linelist n_deceased


"""
    offspring_dist(linelist::DataFrame)

    Compute the offspring distribution for all individuals in `linelist`.

    Optional flag `deceased_only` determines whether or not to include active individuals in distribution.
"""
function offspring_dist(linelist::DataFrame; height=Inf, deceased_only=true)
    dist = fill(0, nrow(linelist))
    deceased = fill(false, nrow(linelist))
    @eachrow! linelist begin
        if :parent_id > 0
            dist[:parent_id] += 1
        end
        deceased[:child_id] = (:t_death < height || !deceased_only) ? true : false
    end
    return dist[deceased]
end


function offspring_dist(proc; deceased_only)
    return offspring_dist(proc.linelist, height=proc.height, deceased_only=deceased_only)
end


function offspring_dist(proc::BDProcess; deceased_only=true)
    @unpack linelist, height = proc
    dist = fill(0, nrow(linelist))
    deceased = fill(false, nrow(linelist))
    @eachrow! linelist begin
        if :parent_id > 0
            dist[:parent_id] += 1
        end
        deceased[:child_id] = (:t_death < height || !deceased_only) ? true : false
    end
    return dist[deceased]
end


"""
    type_dist(proc::BDProcess)

    Compute the distribution of different types in birth-process `proc`.
"""
function type_dist(proc::BDProcess)
    n_types = hasproperty(proc.parms, :n_types) ? proc.parms.n_types : 1
    dist = fill(0, n_types)
    @eachrow! proc.linelist begin
        dist[:child_type] += 1
    end
    return dist
end



function summarize(linelist::DataFrame, t::Float64)
    n = nrow(linelist)
    type = fill(-1, n)
    lifespan = zeros(n)
    offspring = fill(0, n)
    deceased = fill(false, n)
    sampled = fill(false, n)

    type[1] = linelist[1, :child_type]
    deceased[1] = !isinf(linelist[1, :t_death])
    sampled[1] = linelist[1, :t_sam] > 0.
    lifespan[1] = deceased[1] ? linelist[1, :t_death] - linelist[1, :t_birth] : t - linelist[1, :t_birth]

    @eachrow! linelist[2:end, :] begin
        type[:child_id] = :child_type
        deceased[:child_id] = !isinf(:t_death)
        sampled[:child_id] = :t_sam > 0.
        lifespan[:child_id] = deceased[:child_id] ? :t_death - :t_birth : t - :t_birth
        offspring[:parent_id] += 1
    end
    return type, deceased, sampled, lifespan, offspring
end



function summarize(proc::BDProcess)
    n_types = typeof(proc.parms) == MTBDParameters ? proc.parms.n_types : 1
    n = nrow(proc.linelist)
    end_time = proc.traj[1, end]

    sampled = zeros(n_types)
    deceased = zeros(n_types)
    frequency = zeros(n_types)


    # Track cumulative lifespans and offspring
    lifespan = zeros(n_types) .+ eps()
    offspring = zeros(n_types, n_types)
    mutations = zeros(n_types, n_types)

    @eachrow! proc.linelist[2:end, :] begin
        frequency[:child_type] += 1.
        if :event == 1
            offspring[:parent_type, :child_type] += 1
        elseif :event == 2
            mutations[:parent_type, :child_type] += 1
        end

        if !isinf(:t_death)
            deceased[:child_type] += 1
            lifespan[:child_type] += :t_death - :t_birth
        else
            lifespan[:child_type] += end_time - :t_birth
        end

        if :t_sam > 0.
            sampled[:child_type] += 1
        end
    end

    # Calculate expected lifespan and offspring
    δbar = deceased ./ lifespan .+ eps()
    λbar = offspring ./ lifespan
    ψbar = sampled ./ lifespan
    Ribar = λbar ./ δbar
    frequency ./= sum(frequency)
    Rbar = frequency' * Ribar * ones(n_types)

    return  ProcessSummary(δbar, λbar, ψbar, Ribar, Rbar, frequency, sampled)
end


@forward Outbreak.proc summarize