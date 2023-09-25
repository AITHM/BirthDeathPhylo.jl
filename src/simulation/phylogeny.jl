# Benchmarked !!
function constrain_coalescences!(const_lower::Vector{Float64},
                                 const_upper::Vector{Float64},
                                 const_lineages::Vector{Int64},
                                 const_events::Vector{Int64},
                                 lineages::Vector{Int64},
                                 times::Vector{Float64},
                                 leaves::Vector{Int64},
                                 Nₑ::Float64,
                                 bound::Float64,
                                 norm_tol::Float64)::Float64
    # Total the leaves
    total_leaves = sum(leaves)
    # Likelihood increment
    likelihood = 1.0
    # Initial constraints, bound interval
    events = lineages[2] - lineages[1]
    
    # Counter
    c = 1
    for m in 1:events
        const_lower[c] = bound
        const_upper[c] = times[1]
        const_lineages[c] = lineages[2]
        const_events[c] = events
        c += 1
    end

    for k in 2:length(times)
        events = leaves[k-1] + lineages[k+1] - lineages[k]

        for m in 1:events
            const_lower[c] = times[k-1]
            const_upper[c] = times[k]
            const_lineages[c] = lineages[k+1]
            const_events[c] = events
            c += 1
        end
    end

    # Separate coalescence events
    for c in 1:(total_leaves - 1)
        while (const_events[c] > 1)
            events = const_events[c]
            dt = 0.5 * (const_upper[c] - const_lower[c])
            mid_time = const_upper[c] - dt
            prob_norm = homochronous_probability(const_lineages[c], const_lineages[c] - events, 2.0 * dt, Nₑ)
            sig_loss = significance_loss(const_lineages[c], const_lineages[c] - events, dt, Nₑ)
            if sig_loss > norm_tol
                u = rand()
                events_lhs = 0
                events_rhs = events - events_lhs
                prob_rhs = homochronous_probability(const_lineages[c], const_lineages[c] - events_rhs, dt, Nₑ)
                prob_lhs = homochronous_probability(const_lineages[c] - events_rhs, const_lineages[c] - events, dt, Nₑ)
                sum_prob = (prob_lhs * prob_rhs) / prob_norm
                while u > sum_prob
                    events_lhs += 1
                    events_rhs -= 1
                    prob_rhs = homochronous_probability(const_lineages[c], const_lineages[c] - events_rhs, dt, Nₑ)
                    prob_lhs = homochronous_probability(const_lineages[c] - events_rhs, const_lineages[c] - events, dt, Nₑ)
                    sum_prob += (prob_lhs * prob_rhs) / prob_norm
                end
                likelihood *= (prob_lhs * prob_rhs) / prob_norm
            else
                events_lhs = trunc(Int, events / 2)
                events_rhs = events - events_lhs
                likelihood = 0.0
            end

            for m in 0:(events - 1)
                if m < events_lhs
                    const_upper[c + m] = mid_time
                    const_lineages[c + m] -= events_rhs
                    const_events[c + m] = events_lhs
                else
                    const_lower[c + m] = mid_time
                    const_events[c + m] = events_rhs
                end
            end
        end
    end
    return likelihood
end


constrain_coalescences!(const_lower::Vector{Float64},
                        const_upper::Vector{Float64},
                        const_lineages::Vector{Int64},
                        const_event::Vector{Int64},
                        lineages::Vector{Int64},
                        times::Vector{Float64},
                        leaves::Vector{Int64},
                        Nₑ::Float64,
                        bound::Float64)::Float64 = constrain_coalescences!(const_lower::Vector{Float64},
                                                                           const_upper::Vector{Float64},
                                                                           const_lineages::Vector{Int64},
                                                                           const_event::Vector{Int64},
                                                                           lineages::Vector{Int64},
                                                                           times::Vector{Float64},
                                                                           leaves::Vector{Int64},
                                                                           Nₑ::Float64,
                                                                           bound::Float64,
                                                                           1.0e-10)::Float64


# Benchmarked !!!
function homochronous_probability(n_start::Int, 
                                  n_end::Int, 
                                  dt::Float64, 
                                  Nₑ::Float64)::Float64

    # Check that inputs are valid
    (n_start <= 0 || n_end <= 0 || n_start < n_end || dt < 0 || Nₑ <= 0) && return 0.0
    (n_start == 1 && n_end == 1) && return 1.0

    if n_end == 1
        # Initialise total probability
        prob_total = 0.0
        for k in 2:n_start
            dk = convert(Float64, k)
            # Initialise probability increment
            prob_increment = 1.0

            for l in 2:n_start
                dl = convert(Float64, l)
                # Calculate coefficients
                if l != k 
                    prob_increment *= (dl * (dl - 1.0)) /
                                                (dl * (dl - 1.0) - dk * (dk - 1.0))
                end
            end

            # Calculate exponential
            prob_increment *= 1.0 - exp(- ((dk * (dk - 1.0)) / (2.0 * Nₑ)) * dt)

            # Update total probability
            prob_total += prob_increment
        end
    else
        # Initialise total probability
        prob_total = 0.0
        for k in n_end:n_start
            dk = convert(Float64, k)

            # Initialise probability increment
            prob_increment = (dk * (dk - 1.0)) / (convert(Float64, n_end) * (convert(Float64, n_end) - 1.0))
            
            for l in n_end:n_start
                dl = convert(Float64, l)
                # Calculate coefficients
                if l != k
                    prob_increment *= (dl * (dl - 1.0)) /
                                                (dl * (dl - 1.0) - dk * (dk - 1.0))
                end
            end

            # Calculate exponential
            prob_increment *= exp(- ((dk * (dk - 1.0)) / (2.0 * Nₑ)) * dt)

            # Update total probability
            prob_total += prob_increment
        end
    end
    return prob_total
end



# Benchmarked !!!
# Forward algorithm for the bounded coalescent
function forward_algorithm(times::Vector{Float64},
                           leaves::Vector{Int64},
                           Nₑ::Float64,
                           bound::Float64)::Matrix{Float64}

    n_times = length(times)
    sum_leaves = leaves[end]

    # Initiate forward algorithm
    forward_probs = zeros(sum(leaves), n_times+1)
    forward_probs[leaves[end], end] = 1.0

    # Forward recursion through sampling times
    for k in n_times:-1:2
        dt = times[k] - times[k - 1]
        for n_end in 1:sum_leaves
            for n_start in 1:sum_leaves
                transition_prob = homochronous_probability(n_start, n_end, dt, Nₑ)
                forward_probs[n_end + leaves[k-1], k] += transition_prob * forward_probs[n_start, k+1]
            end
        end
        sum_leaves += leaves[k-1]
    end

    # Bound probabilities
    dt = times[1] - bound

    for n_end in 1:sum_leaves
        for n_start in 1:sum_leaves
            transition_prob = homochronous_probability(n_start, n_end, dt, Nₑ)
            forward_probs[n_end, 1] += transition_prob * forward_probs[n_start, 2]
        end
    end
    return forward_probs
end



# Benchmarked !!
# Backward sampler for the bounded coalescent
function backward_sampler!(lineages::Vector{Int64},
                           forward_probs::Matrix{Float64},
                           times::Vector{Float64},
                           leaves::Vector{Int64},
                           Nₑ::Float64,
                           bound::Float64,
                           bound_size::Int64)::Float64
    # Initialise with bound condition
    lineages[1] = bound_size

    # Likelihood of sample
    likelihood = 1.0

    # Total the leaves
    total_leaves = sum(leaves)

    # Initiate backwards sampling
    dt = times[1] - bound

    # Generate random number for sampling
    u = rand()

    # Track cumulative probability for sampling
    sum_prob = 0.0

    # Evaluate smoothed probabilities, sample, and update likelihood
    for n_start in 1:total_leaves
        transition_prob = homochronous_probability(n_start, lineages[1], dt, Nₑ)
        smoothed_prob = (transition_prob * forward_probs[n_start, 2]) /
                            forward_probs[lineages[1], 1]

        sum_prob += smoothed_prob

        if (u < sum_prob)
            lineages[2] = n_start
            likelihood *= smoothed_prob
            break
        end
    end

    # Backwards sampling recursion
    for k in 2:length(times)
        dt = times[k] - times[k - 1]
        # Generate random number for sampling
        u = rand()

        # Track cumulative probability for sampling
        sum_prob = 0.0

        # Evaluate smoothed probabilities, sample, and update likelihood
        for n_start in 1:total_leaves
            # Number of lineages before leaves are added
            n_end = lineages[k] - leaves[k - 1]
            transition_prob = homochronous_probability(n_start, n_end, dt, Nₑ)
            smoothed_prob = (transition_prob * forward_probs[n_start, k + 1]) /
                                forward_probs[lineages[k], k]

            sum_prob += smoothed_prob
            if (u < sum_prob)
                lineages[k + 1] = n_start
                likelihood *= smoothed_prob
                break
            end
        end
    end
    # Return likelihood
    return likelihood
end


backward_sampler!(lineages::Vector{Int64},
                  forward_probs::Matrix{Float64},
                  times::Vector{Float64},
                  leaves::Vector{Int64},
                  Nₑ::Float64,
                  bound::Float64)::Float64 = backward_sampler!(lineages::Vector{Int64},
                                                               forward_probs::Matrix{Float64},
                                                               times::Vector{Float64},
                                                               leaves::Vector{Int64},
                                                               Nₑ::Float64,
                                                               bound::Float64,
                                                               1)::Float64





# Benchmarked !!
# Determine loss of significant figures in probability calculations
function significance_loss(n_start::Int64, 
                           n_end::Int64, 
                           dt::Float64, 
                           Nₑ::Float64)::Float64
    max_increment = 0.0

    # Check that inputs are valid
    (n_start <= 0 || n_end <= 0 || n_start < n_end || dt < 0 || Nₑ <= 0) && return 1.0
    (n_start == 1 && n_end == 1) && return 1.0

    if n_end == 1
        # Initialise total probability
        prob_total = 0.0
        for k in 2:n_start
            dk = convert(Float64, k)
            # Initialise probability increment
            prob_increment = 1.0
            for l in 2:n_start
                dl = convert(Float64, l)
                if l != k
                    prob_increment *= (dl * (dl - 1.0)) /
                                                (dl * (dl - 1.0) - dk * (dk - 1.0))
                end
            end

            # Calculate exponential
            prob_increment *= 1.0 - exp(- ((dk * (dk - 1.0)) / (2.0 * Nₑ)) * dt)

            if max_increment < abs(prob_increment)
                max_increment = abs(prob_increment)
            end

            # Update total probability
            prob_total += prob_increment
        end
    else
        # Initialise total probability
        prob_total = 0.0
        for k in n_end:n_start
            dk = convert(Float64, k)
            # Initialise probability increment
            prob_increment = (dk * (dk - 1.0)) / (convert(Float64, n_end) * (convert(Float64, n_end) - 1.0))
            for l in n_end:n_start
                dl = convert(Float64, l)
                # Calculate coefficients
                if l != k
                    prob_increment *= (dl * (dl - 1.0)) /
                                                (dl * (dl - 1.0) - dk * (dk - 1.0))
                end
            end

            # Calculate exponential
            prob_increment *= exp(- ((dk * (dk - 1.0)) / (2.0 * Nₑ)) * dt)

            if max_increment < abs(prob_increment)
                max_increment = abs(prob_increment)
            end

            # Update total probability
            prob_total += prob_increment
        end
    end
    return prob_total / max_increment
end


# Benchmarked !!
# Double-checked
function sample_bounded_times(times::Vector{Float64},
                              leaves::Vector{Int64},
                              Nₑ::Float64,
                              bound::Float64,
                              n_sam::Int64)::Tuple{Matrix{Float64}, Vector{Float64}}

    # Total the leaves
    total_leaves = sum(leaves)

    # Storage for samples
    sampled_times = Matrix{Float64}(undef, n_sam, (total_leaves - 1))

    # Storage for likelihood
    likelihood = ones(n_sam)

    forward_probs = forward_algorithm(times, leaves, Nₑ, bound)

    # Random numbers for sampling
    u = rand(n_sam, (total_leaves - 1))

    lineages = fill(1,length(times) + 1)
    const_lower = Vector{Float64}(undef, total_leaves - 1)
    const_upper = Vector{Float64}(undef, total_leaves - 1)
    const_lineages = Vector{Int64}(undef, total_leaves - 1)
    const_events = Vector{Int64}(undef, total_leaves - 1)

    for n in 1:n_sam
        likelihood[n] *= backward_sampler!(lineages, forward_probs, times, leaves, Nₑ, bound)

        likelihood[n] *= constrain_coalescences!(const_lower, const_upper, const_lineages, const_events,
                                                    lineages, times, leaves, Nₑ, bound)

        for c in 1:(total_leaves - 1)
            dl = convert(Float64, const_lineages[c])
            z = (Nₑ / (dl - 1.0)) * (1.0 - exp(((dl - 1.0) / Nₑ) *
                    (const_lower[c] - const_upper[c])))

            sampled_times[n, c] = const_upper[c] + (Nₑ / (dl - 1.0)) *
                                    log(1.0 - ((dl - 1.0) / Nₑ) * z * u[n, c])

            likelihood[n] *= (1.0 / z) * exp(((dl - 1.0) / Nₑ) *
                                (sampled_times[n, c] - const_upper[c]))
        end
    end
  return sampled_times, likelihood
end



function sample_bounded_times(times::Vector{Float64},
                              leaves::Vector{Int64},
                              Nₑ::Float64,
                              bound::Float64)::Tuple{Vector{Float64}, Float64}
    sampled_times, likelihood = sample_bounded_times(times, leaves, Nₑ, bound, 1)
    return vec(sampled_times), likelihood[1]
end



function sample_topology(leaf_times::Vector{Float64},
                         leaves::Vector{Int64},
                         coalescence_times::Vector{Float64})::Tuple{Matrix{Int64}, Vector{Float64}, Float64, DataFrame}

    # Total the leaves
    total_leaves = sum(leaves)

    # Internal node ancestors
    edge = Matrix{Int64}(undef, 2*(total_leaves - 1), 2)
    edge_length = Vector{Float64}(undef, 2*(total_leaves - 1))
    # Node ancestors and times
    node_times = fill(0., 2 * total_leaves - 1)

    # Counters
    l_index = total_leaves
    n_index = 2 * total_leaves - 1

    # Current information
    active_nodes = fill(0, 2*total_leaves - 1)
    total_active_nodes = 0

    # Node dataframe
    node_df = DataFrame(t = reduce(vcat, fill.(leaf_times, leaves)),
                        id = 1:total_leaves,
                        left = fill(0, total_leaves),
                        right = fill(0, total_leaves))

    # Random numbers for sampling
    u = rand(2*(total_leaves - 1))

    # Likelihood
    likelihood = 1.0

    k = length(leaf_times)

    for i in 1:leaves[k]
        active_nodes[l_index] = 1
        total_active_nodes += 1
        node_times[l_index] = leaf_times[k]
        l_index -= 1
    end

    k -= 1

    c = length(coalescence_times)

    while (c >= 1)

        if (k < 1 || leaf_times[k] < coalescence_times[c])
            anc_1 = 1
            sum_prob = convert(Float64, active_nodes[anc_1]) / convert(Float64, total_active_nodes)
            while (sum_prob < u[2 * c])
                anc_1 += 1
                sum_prob += convert(Float64, active_nodes[anc_1]) / convert(Float64, total_active_nodes)
            end

            edge[2 * c, 1] = n_index
            edge[2 * c, 2] = anc_1

            edge_length[2 * c] = node_times[anc_1] - coalescence_times[c]

            likelihood *= 2.0 / convert(Float64, total_active_nodes)

            active_nodes[anc_1] = 0
            total_active_nodes -= 1

            anc_2 = 1
            sum_prob = convert(Float64, active_nodes[anc_2]) / convert(Float64, total_active_nodes)
            while (sum_prob < u[2 * c - 1 ])
                anc_2 += 1
                sum_prob += convert(Float64, active_nodes[anc_2]) / convert(Float64, total_active_nodes)
            end
            edge[2 * c - 1, 1] = n_index
            edge[2 * c - 1, 2] = anc_2

            edge_length[2 * c - 1] = node_times[anc_2] - coalescence_times[c]

            likelihood *= 1.0 / convert(Float64, total_active_nodes)

            active_nodes[anc_2] = 0

            active_nodes[n_index] = 1
            node_times[n_index] = coalescence_times[c]

            push!(node_df, [node_times[n_index] n_index anc_1 anc_2])

            c -= 1
            n_index -= 1
        else
            for i in 1:leaves[k]
                active_nodes[l_index] = 1
                total_active_nodes += 1

                node_times[l_index] = leaf_times[k]
                l_index -= 1
            end
            k -= 1
        end
    end
    return edge, edge_length, likelihood, node_df
end



# Benchmarked !!!
function bounded_times_likelihood(leaf_times::Vector{Float64},
                                  leaves::Vector{Int64},
                                  coalescence_times::Vector{Float64},
                                  Nₑ::Float64,
                                  bound::Float64)::Float64
    # Counters
    k = length(leaf_times)
    c = length(coalescence_times)
    # Current time
    current_time = leaf_times[k]
    # Current lineages
    current_lineages = leaves[k]
    k -= 1
    # Next event times
    next_leaf_time = k >= 1 ? leaf_times[k] : bound
    next_coalescence_time = coalescence_times[c]
    likelihood = 1.0
    while current_time > coalescence_times[1]
        if next_leaf_time > next_coalescence_time
            dt = current_time - next_leaf_time
            coef = (convert(Float64, current_lineages) * (convert(Float64,current_lineages) - 1.0)) / (2.0 * Nₑ)
            likelihood *= exp(-coef * dt)
            current_lineages += leaves[k]
            current_time = next_leaf_time
            k -=1
            next_leaf_time = k >= 1 ? leaf_times[k] : bound
        else
            dt = current_time - next_coalescence_time
            coef = (convert(Float64, current_lineages) * (convert(Float64, current_lineages) - 1.0)) / (2.0 * Nₑ)
            likelihood *= coef * exp(-coef * dt)
            current_lineages -= 1
            current_time = next_coalescence_time
            c -=1
            if c >= 1 
                next_coalescence_time = coalescence_times[c]
            end
        end
    end
    total_leaves = sum(leaves)
    # forward_probs = zero(total_leaves, length(leaf_times) + 1)
    forward_probs = forward_algorithm(leaf_times, leaves, Nₑ, bound)
    likelihood /= forward_probs[1, 1]
    return likelihood
end



function sample_wtree(leaf_times::Vector{Float64},
                      leaves::Vector{Int64},
                      Nₑ::Float64,
                      bound::Float64)::Tuple{Matrix{Int64}, Vector{Float64}, Float64, DataFrame}
    coalescence_times, time_likelihood = sample_bounded_times(leaf_times, leaves, Nₑ, bound)
    edge, edge_length, top_likelihood, node_df = sample_topology(leaf_times, leaves, coalescence_times)
    return edge, edge_length, time_likelihood * top_likelihood, node_df
end


"""
    relabel!(wtree::DataFrame)::DataFrame
    Assign unique labels for each node in `wtree`.
"""
function relabel!(wtree::DataFrame)::DataFrame
    new_ids = Dict{Int64, Int64}()
    for node in eachrow(wtree)
        if node.left == 0 && node.right == 0 # Leaf node
            if node.leaf_id == node.host # Sampled leaf
                new_ids[node.id] = CantorPair(node.id, node.host)
                node.id = new_ids[node.id]
            else # Transmission leaf
                new_ids[node.id] = CantorPair(0, node.leaf_id)
                node.id = new_ids[node.id]
            end
        else    # Internal node
            new_ids[node.id] = CantorPair(node.id, node.host)
            node.id = new_ids[node.id]
            node.left = new_ids[node.left]
            node.right = node.right == 0 ? 0 : new_ids[node.right]
        end
    end
    return wtree
end


"""
    normalize!(tree::DataFrame)::DataFrame
    Sequentially label nodes in `tree`.
"""
function normalize!(tree::DataFrame)::DataFrame
    new_labels = Dict(zip(tree.id, 1:nrow(tree)))
    new_labels[0] = 0
    tree.id = 1:nrow(tree)
    tree.left = [new_labels[node] for node in tree.left]
    tree.right = [new_labels[node] for node in tree.right]
    return tree
end


function sample_wtree(host::Host, Nₑ::Float64; relabel=true)
    isnothing(host.infector) && return DataFrame(t = [], id = [], left = [], right = [], leaf_id = [], host = [], type = [])
    leaves = fill(1, length(host.leaf_times))
    tot_leaves = length(leaves)
    if tot_leaves == 1 
        wtree = DataFrame(t = [host.leaf_times[1], host.root], id = [1, 0], left = [0, 1], right = [0, 0])
    else
        _, _, _, wtree = sample_wtree(host.leaf_times, leaves, Nₑ, host.root)
        push!(wtree, (host.root, 0, tot_leaves + 1, 0))
    end
    wtree.leaf_id .= vcat(host.leaf_ids, fill(0, tot_leaves))
    wtree.host .= host.id
    wtree.type .= host.type
    relabel && relabel!(wtree)
    wtree[end, :host] = host.infector
    wtree[end, :type] = host.infector_type
    return wtree
end


function join_wtrees(trees::Vector{DataFrame})::DataFrame
    tree = vcat(trees...)
    tree = tree[(tree.leaf_id .== 0 .|| tree.leaf_id .== tree.host), :] # Remove duplicate nodes
    sort!(tree, :t, rev=true)
    normalize!(tree)
    return tree
end


"""
    get_hosts(linelist::DataFrame) -> Dict{Int64, Host}

    Given a DataFrame `linelist` representing a history of transmission of an outbreak, constructs a dictionary where keys correspond to unique host ids and values are mutable `Host` objects representing each host in the transmission chain. 

    Each row in `linelist` is assumed to represent an infected individual, and should include columns for the infected individual's id (`child_id`), the id of the person that infected them (`parent_id`), the time at which they were infected (`t_birth`) and the time at which they were sampled (`t_sam`), if they were sampled.

    The function processes the rows of `linelist` in reverse order (i.e., from the most recent infections to the earliest), updating the dictionary of `Host` objects to reflect the transmission history. Specifically, for each row in `linelist`, the function checks whether the infectee was sampled (`t_sam > 0.`) or whether the infector already has a record in the `hosts` dictionary (`haskey(hosts, :parent_id)`). If either of these conditions hold, the function creates or updates `Host` objects for both the infector and infectee (as appropriate) by adding or updating information in the `leaf_times` and `leaf_ids` fields. If the infector does not yet have a record in the dictionary, the function creates a new `Host` object with information from the current row of `linelist`.

    The function returns a dictionary of `Host` objects, where each key is an integer representing the id of a unique host in the transmission chain, and each value is a mutable `Host` object with information about that host's place in the transmission chain.

"""
function get_hosts(linelist::DataFrame)
	hosts = Dict{Int64, Host}()
	@eachrow! reverse(linelist) begin
		if :t_sam > 0. || haskey(hosts, :child_id)
			# Check if a record exists for the infector
			if haskey(hosts, :parent_id) # Append existing record
				pushfirst!(hosts[:parent_id].leaf_times, :t_birth)
				pushfirst!(hosts[:parent_id].leaf_ids, :child_id)
			else   # Create new record
				hosts[:parent_id] = Host(:parent_id, :parent_type, nothing, 0, NaN, [:t_birth], [:child_id])
			end

			if :t_sam > 0. && !haskey(hosts, :child_id)	# Sampled infectee
				hosts[:child_id] = Host(:child_id, :child_type, :parent_id, :parent_type, :t_birth, [:t_sam], [:child_id])	# Create a new record for the infectee
			
			elseif haskey(hosts, :child_id)	# Infector
				hosts[:child_id].root = :t_birth	# Update root
                hosts[:child_id].infector = :parent_id  # Update infector
                hosts[:child_id].type = :child_type # Update type
                hosts[:child_id].infector_type = :parent_type
				
				if :t_sam > 0.  # TODO: Generalize to sampling without removal (i.e., the ordered placement of leaves)
					push!(hosts[:child_id].leaf_times, :t_sam)	# Update leaves to include sampling
					push!(hosts[:child_id].leaf_ids, :child_id)
				end
			end
		end
	end
	return hosts
end


function simulate(linelist::DataFrame, Nₑ::Float64)::DataFrame
    hosts = get_hosts(linelist)
    # Generate within-host trees for each sampled host in linelist
    wtrees = [sample_wtree(host, Nₑ) for host in values(hosts) if host.id != 0]
    return join_wtrees(wtrees)
end


function simulate(proc::BDProcess, Nₑ::Float64)::DataFrame
    return simulate(proc.linelist, Nₑ)
end
