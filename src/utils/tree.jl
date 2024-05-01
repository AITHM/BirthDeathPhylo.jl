function isleaf(node::DataFrameRow)::Bool
    return node.left == 0 && node.right == 0
end

function isbinary(node::DataFrameRow)::Bool
    return node.right != 0
end

function isunary(node::DataFrameRow)::Bool
    return node.left != 0 && node.right == 0
end


function isroot(node::DataFrameRow)::Bool
    return node.host == 0
end

function binarize(tree::DataFrame)
    ctree = copy(tree)
    binary_tree = DataFrame()
    for row in eachrow(ctree)
        if isleaf(row)
            push!(binary_tree, row)
        elseif isbinary(row)
            left = row.left
            right = row.right
            while isunary(tree[left, :])
                left = tree[left, :left]
            end
            while isunary(tree[right, :])
                right = tree[right, :left]
            end
            row.left = left
            row.right = right
            push!(binary_tree, row)
        elseif isroot(row)
            left = row.left
            while isunary(tree[left, :])
                left = tree[left, :left]
            end
            row.left = left
            push!(binary_tree, row)
        end
    end
    new_ids = Dict(zip(binary_tree.id, 1:nrow(binary_tree)))
    new_ids[0] = 0
    binary_tree.id = [new_ids[id] for id in binary_tree.id]
    binary_tree.left = [new_ids[id] for id in binary_tree.left]
    binary_tree.right = [new_ids[id] for id in binary_tree.right]
    return binary_tree
end


function update_branch_lengths!(tree::DataFrame)
    lengths = Dict{Int64, Float64}()
    lengths[nrow(tree)] = 0.0
    @eachrow! reverse(tree) begin
        if :left != 0
            lengths[:left] = tree[:left, :].t - :t
        end
        if :right != 0
            lengths[:right] = tree[:right, :].t - :t
        end
    end
    # for row in eachrow(reverse(tree))
    #     if row.left != 0
    #         lengths[row.left] = tree[row.left, :t] - row.t
    #     end
    #     if row.right != 0
    #         lengths[row.right] = tree[row.right, :t] - row.t
    #     end
    # end
    tree.length = [lengths[id] for id in tree.id]
    return
end


function bound_nodes(tree::DataFrame)
    bounds = fill(Inf, nrow(tree))
    @eachrow! tree begin
        if :left == 0
            bounds[:id] = :t
        elseif :right != 0
            bounds[:id] = minimum([bounds[:left], bounds[:right]])
        else
            bounds[:id] = bounds[:left]
        end
    end
    return bounds
end


function isvalid(tree::DataFrame)::Bool
    for node in eachrow(tree)
        if node.left != 0
            node.t > tree[node.left, :t] && return false
            if node.right != 0
                node.t > tree[node.right, :t] && return false
            end
        end
    end
    return true
end


function get_topology(tree::DataFrame)::Matrix{Int64}
    topology = zeros(3, nrow(tree))
    @eachrow! tree begin
        if :left != 0
            topology[1,:left] = :id
            topology[2, :id] = :left
            if :right != 0
                topology[1, :right] = :id
                topology[3, :id] = :right
            end
        end
    end
    return topology
end


function make_tree(topology::Matrix{Int64}, t_leaves::Vector{Float64}, t_nodes::Vector{Float64})
    tree = DataFrame(t = NaN, id = 1:size(topology,2), left = 0, right = 0)
    @eachrow! tree begin
        tree[:id, :].left = topology[2, :id]
        tree[:id, :].right = topology[3, :id]
        tree[:id, :].t = topology[2, :id] == 0 ? popfirst!(t_leaves) : popfirst!(t_nodes)
    end
    add_length!(tree)
    return tree
end


function randomize(tree::DataFrame, scale::Float64)
    rand_tree = copy(tree)
    @eachrow! rand_tree begin
        if :left != 0 && :right != 0
            rand_tree[:id, :].t = minimum([rand_tree[:left, :].t, rand_tree[:right, :].t]) - rand(Exponential(scale))
        end
    end
    rand_tree[end, :t] = rand_tree[rand_tree[end, :left], :t] - rand(Exponential(scale))
    update_branch_lengths!(rand_tree)
    return rand_tree
end


function get_ancestry(tree::DataFrame)
    ancestry = fill(0, n_nodes)
    @eachrow! tree begin
        if :left != 0
            ancestry[:left] = :id
            if :right != 0
                ancestry[:right] = :id
            end
        end
    end
    return ancestry
end


function get_tip_labels(tree::DataFrame)::Vector{String}
    tip_labels = Vector{String}()
    @eachrow! tree begin
        if :leaf_id != 0
            push!(tip_labels, join([:id, "_", :t]))
        end
    end
    return tip_labels
end


function newick(tree::DataFrame)
    function newick(node::DataFrameRow)
        if isleaf(node)
            return "$(node.id)_$(node.t):$(node.length)"
        else
            left_str = newick(tree[node.left, :])
            right_str = newick(tree[node.right, :])
            return "($left_str, $right_str):$(node.length)"
        end
    end
    return readnw(chop(newick(tree[end - 1, :]), tail=4)*";")
end


# TODO: Work appropriate separator between branch length and sequence
function newick(tree::DataFrame, sequences::Dict{Int, String})
    function newick(node::DataFrameRow, sequences::Dict{Int, String})
        if isleaf(node)
            sequence = sequences[node.id]
            return "$(node.id)_$(node.t):$(node.length)$sequence"
        else
            left_str = newick(tree[node.left, :], sequences)
            right_str = newick(tree[node.right, :], sequences)
            return "($left_str, $right_str):$(node.length)"
        end
    end
    return chop(newick(tree[end - 1, :], sequences), tail=4)
end


function nexus(tree::DataFrame, seqs::Dict{Int, String}, filename::AbstractString)
    tip_sequences = Dict{String, String}()
    for node in eachrow(tree)
        if node.left == 0
            tip_sequences["$(node.id)_$(node.t)"] = seqs[node.id]
        end
    end
    num_sequences = length(tip_sequences)
    sequence_length = length(first(values(tip_sequences)))

    # Open the file for writing
    io = open(filename, "w")

    # Write the Nexus header
    write(io, "#NEXUS\n")
    write(io, "BEGIN TAXA;\n")
    write(io, "\tDIMENSIONS NTAX=$num_sequences;\n")
    write(io, "\t$(join(keys(tip_sequences), ' '));")
    write(io, "END;")

    write(io, "BEGIN CHARACTERS;\n")
    write(io, "\tDIMENSIONS NCHAR=$sequence_length;\n")
    write(io, "\tFORMAT\n")
    write(io, "\t\tDATATYPE=DNA\n\n")
    write(io, "\t\tGAP=-\n")
    write(io, "\t\tMISSING=?\n")



    write(io, "BEGIN DATA;\n")
    write(io, "DIMENSIONS NTAX=$num_sequences NCHAR=$sequence_length;\n")
    write(io, "FORMAT DATATYPE=DNA GAP=- MISSING=?;\n")
    write(io, "MATRIX\n")

    # Write sequence data
    for (id, sequence) in tip_sequences
        write(io, "$id $sequence\n")
    end

    write(io, ";\n")
    write(io, "END;\n")

    # Write the tree
    write(io, "BEGIN TREES;\n")
    # write(io, "TRANSLATE\n")
    
    # # Write translation table for node IDs to sequence IDs
    # for (i, id) in enumerate([node.id for node in eachrow(tree) if node.left == 0])
    #     tree[i, :left] == 0 && write(io, " $i $id,\n")
    # end
    
    # write(io, ";\n")
    write(io, "TREE tree = ")
    
    # Generate Newick string with tip sequences
    newick_string = newick(tree)

    write(io, "$newick_string;\n")
    write(io, "END;\n")

    # Close the file
    close(io)
end



function make_tree(df::DataFrame)
    tree = RecursiveTree{OneRoot, String, Dict{String, Any}, Dict{String, Any}, PolytomousBranching, Float64,Dict{String, Any}}()
    for row in eachrow(df)
        lab = string(row.id)
        createnode!(tree, lab)
        setnodedata!(tree, lab, "t", row.t)
        setnodedata!(tree, lab, "host", row.host)
        setnodedata!(tree, lab, "type", row.type)
        row.left != 0 && createbranch!(tree, lab, string(row.left), df[row.left, :t] - row.t)
        row.right != 0 && createbranch!(tree, lab, string(row.right), df[row.right, :t] - row.t)
    end
    return tree
end