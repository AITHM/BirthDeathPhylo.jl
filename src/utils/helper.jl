
"""
    int(x)

    Return ⌊x⌋.
"""
int(x) = convert(Int64, trunc(x))


"""
    logit(x)

    Logistic transform of `x`.

"""
function logit(x)
    return log(x / (1 - x))
end


"""
    inv_logit(x)

    Inverse logistic transform of `x`
"""
function inv_logit(x)
    return 1. / (exp(-x) + 1.)
end


"""
    sortinsert!(v::Vector, c)::Vector

    Insert `c` in the sorted vector `v` such that ordering is preserved 
"""
function sortinsert!(v::Vector, c; rev=false)::Vector
    insert!(v, searchsortedfirst(v, c, rev=rev), c)
    return v
end


"""
    extend!(d::Dict, k, v)::Dict

    Extend dictionary `d` by either adding new `(k, v)` pair, or pushing `v` to `d[k]`
"""
function extend!(d::Dict, k, v)::Dict
    if k in keys(d)
        push!(d[k], v)
    else
        d[k] = [v]
    end
    return d
end

"""
    sortextend!(d::Dict, k, v)::Dict

    Extend dictionary `d` by either adding new `(k,v)` pair, or sort-inserting `v` to `d[k]`
"""
function sortextend!(d::Dict, k, v; rev=false)::Dict
    if k in keys(d)
        sortinsert!(d[k], v, rev=rev)
    else
        d[k] = [v]
    end
    return d
end


"""
    increment!(d::Dict, k, v)::Dict

    Increment key `k` of dictionary `d` by value `v`. If key `k` is not in `d`, create new entry `d[k] = v`.
"""
function increment!(d::Dict, k, v)::Dict
    if k in keys(d)
        d[k] += v
    else
        d[k] = v
    end
    return d
end


"""
    table(x::Vector{Int64})::Vector{Int64}

    Compute frequency distribution of elements in `x`
"""
function table(x::Vector{Int64})::Vector{Int64}
    out = fill(0, maximum(x) + 1)
    for val in x
       out[val + 1] += 1
    end
    return out
end


function table(x::Vector{T}) where T <: Any
    weights = fill(1, length(x))
    for idx in eachindex(x)
        for jdx in 1:idx-1
            if x[jdx] == x[idx]
                weights[jdx] += 1
                weights[idx] -= 1
            end
        end
    end
    return weights
end


"""
    CantorPair(x,y)
    Calculate the Cantor Pair mapping of `x` and `y`.
"""
function CantorPair(x, y)
    y == 0 && return 0
    return trunc(Int, (x^2 + 3 * x + 2 * x * y + y + y^2) / 2)
end