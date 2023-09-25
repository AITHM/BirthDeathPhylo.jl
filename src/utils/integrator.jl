using FastBroadcast
using StaticArrays

# function rk4_general!(k, y, t, f, h)
#     k[1] .= f(t, y)
#     k[2] .= f(t + c2 * h, y .+ h * (a21 * k[1]))
#     k[3] .= f(t + c3 * h, y .+ h * (a31 * k[1] .+ a32 * k[2]))
#     k[4] .= f(t + c4 * h, y .+ h * (a41 * k[1] .+ a42 * k[2] .+ a43 * k[3]))
#     k[5] .= f(t + c5 * h, y .+ h * (a51 * k[1] .+ a52 * k[2] .+ a53 * k[3] .+ a54 * k[4]))
#     return k
# end


function rungekutta4(f, y0, t)
    n = length(t)
    y = zeros((n, length(y0)))
    y[1,:] = y0
    for i in 1:n-1
        h = t[i+1] - t[i]
        k1 = f(y[i,:], t[i])
        k2 = f(y[i,:] + k1 * h/2, t[i] + h/2)
        k3 = f(y[i,:] + k2 * h/2, t[i] + h/2)
        k4 = f(y[i,:] + k3 * h, t[i] + h)
        y[i+1,:] = y[i,:] + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
    end
    return y
end



@inline function dp_coeffs!(k, y, t, f, h)
    @views begin
        @.. k[:, 1] = f(t, y)
        k[:, 2] .= f(t + 1/5 * h, y .+ h * (1/5 * k[:, 1]))
        @.. k[:, 3] = f(t + 3/10 * h, y + h * (3/40 * k[:, 1] + 9/40 * k[:, 2]))
        @.. k[:, 4] = f(t + 4/5 * h, y + h * (44/45 * k[:, 1] - 56/15 * k[:, 2] + 32/9 * k[:, 3]))
        @.. k[:, 5] = f(t + 8/9 * h, y + h * (19372/6561 * k[:, 1] - 25360/2187 * k[:, 2] + 64448/6561 * k[:, 3] - 212/729 * k[:, 4]))
        @.. k[:, 6] = f(t + h, y + h * (9017/3168 * k[:, 1] - 355/33 * k[:, 2] + 46732/5247 * k[:, 3] + 49/176 * k[:, 4] - 5103/18656 * k[:, 5]))
        @.. k[:, 7] = f(t + h, y + h * (35/384 * k[:, 1] + 500/1113 * k[:, 3] + 125/192 * k[:, 4] - 2187/6784 * k[:, 5] + 11/84 * k[:, 6]))
    end
    return k
end


@inline function dp_update(y_next, y, t, f, h, k, abs_tol, rel_tol, b_high = [35/384, 0., 500/1113, 125/192, -2187/6784, 11/84, 0.], b_low = [5179/57600, 0., 7571/16695, 393/640, -92097/339200, 187/2100, 1/40])
    
    dp_coeffs!(k, y, t, f, h)
    t_next = t + h
    y_next .= y .+ h * k * b_high
    
    abs_err = h * k * (b_high .- b_low)

    # scale = abs_tol + rel_tol * maximum(norm.([y_next, y]))
    scale = abs_tol + rel_tol * my_norm2(y_next)
    err = my_norm(abs_err / scale)
    # err = err > 1e-10 ? err : 1e-10

    h *= (1. / err) ^ 0.2

    if err ≤ 1
        @.. y = y_next
        return t_next, h
    end

    dp_update(y_next, y, t, f, h, k, abs_tol, rel_tol)

end


function integrate(y, f, t_init, t_final, h, abs_tol=1e-5, rel_tol=1e-5)
    k = zeros(length(y), 7)
    y_next = similar(y)
    t = t_init
    h_max = t_final - t
    while h_max > 0.
        h = h < h_max ? h : h_max
        t, h = dp_update(y_next, y, t, f, h, k, abs_tol, rel_tol)
        h_max = t_final - t
    end
    return y, t
end


function integrate_sa(y, f, t_init, t_final, h, abs_tol=1e-5, rel_tol=1e-5)
    k = @MMatrix zeros(length(y), 7)
    y_next = @MVector zeros(length(y))
    t = t_init
    h_max = t_final - t
    while h_max > 0.
        h = h < h_max ? h : h_max
        t, h = dp_update(y_next, y, t, f, h, k, abs_tol, rel_tol)
        h_max = t_final - t
    end
    return y, t
end


@inline function my_norm(u::Array{T}) where T <: AbstractFloat
    x = zero(T)
    @inbounds @fastmath for ui in u
        x += abs2(ui)
    end
    Base.FastMath.sqrt_fast(x / length(u))
end


@inline function my_norm2(u::Array{T}) where T <: AbstractFloat
    x = zero(T)
    @inbounds @fastmath for ui in u
        x += abs2(ui)
    end
    Base.FastMath.sqrt_fast(x)
end


@inline function my_norm(u::SVector)
    x = 0.
    @inbounds @fastmath for ui in u
        x += abs2(ui)
    end
    Base.FastMath.sqrt_fast(x / length(u))
end




@inline function ODE_DEFAULT_NORM(u::Array{T}, t) where {T <: Union{AbstractFloat, Complex}}
    x = zero(T)
    @inbounds @fastmath for ui in u
        x += abs2(ui)
    end
    Base.FastMath.sqrt_fast(real(x) / max(length(u), 1))
end



function k_vecs!(k1, k2, k3, k4)
    k1 = rand(10_000)
    k2 = 2. * k1
    k3 = 2. * k1 .+ k2
    k4 = 1. * k1 .+ k2 .+ k3
    return k1, k2, k3, k4
end


# function k_mat!(k)
#     k[:,1] .= rand(10_000)
#     k[:,2] .= 2. * view(k, :, 2)
#     k[:,3] .= 2. * view(k, :, 1) .+ view(k,:,2)
#     k[:,4] .= 2. * view(k,:,1) .+ view(k,:,2) .+ view(k,:,3)
#     return k
# end

function k_mat!(k)
    @views begin
    k[:,1] .= rand(10_000)
    k[:,2] .= 2. * k[:,2]
    k[:,3] .= 2. * k[:,1] .+ k[:,2]
    k[:,4] .= 2. * k[:,1] .+ k[:,2] .+ k[:,3]
    end
    return k
end