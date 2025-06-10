"""
    bottom_bump!(p::Problem, h, λ, x₀ = 0)

Calculate `β̂` coefficients of a problem `p` for a bottom bump of height `h`
and characteristic length `λ` at a position `x₀`.
"""
function bottom_bump!(p::Problem, h, λ, x₀ = 0)
    β̂, κ, ℓ = p.β̂, p.κ, p.ℓ
    β̂[:] = @. h * λ * √(2π) / 4ℓ * exp(- λ^2 * κ^2 / 32) * exp(-im * κ * x₀)
    return nothing
end

"""
    bottom_step!(p::Problem, h)

Calculate `β̂` coefficients of a wave problem `p` for a bottom step of height `h`.
"""
function bottom_step!(p::Problem, h)
    β̂, κ, ℓ, ℐ = p.β̂, p.κ, p.ℓ, p.ℐ
    β̂[:] = @. -h / 2 * sinc(κ * ℓ / 4π)
    β̂[ℐ+1] = h / 2
    return nothing
end

"""
    bottom_slope!(p::Problem, h)

Calculate `β̂` coefficients of a wave problem `p` for a bottom slope of height `h`.

"""
function bottom_slope!(p::Problem, h)
    β̂, κ, ℓ, ℐ = p.β̂, p.κ, p.ℓ, p.ℐ
    β̂[:] = @. -4h / 3 * sinc(κ * ℓ / 3π)^2
    β̂[ℐ+1] = 2h / 3
    return nothing
end

"""
    bottom_bar!(p::Problem, h)

Calculate `β̂` coefficients of a wave problem `p` for a bottom bar of height `h`.

"""
function bottom_bar!(p::Problem, h)
    β̂, κ, ℓ, ℐ = p.β̂, p.κ, p.ℓ, p.ℐ
    β̂[:] = @. -4h / 3 * sinc(κ * ℓ / 3π)^2
    β̂[ℐ+1] = 2h / 3
    return nothing
end

"""
    bottom_vector!(p::Problem, x::AbstractRange{<:Real}, β::Vector{<:Real})

Calculate `β̂` coefficients of a wave problem `p` for a bottom profile `β(x)`.

"""
function bottom_vector!(p::Problem, x::AbstractRange{<:Real}, β::Vector{<:Real})
    β̂, κ = p.β̂, p.κ
    for i in eachindex(κ)
        β̂[i] = fourier_transform(β, κ[i], x)
    end
    return nothing
end

"""
    moving_bottom_bump!(p::Problem, h, λ, u, x₀ = 0)

Calculate `β̂` coefficients of a wave problem `p` for a moving bottom bump of height `h`,
characteristic length `λ`, and velocity `u` at initial position `x₀`.
"""
function moving_bottom_bump!(p::Problem, h, λ, u, x₀ = 0)
    β̂, β̇, κ, ℓ, t, O, N = p.β̂, p.β̇, p.κ, p.ℓ, p.t, p.O, p.N
    for n in O:N+O-1
        β̂[:, n] = @. h * λ * √(2π) / 4ℓ * exp(- λ^2 * κ^2 / 32) * exp(-im * κ * x₀) * exp(-im * κ * u * t[n-O+1])
        β̇[:, n] = @. -im * κ * u * β̂[:, n]
    end
    return nothing
end
