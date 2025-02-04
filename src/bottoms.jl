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
    β̂, κ, ℓ = p.β̂, p.κ, p.ℓ
    β̂[:] = @. -h / 2 * sinc(κ * ℓ / 4π)
    β̂[ℐ+1] = h / 2
    return nothing
end

"""
    bottom_slope!(p::Problem, h)

Calculate `β̂` coefficients a wave problem `p` for a bottom slope of height `h`.

"""
function bottom_slope!(p::Problem, h)
    β̂, κ, ℓ = p.β̂, p.κ, p.ℓ
    β̂[:] = @. -4h / 3 * sinc(κ * ℓ / 3π)^2
    β̂[ℐ+1] = 2h / 3
    return nothing
end
