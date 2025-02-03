function bottom_bump!(p::Problem, h, λ, x₀ = 0)
    β̂, κ, ℓ = p.β̂, p.κ, p.ℓ
    β̂[:] = @. h * λ * √(2π) / 4ℓ * exp(- λ^2 * κ^2 / 32) * exp(-im * κ * x₀)
    return nothing
end
