function water_surface(p::Problem, x::Real, n::Integer)
    η̂, κ, O = p.η̂, p.κ, p.O
    η = inverse_fourier_transform(η̂[:, n+O-1], κ, x)
    return η
end

function bottom_surface(p::Problem, x::Real, n=1)
    β̂, κ, O = p.β̂, p.κ, p.O
    β = inverse_fourier_transform(β̂[:, n+O-1], κ, x)
    return β
end
