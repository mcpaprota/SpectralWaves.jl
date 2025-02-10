function water_surface(p::Problem, x::Real, n::Integer)
    η̂, κ, O = p.η̂, p.κ, p.O
    η = inverse_fourier_transform(η̂[:, n+O-1], κ, x)
    return η
end

function bottom_surface(p::Problem, x::Real, n=1)
    β̂, κ = p.β̂, p.κ
    β = inverse_fourier_transform(β̂[:, n], κ, x)
    return β
end

function water_velocity(p::Problem, x::Real, z::Real, n::Integer, c::Integer)
    ϕ̂, ψ̂, κ, ℐ, O, d = p.ϕ̂, p.ϕ̃, p.ψ̂, p.ψ̃, p.κ, p.ℐ, p.O, p.d
    if c==1
        û = @. ϕ̂[:, n+O-1] * im * κ * cosh(κ * (z + d)) / cosh(κ * d) +
            ψ̂[:, n+O-1] * im * sinh(κ * z) / cosh(κ * d)
        u = inverse_fourier_transform(û, κ, x)
        return u
    elseif c==3
        ŵ = @. ϕ̃[:, n+O-1] * im * κ * sinh(κ * (z + d)) / cosh(κ * d) +
            ψ̃[:, n+O-1] * im * cosh(κ * z) / cosh(κ * d)
        w = inverse_fourier_transform(ŵ, κ, x)
        return w
    end
end
