# SPDX-License-Identifier: MIT

"""
    Problem(ℓ::Number, d::Number, ℐ::Integer, t::AbstractRange{<:Number}; O = 4, M_s = 0, M_b = 0, static_bottom = true)

Construct an IBV Problem object corresponding to a fluid domain of length `ℓ` and depth `d`
with `ℐ` harmonics and `N` time steps.

Output is a Problem object with fields:
- `ℓ` is the fluid domain length (m),
- `d` is the water depth (m),
- `ℐ` is the number of harmonics,
- `t` is the time range,
- `Δt` is the time step (s),
- `N` is the number of time steps,
- `O` is the order of the time-stepping scheme,
- `κ` are wave numbers (rad/m),
- `η̂` are free-surface elevation amplitudes (m),
- `η̇` are free-surface vertical velocity amplitudes (m/s),
- `β̂` are bottom-surface elevation amplitudes (m),
- `β̇` are bottom-surface vertical velocity amplitudes (m/s),
- `ϕ̂` are flat-bottom velocity potential amplitudes (m²/s),
- `ϕ̇` are flat-bottom acceleration potential amplitudes (m²/s²),
- `ψ̂` are uneven-bottom velocity potential amplitudes (m²/s),
- `ψ̇` are uneven-bottom acceleration potential amplitudes (m²/s²),
- `p̂` are surface pressure head amplitudes (m),
- `χ` is wavemaker paddle displacement (m),
- `ξ` is wavemaker paddle velocity (m/s),
- `ζ` is wavemaker paddle acceleration (m/s²),
- `𝒯` are hyperbolic tangent lookup values,
- `𝒮` are hyperbolic secant lookup values,
- `static_bottom` is a boolean flag to indicate if the bottom is static.

"""
struct Problem
    ℓ::Number
    d::Number
    ℐ::Integer
    t::AbstractRange{<:Number}
    Δt::Number
    N::Integer
    O::Integer
    M_s::Integer
    M_b::Integer
    F::Vector{Number}
    κ::AbstractRange{<:Number}
    η̂::Matrix{Complex}
    η̇::Matrix{Complex}
    β̂::Matrix{Complex}
    β̇::Matrix{Complex}
    ϕ̂::Matrix{Complex}
    ϕ̇::Matrix{Complex}
    ψ̂::Matrix{Complex}
    ψ̇::Matrix{Complex}
    p̂::Matrix{Complex}
    χ::Vector{Number}
    ξ::Vector{Number}
    ζ::Vector{Number}
    𝒯::Vector{Number}
    𝒮::Vector{Number}
    static_bottom::Bool
    function Problem(ℓ::Number, d::Number, ℐ::Integer, t::AbstractRange{<:Number}; O = 4, M_s = 0, M_b = 0, static_bottom = true)
        N = length(t) - 1
        Δt = step(t)
        F = factorial_lookup(max(M_s, M_b))
        κ = 2π / ℓ * (-ℐ:ℐ)
        η̂ = complex(zeros(2ℐ + 1, N+O))
        η̇ = complex(zeros(2ℐ + 1, N+O))
        if static_bottom
            β̂ = complex(zeros(2ℐ + 1, 1))
        else
            β̂ = complex(zeros(2ℐ + 1, N+O))
        end
        β̇ = complex(zeros(2ℐ + 1, N+O))
        ϕ̂ = complex(zeros(2ℐ + 1, N+O))
        ϕ̇ = complex(zeros(2ℐ + 1, N+O))
        ψ̂ = complex(zeros(2ℐ + 1, N+O))
        ψ̇ = complex(zeros(2ℐ + 1, N+O))
        p̂ = complex(zeros(2ℐ + 1, N+O))
        χ = zeros(N+O)
        ξ = zeros(N+O)
        ζ = zeros(N+O)
        𝒯 = tanh.(κ * d)
        𝒮 = sech.(κ * d)
        new(ℓ, d, ℐ, t, Δt, N, O, M_s, M_b, F, κ, η̂, η̇, β̂, β̇, ϕ̂, ϕ̇, ψ̂, ψ̇, p̂, χ, ξ, ζ, 𝒯, 𝒮, static_bottom)
    end
end

"""
    init_nonlinear_surface_boundary_condition(κ, 𝒯, 𝒮, ℐ, M)

Initialize expansion coefficients for nonlinear free-surface boundary conditions
for eigenvalues `κ`, hyperbolic tangent `𝒯` and secant `𝒮` values,
number of harmonics `ℐ` and order of nonlinear expansion `M`.

Output is a tuple `(Φ̇′, Φ̇″, Φ̂′, Φ̂″, Φ̃′, Φ̃″)`, where:
- `Φ̇′` are surface-potential-amplitude dependent expansion coefficients and
- `Φ̇″` are bottom-potential-amplitude dependent expansion coefficients for computing
surface acceleration potential amplitudes and its vertical gradients,
- `Φ̂′` are surface-potential-amplitude dependent expansion coefficients and
- `Φ̂″` are bottom-potential-amplitude dependent expansion coefficients for computing
surface velocity potential amplitudes and its vertical gradients,
- `Φ̃′` are surface-potential-amplitude dependent expansion coefficients and
- `Φ̃″` are bottom-potential-amplitude dependent expansion coefficients for computing
surface horizontal velocity potential amplitudes and its vertical gradients.

"""
function init_nonlinear_surface_boundary_condition(κ, 𝒯, 𝒮, ℐ, M)
    Φ̇′ = zeros(2ℐ + 1, M + 1)
    Φ̇″ = zeros(2ℐ + 1, M + 1)
    Φ̂′ = zeros(2ℐ + 1, M + 1)
    Φ̂″ = zeros(2ℐ + 1, M + 1)
    Φ̃′ = complex(zeros(2ℐ + 1, M + 1))
    Φ̃″ = complex(zeros(2ℐ + 1, M + 1))
    for m in 0:M
        Φ̇′[:, m + 1] = iseven(m) ? κ .^ m : κ .^ m .* 𝒯
        Φ̇″[:, m + 1] = iseven(m) ? zero(κ) : κ .^ (m - 1) .* 𝒮
        Φ̂′[:, m + 1] = iseven(m) ? κ .^ (m + 1) .* 𝒯 : κ .^ (m + 1)
        Φ̂″[:, m + 1] = iseven(m) ? κ .^ m .* 𝒮 : zero(κ)
        Φ̃′[:, m + 1] = (iseven(m) ? κ .^ (m + 1) : κ .^ (m + 1) .* 𝒯) * im
        Φ̃″[:, m + 1] = (iseven(m) ? zero(κ) : κ .^ m .* 𝒮) * im
    end
    return Φ̇′, Φ̇″, Φ̂′, Φ̂″, Φ̃′, Φ̃″
end

"""
    init_nonlinear_bottom_boundary_condition(κ, 𝒯, 𝒮, ℐ, M)

Initialize expansion coefficients for nonlinear bottom boundary conditions
for eigenvalues `κ`, hyperbolic tangent `𝒯` and secant `𝒮` values,
number of harmonics `ℐ` and order of nonlinear expansion `M`.

Output is a tuple `(Ψ̂′, Ψ̂″, Ψ̃′, Ψ̃″)`, where:
- `Ψ̂′` are surface-potential-amplitude dependent expansion coefficients and
- `Ψ̂″` are bottom-potential-amplitude dependent expansion coefficients for computing bottom velocity potential amplitudes and its vertical gradients,
- `Ψ̃′` are surface-potential-amplitude dependent expansion coefficients and
- `Ψ̃″` are bottom-potential-amplitude dependent expansion coefficients for computing bottom horizontal velocity potential amplitudes and its vertical gradients.

"""
function init_nonlinear_bottom_boundary_condition(κ, 𝒯, 𝒮, ℐ, M)
    A′ = complex(zeros(2ℐ + 1, 2ℐ + 1))
    A″ = complex(zeros(2ℐ + 1, 2ℐ + 1))
    Ψ̂′ = zeros(2ℐ + 1, M + 1)
    Ψ̂″ = zeros(2ℐ + 1, M + 1)
    Ψ̃′ = complex(zeros(2ℐ + 1, M + 1))
    Ψ̃″ = complex(zeros(2ℐ + 1, M + 1))
    w′ = complex(zeros(4ℐ + 1))
    for m in 0:M
        Ψ̂′[:, m+1] = iseven(m) ? zero(κ) : κ .^ (m + 1) .* 𝒮
        Ψ̂″[:, m+1] = iseven(m) ? κ .^ m : -(κ .^ m) .* 𝒯
        Ψ̃′[:, m+1] = (iseven(m) ? κ .^ (m + 1) .* 𝒮 : zero(κ)) * im
        Ψ̃″[:, m+1] = (iseven(m) ? -(κ .^ m) .* 𝒯 : κ .^ m) * im
    end
    return A′, A″, Ψ̂′, Ψ̂″, Ψ̃′, Ψ̃″, w′
end
