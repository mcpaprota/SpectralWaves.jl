"""
    init_problem(ℓ::Number, d::Number, ℐ::Integer, N::Integer; O = 4)

Initialize IBVP wave problem corresponding to a fluid domain of length `ℓ` and depth `d`
with `ℐ` harmonics and `N` time steps.

Output is a tuple `(η̂, η̇, β̂, β̇, ϕ̂, ϕ̇, ψ̂, ψ̇, p̂, χ, ξ, ζ, 𝒯, 𝒮)`, where:
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
- `𝒮` are hyperbolic secant lookup values.
"""
function init_problem(ℓ::Number, d::Number, ℐ::Integer, N::Integer; O = 4)
    κ = 2π / ℓ * (-ℐ:ℐ)
    η̂ = complex(zeros(2ℐ + 1, N))
    η̇ = complex(zeros(2ℐ + 1, N))
    β̂ = complex(zeros(2ℐ + 1, N))
    β̇ = complex(zeros(2ℐ + 1, N))
    ϕ̂ = complex(zeros(2ℐ + 1, N))
    ϕ̇ = complex(zeros(2ℐ + 1, N))
    ψ̂ = complex(zeros(2ℐ + 1, N))
    ψ̇ = complex(zeros(2ℐ + 1, N))
    p̂ = complex(zeros(2ℐ + 1, N))
    χ = zeros(N)
    ξ = zeros(N)
    ζ = zeros(N)
    𝒯 = tanh.(κ * d)
    𝒮 = sech.(κ * d)
    return κ, η̂, η̇, β̂, β̇, ϕ̂, ϕ̇, ψ̂, ψ̇, p̂, χ, ξ, ζ, 𝒯, 𝒮, O
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
        Φ̃″[:, m + 1] = (iseven(m) ? zero(κ) : κ .^ m .* S) * im
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
- `Ψ̂″` are bottom-potential-amplitude dependent expansion coefficients for computing
bottom velocity potential amplitudes and its vertical gradients,
- `Ψ̃′` are surface-potential-amplitude dependent expansion coefficients and
- `Ψ̃″` are bottom-potential-amplitude dependent expansion coefficients for computing
bottom horizontal velocity potential amplitudes and its vertical gradients.

"""
function init_nonlinear_bottom_boundary_condition(κ, 𝒯, 𝒮, ℐ, M)
    Ψ̂′ = zeros(2ℐ + 1, M + 1)
    Ψ̂″ = zeros(2ℐ + 1, M + 1)
    Ψ̃′ = complex(zeros(2ℐ + 1, M + 1))
    Ψ̃″ = complex(zeros(2ℐ + 1, M + 1))
    A′ = zeros(2ℐ + 1, 2ℐ + 1)
    A″ = zeros(2ℐ + 1, 2ℐ + 1)
    for m in 0:M
        Ψ̂′[:, m + 1] = iseven(m) ? zero(κ) : κ .^ (m + 1) .* 𝒮
        Ψ̂″[:, m + 1] = iseven(m) ? κ .^ m : -(κ .^ m) .* 𝒯
        Ψ̃′[:, m + 1] = (iseven(m) ? κ .^ (m + 1) .* 𝒮 : zero(κ)) * im
        Ψ̃″[:, m + 1] = (iseven(m) ? -(κ .^ m) .* 𝒯 : κ .^ m) * im
    end
    return Ψ̂′, Ψ̂″, Ψ̃′, Ψ̃″, A′, A″
end
