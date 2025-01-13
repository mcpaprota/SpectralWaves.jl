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
function init_problem(ℓ::Number, d::Number, ℐ::Integer, N::Integer)
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
    return κ, η̂, η̇, β̂, β̇, ϕ̂, ϕ̇, ψ̂, ψ̇, p̂, χ, ξ, ζ, 𝒯, 𝒮
end

function init_nonlinear_surface_problem(κ, 𝒯, 𝒮, ℐ, M)
    Φ̇′ = zeros(2ℐ + 1, M + 1)
    Φ̇″ = zeros(2ℐ + 1, M + 1)
    Φ̂′ = zeros(2ℐ + 1, M + 1)
    Φ̂″ = zeros(2ℐ + 1, M + 1)
    Φ̃′ = complex(zeros(2ℐ + 1, M + 1))
    Φ̃″ = complex(zeros(2ℐ + 1, M + 1))
    for m in 0:M
        Φ̇′[:, m + 1] = iseven(m) ? zero(κ) : κ .^ (m - 1) .* 𝒮
        Φ̇″[:, m + 1] = iseven(m) ? κ .^ m : κ .^ m .* 𝒯
        Φ̂′[:, m + 1] = iseven(m) ? κ .^ m .* 𝒮 : zero(κ)
        Φ̂″[:, m + 1] = iseven(m) ? κ .^ (m + 1) .* 𝒯 : κ .^ (m + 1)
        Φ̃′[:, m + 1] = (iseven(m) ? zero(κ) : κ .^ m .* S) * im
        Φ̃″[:, m + 1] = (iseven(m) ? κ .^ (m + 1) : κ .^ (m + 1) .* 𝒯) * im
    end
    return Φ̇′, Φ̇″, Φ̂′, Φ̂″, Φ̃′, Φ̃″
end

function init_nonlinear_bottom_problem(κ, 𝒯, 𝒮, ℐ, M)
    Ψ̂′ = zeros(2ℐ + 1, M + 1)
    Ψ̂″ = zeros(2ℐ + 1, M + 1)
    Ψ̃′ = complex(zeros(2ℐ + 1, M + 1))
    Ψ̃″ = complex(zeros(2ℐ + 1, M + 1))
    for m in 0:M
        Ψ̂′[:, m + 1] = iseven(m) ? κ .^ m : -(κ .^ m) .* 𝒯
        Ψ̂″[:, m + 1] = iseven(m) ? zero(κ) : κ .^ (m + 1) .* 𝒮
        Ψ̃′[:, m + 1] = (iseven(m) ? -(κ .^ m) .* 𝒯 : κ .^ m) * im
        Ψ̃″[:, m + 1] = (iseven(m) ? κ .^ (m + 1) .* 𝒮 : zero(κ)) * im
    end
    return Ψ̂′, Ψ̂″, Ψ̃′, Ψ̃″
end
