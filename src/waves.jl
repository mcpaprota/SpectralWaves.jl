"""
    linear_regular_wave!(p::Problem, H, ω)

Initialize `p` with the free-surface elevation, velocity, and potential coefficients
for a regular linear wave of height `H` and angular frequency `ω` (rad/s).
"""
function linear_regular_wave!(p::Problem, H, ω)
    η̂, η̇, ϕ̂, ϕ̇ = p.η̂, p.η̇, p.ϕ̂, p.ϕ̇
    ℐ, Δt, O = p.ℐ, p.Δt, p.O
    a = H / 2 # wave amplitude (m)
    η̂[ℐ, O] = a / 2
    η̂[ℐ+2, O] = a / 2
    [η̇[ℐ, i] = im * ω * a / 2 * exp((i - O) * im * ω * Δt) for i in 1:O]
    [η̇[ℐ+2, i] = -im * ω * a / 2 * exp(-(i - O) * im * ω * Δt) for i in 1:O]
    ϕ̂[ℐ, O] = im * a * g / ω / 2
    ϕ̂[ℐ+2, O] = -im * a * g / ω / 2
    [ϕ̇[ℐ, i] = -g * a / 2 * exp((i - O) * im * ω * Δt) for i in 1:O]
    [ϕ̇[ℐ+2, i] = -g * a / 2 * exp(-(i - O) * im * ω * Δt) for i in 1:O]
    return nothing
end

"""
    linear_wavemaker!(p::Problem, H, T, L, nT₀)

Calculate wavemaker paddle displacement `χ`, velocity `ξ`, and acceleration `ζ`
for linear waves of height `H`, period `T`, and wavelength `L`. The paddle motion
is ramped over `nT₀` periods.
"""
function linear_wavemaker!(p::Problem, H, T, L, nT₀)
    χ, ξ, ζ = p.χ, p.ξ, p.ζ
    d, O = p.d, p.O
    t = p.t
    ω = 2π / T # angular wave frequency (rad/s)
    t₀ = nT₀ * T # ramped period (s)
    a = H / 2 # wave amplitude (m)
    k = 2π / L # wave number (rad/m)
    s = a / (cosh(2k * d) - 1) * (sinh(2k * d) + 2k * d) / 2 # paddle motion amplitude (m)
    h(t) = s * sin(ω * t)
    h′(t) = s * ω * cos(ω * t)
    h″(t) = - s * ω^2 * sin(ω * t)
    r(t) = t < t₀ ? 1 / 2 - cos(2π / t₀ * t / 2) / 2 : 1
    r′(t) = t < t₀ ? π / t₀ * sin(2π / t₀ * t / 2) / 2 : 0
    r″(t) = t < t₀ ? (π / t₀)^2 * cos(2π / t₀ * t / 2) / 4 : 0
    χ[O:end] = @. h(t) * r(t)
    ξ[O:end] = @. h′(t) * r(t) + h(t) * r′(t)
    ζ[O:end] = @. h″(t) * r(t) + 2 * h′(t) * r′(t) + h(t) * r″(t)
    return nothing
end

"""
    solitary_wavemaker!(p::Problem, H; centered = false)

Calculate wavemaker paddle displacement `χ`, velocity `ξ`, and acceleration `ζ`
for a solitary wave of height `H`.
"""
function solitary_wavemaker!(p::Problem, H)
    χ, ξ, ζ = p.χ, p.ξ, p.ζ
    d, O, t = p.d, p.O, p.t

    # Goring (1978), Goring and Reichlen (1980) solitary wave parameters
    c = sqrt(g * (d + H)) # solitary wave celerity
    beta = sqrt(3 * H / d^3) # wavenumber parameter
    alpha = atanh(0.999) # cutoff factor
    tau = 4 / beta / c * (alpha + H / d) # motion duration
    S = 4 * sqrt(H * d / 3) # total stroke length
    S₀ = S / 2 # stroke amplitude

    θ(ti) = 2 * alpha * (ti / tau - 0.5)

    # paddle motion profiles
    χ_func(ti) = ti < tau ? S₀ * (tanh(θ(ti)) - tanh(-alpha)) : S₀ * (tanh(alpha) - tanh(-alpha))
    ξ_func(ti) = ti < tau ? S₀ * (2 * alpha / tau) * (1 - tanh(θ(ti))^2) : 0.0
    ζ_func(ti) = ti < tau ? -2 * S₀ * (2 * alpha / tau)^2 * tanh(θ(ti)) * (1 - tanh(θ(ti))^2) : 0.0

    # assign paddle motion to problem
    χ[O:end] = χ_func.(t)
    ξ[O:end] = ξ_func.(t)
    ζ[O:end] = ζ_func.(t)

    return nothing
end

"""
    surface_bump!(p::Problem, h, λ, x₀ = 0)

Calculate `η̂` coefficients of a problem `p` for a surface bump of height `h` and characteristic length `λ` at position `x₀`.

"""
function surface_bump!(p::Problem, h, λ, x₀ = 0)
    η̂, κ, ℓ, O = p.η̂, p.κ, p.ℓ, p.O
    for n = 1:O
        η̂[:, n] = @. h * λ * √(2π) / 4ℓ * exp(- λ^2 * κ^2 / 32) * exp(-im * κ * x₀)
    end
    return nothing
end

"""
    moving_pressure_bump!(p::Problem, h, λ, u, x₀ = 0)

Calculate `p̂` coefficients of a problem `p` for a moving pressure bump of height `h`, characteristic length `λ`, and velocity `u` at initial position `x₀`.

"""
function moving_pressure_bump!(p::Problem, h, λ, u, x₀ = 0)
    p̂, κ, ℓ, t, O, N = p.p̂, p.κ, p.ℓ, p.t, p.O, p.N
    for n in O:N+O-1
        p̂[:, n] = @. h * λ * √(2π) / 4ℓ * exp(- λ^2 * κ^2 / 32) * exp(-im * κ * x₀) * exp(-im * κ * u * t[n-O+1])
    end
    return nothing
end
