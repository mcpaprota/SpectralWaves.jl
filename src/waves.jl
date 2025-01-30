function linear_regular_wave!(p::Problem, H, L, nΔt)
    η̂, η̇, ϕ̂, ϕ̇ = p.η̂, p.η̇, p.ϕ̂, p.ϕ̇
    d, ℐ, O = p.d, p.ℐ, p.O
    a = H / 2 # wave amplitude (m)
    k = 2π / L # wave number (rad/m)
    ω = √(g * k * tanh(k * d)) # angular wave frequency (rad/s)
    T = 2π / ω # wave period (s)
    Δt = T / nΔt # time step (s)
    η̂[ℐ, O] = a / 2
    η̂[ℐ+2, O] = a / 2
    [η̇[ℐ, i] = im * ω * a / 2 * exp((i - O) * im * ω * Δt) for i in 1:O]
    [η̇[ℐ+2, i] = -im * ω * a / 2 * exp(-(i - O) * im * ω * Δt) for i in 1:O]
    ϕ̂[ℐ, O] = im * a * g / ω / 2
    ϕ̂[ℐ+2, O] = -im * a * g / ω / 2
    [ϕ̇[ℐ, i] = -g * a / 2 * exp((i - O) * im * ω * Δt) for i in 1:O]
    [ϕ̇[ℐ+2, i] = -g * a / 2 * exp(-(i - O) * im * ω * Δt) for i in 1:O]
    return T, Δt
end

function linear_wavemaker!(χ, ξ, ζ, H, L, d, nΔt, nT, nT₀, O)
    a = H / 2 # wave amplitude (m)
    k = 2π / L # wave number (rad/m)
    ω = √(g * k * tanh(k * d)) # angular wave frequency (rad/s)
    T = 2π / ω # wave period (s)
    Δt = T / nΔt # time step (s)
    s = a / (cosh(2k * d) - 1) * (sinh(2k * d) + 2k * d) / 2 # paddle motion amplitude (m)
    t = range(0, nT * T, step = Δt) # time vector (s)
    h(t) = s * sin(ω * t)
    h′(t) = s * ω * cos(ω * t)
    h″(t) = - s * ω^2 * sin(ω * t)
    r(t) = t < nT₀ * T ? 1 / 2 - cos(ω / nT₀ * t / 2) / 2 : 1
    r′(t) = t < nT₀ * T ? ω / 2 / nT₀ * sin(ω / nT₀ * t / 2) / 2 : 0
    r″(t) = t < nT₀ * T ? ω^2 / 4 / nT₀^2 * cos(ω / nT₀ * t / 2) / 2 : 0
    χ[O:end] = @. h(t) * r(t)
    ξ[O:end] = @. h′(t) * r(t) + h(t) * r′(t)
    ζ[O:end] = @. h″(t) * r(t) + 2 * h′(t) * r′(t) + h(t) * r″(t)
    return T, Δt, t
end
