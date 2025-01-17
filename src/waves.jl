function linear_regular_wave!(η̂, η̇, ϕ̂, ϕ̇, H, L, d, ℐ, nΔt, O)
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
