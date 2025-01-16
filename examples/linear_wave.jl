# Example 1 - propagation of a linear regular wave

using SpectralWaves
using CairoMakie

# Define fluid domain and wave parameters
d = 1.0 # water depth (m)
H = 0.1 # wave height (m)
L = 2.0 # wavelength (m)
ℓ = L # fluid domain length (m)
a = H / 2 # wave amplitude (m)
k = 2π / L # wave number (rad/m)
ω = √(g * k * tanh(k * d)) # angular wave frequency (rad/s)
T = 2π / ω # wave period (s)

# Define numerical model parameters
M_s = 0 # FSBC Taylor series order (linear wave)
M_b = 0 # BBC Taylor series order (horizontal bottom)
ℐ = 10 # number of harmonics
nΔt = 100 # number of time steps per wave period
Δt = T / nΔt # time step (s)
nT = 1 # number of periods
N = nΔt * nT # number of time steps

# Initialize wave problem
κ, η̂, η̇, β̂, β̇, ϕ̂, ϕ̇, ψ̂, ψ̇, p̂, χ, ξ, ζ, 𝒯, 𝒮, O = init_problem(ℓ, d, ℐ, N)

# Define initial conditions
η̂[ℐ:ℐ+2, O] = [a / 2, 0, a / 2]
[η̇[ℐ, i] = im * ω * a / 2 * exp((i - O) * im * ω * Δt) for i in 1:O]
[η̇[ℐ+2, i] = - im * ω * a / 2 * exp(-(i - O) * im * ω * Δt) for i in 1:O]
ϕ̂[ℐ:ℐ+2, O] = [im * a * g / ω / 2, 0, -im * a * g / ω / 2]
[ϕ̇[ℐ, i] = -g * a / 2 * exp((i - O) * im * ω * Δt) for i in 1:O]
[ϕ̇[ℐ+2, i] = -g * a / 2 * exp(-(i - O) * im * ω * Δt) for i in 1:O]

# Solve wave problem
solve_problem!(η̂, η̇, ϕ̂, ϕ̇, ψ̂, ψ̇, β̂, β̇, p̂, κ, 𝒯, 𝒮, ℐ, M_s, M_b, Δt, O, N, χ, ξ, ζ, ℓ, d)
