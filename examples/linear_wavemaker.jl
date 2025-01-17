# Example 2 - propagation of a regular waves in a wave flume

using SpectralWaves
using CairoMakie

# Define fluid domain and wave parameters
d = 1.0 # water depth (m)
H = 0.1 # wave height (m)
L = 2.0 # wavelength (m)
ℓ = 10.0 # fluid domain length (m)

# Define numerical model parameters
M_s = 0 # FSBC Taylor series order (linear wave)
M_b = 0 # BBC Taylor series order (horizontal bottom)
ℐ = 25 # number of harmonics
nΔt = 200 # number of time steps per wave period
nT = 10 # number of wave periods
nT₀ = 1 # number of ramped wave periods
N = nΔt * nT # number of time steps

# Initialize wave problem
κ, η̂, η̇, β̂, β̇, ϕ̂, ϕ̇, ψ̂, ψ̇, p̂, χ, ξ, ζ, 𝒯, 𝒮, O = init_problem(ℓ, d, ℐ, N)

# Define wavemaker motion
T, Δt, t = linear_wavemaker!(χ, ξ, ζ, H, L, d, nΔt, nT, nT₀, O)

# Solve wave problem
solve_problem!(η̂, η̇, ϕ̂, ϕ̇, ψ̂, ψ̇, β̂, β̇, p̂, κ, 𝒯, 𝒮, ℐ, M_s, M_b, Δt, O, N, χ, ξ, ζ, ℓ, d)

# Define free-surface elevation
η₁(x, n) = inverse_fourier_transform(η̂[:, n], κ, x)
x = range(0, ℓ, length = 100)
η(n) = η₁.(x, n)
η₀ = Observable(η(O))

# animate free-surface elevation
fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1], xlabel = "x (m)", ylabel = "η (m)")
lines!(ax, x, η₀[], color = :blue, linewidth = 2)
limits!(ax, 0, ℓ, -H, H)

record(fig, "linear_wavemaker.mp4", O:N, framerate = 10) do n
    η₀[] = η(n)
end
