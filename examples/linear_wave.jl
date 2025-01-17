# Example 1 - propagation of a linear regular wave

using SpectralWaves
using CairoMakie

# Define fluid domain and wave parameters
d = 1.0 # water depth (m)
H = 0.1 # wave height (m)
L = 2.0 # wavelength (m)
ℓ = L # fluid domain length (m)

# Define numerical model parameters
M_s = 0 # FSBC Taylor series order (linear wave)
M_b = 0 # BBC Taylor series order (horizontal bottom)
ℐ = 1 # number of harmonics
nΔt = 200 # number of time steps per wave period
nT = 1 # number of periods
N = nΔt * nT # number of time steps

# Initialize wave problem
κ, η̂, η̇, β̂, β̇, ϕ̂, ϕ̇, ψ̂, ψ̇, p̂, χ, ξ, ζ, 𝒯, 𝒮, O = init_problem(ℓ, d, ℐ, N; O=4)

# Define initial conditions
T, Δt = linear_regular_wave!(η̂, η̇, ϕ̂, ϕ̇, H, L, d, ℐ, nΔt, O)

# Solve wave problem
solve_problem!(η̂, η̇, ϕ̂, ϕ̇, ψ̂, ψ̇, β̂, β̇, p̂, κ, 𝒯, 𝒮, ℐ, M_s, M_b, Δt, O, N, χ, ξ, ζ, ℓ, d)

# Plot results
t = range(0, N*Δt, step = Δt)
set_theme!(theme_latexfonts())
fig = Figure(size = (400, 300))
ax = Axis(fig[1, 1], xlabel = L"t/T", ylabel = L"4η̂/H", xticks = 0:0.1:N, yticks = -1:0.5:1)
lines!(ax, t / T, 4 * abs.(η̂[1, O:end]) / H, label = L"|η̂|")
lines!(ax, t / T, 4 * real.(η̂[1, O:end]) / H, label = L"Re(η̂)")
lines!(ax, t / T, 4 * imag.(η̂[1, O:end]) / H, label = L"Im(η̂)")
axislegend(ax, position = :lb)
display(fig)
