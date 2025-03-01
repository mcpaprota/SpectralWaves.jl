# Example 1 - propagation of a linear regular wave

using SpectralWaves
using CairoMakie

# Define fluid domain parameters
d = 1.0 # water depth (m)
ℓ = L # fluid domain length (m)

# Define wave parameters
H = 0.1 # wave height (m)
L = 2.0 # wavelength (m)
k = 2π / L # wave number (rad/m)
ω = sqrt(g * k * tanh(k * d)) # angular wave frequency (rad/s)
T = 2π / ω # wave period (s)

# Define numerical model parameters
ℐ = 1 # number of harmonics
nΔt = 200 # number of time steps per wave period
Δt = T / nΔt # time step (s)
nT = 1 # number of periods
t₀ = 0.0 # initial time (s)
τ = nT * T # total simulation time (s)
t = range(start = t₀, stop = τ, step = Δt) # time range

# Initialize wave problem
p = Problem(ℓ, d, ℐ, t)

# Define initial conditions
linear_regular_wave!(p, H, ω)

# Solve wave problem
solve_problem!(p)

# Plot results
set_theme!(theme_latexfonts())
fig = Figure(size = (400, 300))
ax = Axis(fig[1, 1], xlabel = L"t/T", ylabel = L"4η̂/H", xticks = 0:0.1:N, yticks = -1:0.5:1)
lines!(ax, t / T, 4 * abs.(p.η̂[1, p.O:end]) / H, label = L"|η̂|")
lines!(ax, t / T, 4 * real.(p.η̂[1, p.O:end]) / H, label = L"Re(η̂)")
lines!(ax, t / T, 4 * imag.(p.η̂[1, p.O:end]) / H, label = L"Im(η̂)")
axislegend(ax, position = :lb)
display(fig)
