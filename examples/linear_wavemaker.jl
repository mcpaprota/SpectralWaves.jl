# Example 2 - propagation of regular waves in a wave flume

using SpectralWaves
using GLMakie

# Define fluid domain
d = 1.0 # water depth (m)
ℓ = 100.0 # fluid domain length (m)

# Define wave parameters
H = 0.2 # wave height (m)
L = 2 # wave length (m)
k = 2π / L # wave number (rad/m)
ω = sqrt(g * k * tanh(k * d)) # angular wave frequency (rad/s)
T = 2π / ω # wave period (s)

# Define numerical model parameters
ℐ = 100 # number of harmonics
nT = 40 # number of simulated wave periods
nT₀ = 5 # number of ramped wave periods
nΔt = 100 # number of time steps per wave period
Δt = T / nΔt # time step (s)
t₀ = 0.0 # initial time (s)
τ = nT * T # total simulation time (s)
t = range(start = t₀, stop = τ, step = Δt) # time range
x = range(0, ℓ, length = 1001) # spatial range
βₙ = @. - 0.3 * cos(1.5 * x) * (cos(x) - 1) # bottom profile

# Initialize wave problem
p = Problem(ℓ, d, ℐ, t)
p2 = Problem(ℓ, d, ℐ, t; M_b = 30)

# Define wavemaker motion
linear_wavemaker!(p, H, T, L, nT₀)
linear_wavemaker!(p2, H, T, L, nT₀)

# Define bottom profile
#bottom_vector!(p2, x, βₙ)
bottom_slope!(p2, 0.9d)

# Solve wave problem
solve_problem!(p)
solve_problem!(p2)

# plot results
η(x) = water_surface(p, x, lastindex(t))
β(x) = bottom_surface(p, x)
set_theme!(theme_latexfonts())
fig = Figure()
ax = Axis(fig[1, 1], xlabel = L"$x$ (m)", ylabel = L"$z$ (m)")
band!(ax, x, η.(x), β.(x) .- d, color=:azure) # water bulk
band!(ax, x, β.(x) .- d, -1.1d, color=:wheat) # bottom bulk
lines!(ax, x, η.(x), color=:black, linewidth = 0.7) # free surface
lines!(ax, x, β.(x) .- d, color=:black, linewidth = 0.7) # bottom surface
limits!(ax, x[1], x[end], -1.1d, d)
display(fig)
