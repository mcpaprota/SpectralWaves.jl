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
ℐ = 200 # number of harmonics
nT = 20 # number of simulated wave periods
nT₀ = 5 # number of ramped wave periods
nΔt = 200 # number of time steps per wave period
Δt = T / nΔt # time step (s)
t₀ = 0.0 # initial time (s)
τ = nT * T # total simulation time (s)
t = range(start = t₀, stop = τ, step = Δt) # time range

# Initialize wave problem
p = Problem(ℓ, d, ℐ, t)

# Define wavemaker motion
linear_wavemaker!(p, H, T, L, nT₀)

# Solve wave problem
solve_problem!(p)

# Define free-surface elevation
x = range(0, ℓ / 2, length = 500)
η(x, n) = water_surface(p, x, n)
η₀  = Observable(η.(x, p.O))

# animate free-surface elevation
fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1], xlabel = "x (m)", ylabel = "η (m)")
lines!(ax, x, η₀, color = :blue, linewidth = 2)
limits!(ax, 0, ℓ / 2, -H, H)
display(fig)

for n in p.O:10:p.N
    η₀[] = η.(x, n)
    sleep(0.001)
end
