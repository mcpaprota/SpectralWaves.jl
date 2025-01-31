# Example 3 - propagation of regular waves in a wave flume

using SpectralWaves
using GLMakie

# Define fluid domain
d = 1.0 # water depth (m)
ℓ = 100.0 # fluid domain length (m)

# Define wave parameters
H = 0.02 # wave height (m)
L = 5.0 # wavelength (m)
k = 2π / L # wave number (rad/m)
ω = sqrt(g * k * tanh(k * d)) # angular wave frequency (rad/s)
T = 2π / ω # wave period (s)

# Define numerical model parameters
M_b = 30 # BBC Taylor series order (horizontal bottom)
ℐ = 200 # number of harmonics
nΔt = 200 # number of time steps per wave period
Δt = T / nΔt # time step (s)
nT = 30 # number of wave periods
nT₀ = 5 # number of ramped wave periods
t₀ = 0.0 # initial time (s)
τ = nT * T # total simulation time (s)
t = range(start = t₀, stop = τ, step = Δt) # time range

# Initialize wave problem
p = Problem(ℓ, d, ℐ, t; M_b=M_b)

# Define wavemaker motion
linear_wavemaker!(p, H, T, L, nT₀)

# Define bathymetry - slope
h = 0.97 # slope height (m)
p.β̂[:, 1] = @. -4h / 3 * sinc(p.κ * ℓ / 3π)^2
p.β̂[ℐ+1, 1] = 2h / 3

# Solve wave problem
solve_problem!(p)

# Define free-surface elevation
η(x, n) = water_surface(p, x, n)
β(x) = bottom_surface(p, x)
x = range(0, ℓ / 2, length = 1000)
η₀  = Observable(η.(x, p.O))

# animate free-surface elevation
fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1], xlabel = "x (m)", ylabel = "η (m)")
lines!(ax, x, η₀, color = :blue, linewidth = 2)
lines!(ax, x, β.(x) .- d, color = :yellow, linewidth = 2)
limits!(ax, 0, ℓ / 2, -d, 2H)
display(fig)

for n in p.O:10:p.N
    η₀[] = η.(x, n)
    sleep(0.001)
end
