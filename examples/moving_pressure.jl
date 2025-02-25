# Example 4 - waves excited by a moving bottom

using SpectralWaves
using GLMakie

# Define fluid domain parameters
d = 1.0 # water depth (m)
ℓ = 50.0 # fluid domain length (m)

# Define numerical model parameters
M_s = 0 # FSBC Taylor series order
M_b = 4 # BBC Taylor series order
ℐ = 100 # number of harmonics

# Define moving bottom parameters
h = 0.1 # bottom obstacle height (m)
μ = 0.2 # bottom obstacle relative length
Fr = 0.5 # bottom obstacle Froude number
λ = d / μ # bottom obstacle characteristic length (m)
T = λ / sqrt(g * d) # bottom obstacle characteristic period (s)
nΔt = 100 # number of time steps per wave period
nT = 10 # number of simulated obstacle periods
Δt = T / nΔt # time step (s)
t₀ = 0.0 # initial time (s)
τ = nT * T # total simulation time (s)
t = range(start = t₀, stop = τ, step = Δt) # time range

# Initialize wave problem
p = Problem(ℓ, d, ℐ, t)

# Bottom topography
moving_pressure_bump!(p, h, λ, Fr * sqrt(g * d))


# Solve wave problem
solve_problem!(p)

# Define free-surface elevation
η(x, n) = water_surface(p, x, n)
β(x) = bottom_surface(p, x)
pres(x, n) = surface_pressure(p, x, n)
x = range(- ℓ / 2, ℓ / 2, length = 1000)
η₀  = Observable(η.(x, p.O))
p₀  = Observable(pres.(x, p.O))

# animate free-surface elevation
fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1], xlabel = "x (m)", ylabel = "η (m)")
lines!(ax, x, η₀, color = :blue, linewidth = 2)
lines!(ax, x, p₀, color = :red, linewidth = 2)
lines!(ax, x, β.(x) .- d, color = :yellow, linewidth = 2)
limits!(ax, - ℓ / 2, ℓ / 2, -d, d)
display(fig)

for n in p.O:10:p.N
    η₀[] = η.(x, n)
    p₀[] = pres.(x, n)
    sleep(0.001)
end
