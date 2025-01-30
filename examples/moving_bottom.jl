# Example 4 - waves excited by a moving bottom

using SpectralWaves
using GLMakie

# Define fluid domain parameters
d = 1.0 # water depth (m)
ℓ = 50.0 # fluid domain length (m)

# Define numerical model parameters
M_s = 2 # FSBC Taylor series order
M_b = 4 # BBC Taylor series order
ℐ = 100 # number of harmonics

# Define moving bottom parameters
h = 0.1 # bottom obstacle height (m)
μ = 0.2 # bottom obstacle relative length
Fr = 0.5 # bottom obstacle Froude number
λ = d / μ # bottom obstacle characteristic length (m)
T = λ / sqrt(g * d) # bottom obstacle characteristic period (s)
nΔt = 100 # number of time steps per wave period
nT = 1 # number of simulated obstacle periods
Δt = T / nΔt # time step (s)
t₀ = 0.0 # initial time (s)
τ = nT * T # total simulation time (s)
t = range(start = t₀, stop = τ, step = Δt) # time range

# Initialize wave problem
p = Problem(ℓ, d, ℐ, t; static_bottom = false, M_s = M_s, M_b = M_b)

# Bottom topography
for n in p.O:p.N+p.O
    p.β̂[:, n] = @. h * λ * sqrt(2π) / 4ℓ * exp(-λ^2 * p.κ^2 / 32) * exp(-im * (n - p.O) * Δt * Fr * sqrt(g * d) * p.κ)
    p.β̇[:, n] = @. -im * Fr * sqrt(g * d) * p.κ * p.β̂[:, n]
end

# Solve wave problem
solve_problem!(p)

# Define free-surface elevation
η₁(x, n) = inverse_fourier_transform(p.η̂[:, n], p.κ, x)
β₁(x, n) = inverse_fourier_transform(p.β̂[:, n], p.κ, x)
x = range(- ℓ / 2, ℓ / 2, length = 1000)
η(n) = η₁.(x, n)
β(n) = β₁.(x, n)
η₀  = Observable(η(p.O))
β₀  = Observable(β(p.O) .- d)

# animate free-surface elevation
fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1], xlabel = "x (m)", ylabel = "η (m)")
lines!(ax, x, η₀, color = :blue, linewidth = 2)
lines!(ax, x, β₀, color = :yellow, linewidth = 2)
limits!(ax, - ℓ / 2, ℓ / 2, -d, d)
display(fig)

for n in p.O:10:p.N
    η₀[] = η(n)
    β₀[] = β(n) .- d
    sleep(0.001)
end
