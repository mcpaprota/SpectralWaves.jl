# Example 4 - waves excited by a moving bottom

using SpectralWaves
using GLMakie

# Define fluid domain parameters
d = 1.0 # water depth (m)
ℓ = 100.0 # fluid domain length (m)

# Define numerical model parameters
M_s = 2 # FSBC Taylor series order
M_b = 4 # BBC Taylor series order
ℐ = 100 # number of harmonics

# Define moving bottom parameters
h = 0.1 # bottom obstacle height (m)
μ = 0.2 # bottom obstacle relative length
Fr = 1.0 # bottom obstacle Froude number
λ = d / μ # bottom obstacle characteristic length (m)
T = λ / sqrt(g * d) # bottom obstacle characteristic period (s)
nΔt = 200 # number of time steps per wave period
nT = 5 # number of simulated obstacle periods
Δt = T / nΔt # time step (s)
N = nΔt * nT # number of time steps

# Initialize wave problem
κ, η̂, η̇, β̂, β̇, ϕ̂, ϕ̇, ψ̂, ψ̇, p̂, χ, ξ, ζ, 𝒯, 𝒮, O = init_problem(ℓ, d, ℐ, N)

# Bottom topography
for n in O:N+O
    β̂[:, n] = @. h * λ * sqrt(2π) / 4ℓ * exp(-λ^2 * κ^2 / 32) * exp(-im * (n - O) * Δt * Fr * sqrt(g * d) * κ)
    β̇[:, n] = @. -im * Fr * sqrt(g * d) * κ * β̂[:, n]
end

# Solve wave problem
solve_problem!(η̂, η̇, ϕ̂, ϕ̇, ψ̂, ψ̇, β̂, β̇, p̂, κ, 𝒯, 𝒮, ℐ, M_s, M_b, Δt, O, N, χ, ξ, ζ, ℓ, d; static_bottom=false)

# Define free-surface elevation
η₁(x, n) = inverse_fourier_transform(η̂[:, n], κ, x)
β₁(x, n) = inverse_fourier_transform(β̂[:, n], κ, x)
x = range(- ℓ / 2, ℓ / 2, length = 1000)
η(n) = η₁.(x, n)
β(n) = β₁.(x, n)
η₀  = Observable(η(O))
β₀  = Observable(β(O) .- d)

# animate free-surface elevation
fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1], xlabel = "x (m)", ylabel = "η (m)")
lines!(ax, x, η₀, color = :blue, linewidth = 2)
lines!(ax, x, β₀, color = :yellow, linewidth = 2)
limits!(ax, - ℓ / 2, ℓ / 2, -d, d)
display(fig)

for n in O:10:N
    η₀[] = η(n)
    β₀[] = β(n) .- d
    sleep(0.001)
end
