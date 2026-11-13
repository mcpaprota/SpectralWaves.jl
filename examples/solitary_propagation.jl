# Example - solitary wave propagation

using SpectralWaves
using CairoMakie # plotting package

# Define fluid domain parameters
H = 0.2 # wave height (m)
d = 1.0 # water depth (m)
ℓ = 100.0 # fluid domain length (m)

# Define numerical model parameters
ℐ = 60 # number of harmonics
Δt = 0.01 # time step (s)
t₀ = 0.0 # initial time (s)
τ = 10.0 # total simulation time (s)
t = range(start = t₀, stop = τ, step = Δt) # time range

# Initialize wave problem and define wavemaker motion
p = Problem(ℓ, d, ℐ, t; M_s = 4)
solitary_wavemaker!(p, H)

# Solve wave problem
solve_problem!(p; msg_flag = false)

# Plot the solitary wave profile at the final time
η(x, n) = water_surface(p, x, n)
β(x) = bottom_surface(p, x)
n = length(t)
x = range(start = 0, stop = ℓ / 2, length = 501) # spatial range

set_theme!(theme_latexfonts())
fig = Figure(size = (700, 300))
ax = Axis(fig[1, 1],
        xlabel = L"$x$ (m)",
        ylabel = L"$z$ (m)")
band!(ax, x, η.(x, n), β.(x) .- d,
        color = :azure)
lines!(ax, x, η.(x, n),
        color = :black,
        linewidth = 1)
band!(ax, x, β.(x) .- d, -1.1d,
        color = :wheat)
lines!(ax, x, β.(x) .- d,
        color = :black,
        linewidth = 1)
limits!(ax, x[1], x[end], -1.1d, 2H)
display(fig)
