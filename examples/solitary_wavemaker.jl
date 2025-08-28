# Example: propagation of a solitary wave using a wavemaker

using SpectralWaves
using CairoMakie

# Define fluid domain
d = 0.5 # water depth (m)
ℓ = 100.0 # fluid domain length (m)

# Define solitary wave parameters
H = 0.05 # solitary wave height (m)

# Define numerical model parameters
ℐ = 100 # number of harmonics
nΔt = 1000 # number of time steps per wave period
τ = 15.0 # total simulation time (s), adjust as needed
Δt = τ / nΔt # time step (s)
t₀ = 0.0 # initial time (s)
t = range(start = t₀, stop = τ, step = Δt) # time range
x = range(0, ℓ/ 2, length = 1001) # spatial range

# Initialize wave problem
p = Problem(ℓ, d, ℐ, t; M_s = 2)

# Define wavemaker motion
solitary_wavemaker!(p, H)

# Solve wave problem
solve_problem!(p)

# Plot the last position of the free surface in one figure
η(x) = water_surface(p, x, lastindex(t))
β(x) = bottom_surface(p, x)

set_theme!(theme_latexfonts())
fig = Figure()
ax = Axis(fig[1, 1], xlabel = L"$x$ (m)", ylabel = L"$z$ (m)")

band!(ax, x, η.(x), β.(x) .- d, color=:azure) # water bulk
band!(ax, x, β.(x) .- d, -1.1d, color=:wheat) # bottom bulk
lines!(ax, x, η.(x), color=:black, linewidth = 0.7) # free surface
lines!(ax, x, β.(x) .- d, color=:black, linewidth = 0.7) # bottom surface
limits!(ax, x[1], x[end], -1.1d, 4H)

# Annotate maximum free-surface elevation
η_vals = η.(x)
imax = argmax(η_vals)
xmax = x[imax]
ηmax = η_vals[imax]

scatter!(ax, [xmax], [ηmax], color=:red, markersize=12, label="max η")
text!(ax, "Maximum = $(round(ηmax, digits=4))", position = (xmax, ηmax + 0.02), align = (:center, :bottom), color=:red)

display(fig)
