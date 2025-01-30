# Example 3 - propagation of regular waves in a wave flume

using SpectralWaves
using GLMakie

# Define fluid domain and wave parameters
d = 1.0 # water depth (m)
H = 0.02 # wave height (m)
L = 5.0 # wavelength (m)
ℓ = 100.0 # fluid domain length (m)

# Define numerical model parameters
M_s = 0 # FSBC Taylor series order (linear wave)
M_b = 30 # BBC Taylor series order (horizontal bottom)
ℐ = 100 # number of harmonics
nΔt = 200 # number of time steps per wave period
nT = 20 # number of wave periods
nT₀ = 5 # number of ramped wave periods
N = nΔt * nT # number of time steps

# Initialize wave problem
p = Problem(ℓ, d, ℐ, N)

# Define wavemaker motion
T, Δt, t = linear_wavemaker!(p, H, L, nΔt, nT, nT₀)

# Define bathymetry - slope
h = 0.8 # slope height (m)
p.β̂[:, 1] = @. -4h / 3 * sinc(p.κ * ℓ / 3π)^2
p.β̂[ℐ+1, 1] = 2h / 3

# Solve wave problem
solve_problem!(p, M_s, M_b, Δt)

# Define free-surface elevation
η₁(x, n) = inverse_fourier_transform(p.η̂[:, n], p.κ, x)
β₁(x) = inverse_fourier_transform(p.β̂[:, 1], p.κ, x)
x = range(0, ℓ / 2, length = 1000)
η(n) = η₁.(x, n)
β = β₁.(x)
η₀  = Observable(η(p.O))

# animate free-surface elevation
fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1], xlabel = "x (m)", ylabel = "η (m)")
lines!(ax, x, η₀, color = :blue, linewidth = 2)
lines!(ax, x, β .- d, color = :yellow, linewidth = 2)
limits!(ax, 0, ℓ / 2, -d, 2H)
display(fig)

for n in p.O:10:N
    η₀[] = η(n)
    sleep(0.001)
end
