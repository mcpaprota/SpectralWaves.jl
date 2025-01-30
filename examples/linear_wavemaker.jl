# Example 2 - propagation of regular waves in a wave flume

using SpectralWaves
using GLMakie

# Define fluid domain and wave parameters
d = 1.0 # water depth (m)
H = 0.2 # wave height (m)
L = 5.0 # wavelength (m)
ℓ = 200.0 # fluid domain length (m)

# Define numerical model parameters
M_s = 0 # FSBC Taylor series order (linear wave)
M_b = 0 # BBC Taylor series order (horizontal bottom)
ℐ = 200 # number of harmonics
nΔt = 200 # number of time steps per wave period
nT = 21 # number of wave periods
nT₀ = 5 # number of ramped wave periods
N = nΔt * nT # number of time steps

# Initialize wave problem
p = Problem(ℓ, d, ℐ, N)

# Define wavemaker motion
T, Δt, t = linear_wavemaker!(p, H, L, nΔt, nT, nT₀)

# Solve wave problem
solve_problem!(p, M_s, M_b, Δt)

# Define free-surface elevation
η₁(x, n) = inverse_fourier_transform(p.η̂[:, n], p.κ, x)
x = range(0, ℓ / 2, length = 500)
η(n) = η₁.(x, n)
η₀  = Observable(η(p.O))

# animate free-surface elevation
fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1], xlabel = "x (m)", ylabel = "η (m)")
lines!(ax, x, η₀, color = :blue, linewidth = 2)
limits!(ax, 0, ℓ / 2, -H, H)
display(fig)

for n in p.O:10:N
    η₀[] = η(n)
    sleep(0.001)
end
