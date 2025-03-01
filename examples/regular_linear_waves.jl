using SpectralWaves
using CairoMakie # plotting package

L = 2.0 # wavelength (m)
H = 0.1 # wave height (m)
d = 1.0 # water depth (m)
ℓ = L # fluid domain length (m) - one wave

k = 2π / L # wave number (rad/m)
ω = sqrt(g * k * tanh(k * d)) # angular wave frequency (rad/s)
T = 2π / ω # wave period (s)

ℐ = 1 # number of harmonics
nΔt = 200 # number of time increments per wave period
Δt = T / nΔt # time step (s)
nT = 1 # number of periods
t₀ = 0.0 # initial time (s)
τ = nT * T # total simulation time (s)
t = range(start = t₀, stop = τ, step = Δt) # time range

p = Problem(ℓ, d, ℐ, t)

linear_regular_wave!(p, H, ω)

solve_problem!(p)

η(x, n) = water_surface(p, x, n)
u(x, z, n) = (z < η(x, n)) * water_velocity(p, x, z, n, :x)
w(x, z, n) = (z < η(x, n)) * water_velocity(p, x, z, n, :z)
v(x, z, n) = sqrt.(u(x, z, n)^2 + w(x, z, n)^2)

x = range(start = 0, stop = ℓ, length = 21)
z = range(start = -d, stop = H, length = 11)

n = 101
set_theme!(theme_latexfonts()) # set latex fonts
fig = Figure(size = (700, 300))
ax = Axis(fig[1, 1], xlabel = L"$x$ (m)", ylabel = L"$\eta$ (m)")
band!(ax, x, η.(x, n), -d,
        color=:azure) # plot water bulk
lines!(ax, x, η.(x, n),
        color=:black,
        linewidth = 1) # plot free surface
band!(ax, x, -1.1d, - d,
        color=:wheat) # plot bottom bulk
hlines!(ax, -d,
        color=:black,
        linewidth = 0.7) # plot bottom line
arrows!(ax, x, z, u.(x, z', n), w.(x, z', n);
    lengthscale = 0.3,
    arrowsize = 10 * vec(v.(x, z', n)/maximum(v.(x, z', n)))) # plot velocity vectors
limits!(ax, 0, L, -1.1d, H)
fig
