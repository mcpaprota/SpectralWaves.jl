using SpectralWaves
using CairoMakie # plotting package

L = 5.0 # wavelength (m)
H = 0.05 # wave height (m)
d = 1.0 # water depth (m)
ℓ = 200.0 # fluid domain length (m)

k = 2π / L # wave number (rad/m)
ω = sqrt(g * k * tanh(k * d)) # angular wave frequency (rad/s)
T = 2π / ω # wave period (s)

ℐ = 600 # number of harmonics
nT = 30 # number of simulated wave periods
nT₀ = 3 # number of ramped wave periods
nΔt = 200 # number of time steps per wave period
Δt = T / nΔt # time step (s)
t₀ = 0.0 # initial time (s)
τ = nT * T # total simulation time (s)
t = range(start = t₀, stop = τ, step = Δt) # time range

p = Problem(ℓ, d, ℐ, t; M_b=60, M_s=2)

linear_wavemaker!(p, H, T, L, nT₀)

h = 0.9d
bottom_slope!(p, h)

solve_problem!(p)

η(x, n) = water_surface(p, x, n)
β(x) = bottom_surface(p, x)

x = range(start = 0, stop = ℓ / 3, length = 501) # spatial range
η₀ = Observable(η.(x, firstindex(t))) # set free-surface observable for p
set_theme!(theme_latexfonts()) # set latex fonts
fig = Figure(size = (700, 300)) # initialize a figure
ax = Axis(fig[1, 1],
        xlabel = L"$x$ (m)",
        ylabel = L"$z$ (m)") # define axis with labels and title
band!(ax, x, η₀, β.(x) .-d,
        color=:azure) # plot water bulk
lines!(ax, x, η₀,
        color=:black,
        linewidth = 1) # plot free surface line
band!(ax, x, β.(x) .- d, - 1.1d,
        color=:wheat) # plot bottom bulk
lines!(ax, x, β.(x) .- d,
        color=:black,
        linewidth = 1) # plot bottom line
limits!(ax, x[1], x[end], -1.1d, 2H) # set limits

# animate free surface
record(fig, "nonlinear_shoaling.mp4", lastindex(t)-nΔt+1:lastindex(t);
        framerate = nΔt) do n
    η₀[] = η.(x, n)
end
