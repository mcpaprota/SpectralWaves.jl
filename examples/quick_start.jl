using SpectralWaves
using CairoMakie # plotting package

d = 1.0 # water depth (m)
ℓ = 10.0 # fluid domain length (m)

ℐ = 40 # number of harmonics
Δt = 0.01 # time step (s)
t₀ = 0.0 # initial time (s)
τ = 2.0 # total simulation time (s)
t = range(start = t₀, stop = τ, step = Δt) # time range

p₀ = Problem(ℓ, d, ℐ, t)
p₁ = Problem(ℓ, d, ℐ, t; M_b=40)

h = 0.4d # surface bump height (m)
λ = 0.1ℓ # surface bump length (m)
surface_bump!(p₀, h, λ)
surface_bump!(p₁, h, λ)
h₁ = 0.9d # bottom bump height (m)
λ₁ = 0.5ℓ # bottom bump length (m)
bottom_bump!(p₁, h₁, λ₁)

solve_problem!(p₀)
solve_problem!(p₁)

η₀(x, n) = water_surface(p₀, x, n)
η₁(x, n) = water_surface(p₁, x, n)
β(x) = bottom_surface(p₁, x)

x = range(start = - ℓ / 2, stop = ℓ / 2, length = 1001) # spatial range
o₀ = Observable(η₀.(x, firstindex(t))) # set free-surface observable for p₀
o₁ = Observable(η₁.(x, firstindex(t))) # set free-surface observable for p₁
title = Observable(L"t = %$(round(t[1], digits=1))\,\mathrm{s}") # set string observable
set_theme!(theme_latexfonts()) # set latex fonts
fig = Figure(size = (700, 300)) # initialize a figure

# left plot p₀
ax0 = Axis(fig[1, 1],
        xlabel = L"$x$ (m)",
        ylabel = L"$z$ (m)") # define axis with labels
band!(ax0, x, o₀, -d,
        color=:azure) # plot water bulk
lines!(ax0, x, o₀,
        color=:black,
        linewidth = 1) # plot free surface line
band!(ax0, x, -1.1d, - d,
        color=:wheat) # plot bottom bulk
hlines!(ax0, -d,
        color=:black,
        linewidth = 0.7) # plot bottom line
limits!(ax0, x[1], x[end], -1.1d, d) # set limits

# right plot p₁
ax1 = Axis(fig[1, 2],
        xlabel = L"$x$ (m)") # define axis with labels
band!(ax1, x, o₁, β.(x) .-d,
        color=:azure) # plot water bulk
lines!(ax1, x, o₁,
        color=:black,
        linewidth = 1) # plot free surface line
band!(ax1, x, β.(x) .- d, - 1.1d,
        color=:wheat) # plot bottom bulk
lines!(ax1, x, β.(x) .- d,
        color=:black,
        linewidth = 1) # plot bottom line
limits!(ax1, x[1], x[end], -1.1d, d) # set limits
Label(fig[0, :], text=title)

# animate free surface
record(fig, "animation.mp4", 1:lastindex(t);
        framerate = 30) do n
    o₀[] = η₀.(x, n)
    o₁[] = η₁.(x, n)
    title[] = L"t = %$(round(t[n], digits=1))\,\mathrm{s}"
end
