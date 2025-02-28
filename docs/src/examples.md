# Examples

## Regular linear waves

We are modelling linear and regular waves of length `L` and height `H` propagating in water of constant depth `d`. We consider one wave along the length of the domain `ℓ`. 

```@example 1
using SpectralWaves
using CairoMakie # plotting package

L = 2.0 # wavelength (m)
H = 0.1 # wave height (m)
d = 1.0 # water depth (m)
ℓ = L # fluid domain length (m) - one wave
nothing # hide
```

We need a wave period `T`.

```@example 1
k = 2π / L # wave number (rad/m)
ω = sqrt(g * k * tanh(k * d)) # angular wave frequency (rad/s)
T = 2π / ω # wave period (s)
nothing # hide
```

We define a number of numerical model parameters.

```@example 1
ℐ = 1 # number of harmonics
nΔt = 200 # number of time increments per wave period
Δt = T / nΔt # time step (s)
nT = 1 # number of periods
t₀ = 0.0 # initial time (s)
τ = nT * T # total simulation time (s)
t = range(start = t₀, stop = τ, step = Δt) # time range
nothing # hide
```

We initialize wave problem using a struct [`Problem`](@ref).

```@example 1
p = Problem(ℓ, d, ℐ, t)
nothing # hide
```

Initial condition values of ``\hat{\eta}``, ``\dot{\eta}``, ``\hat{\phi}``, and ``\dot{\phi}`` are computed and inserted into vectors `η̂`, `η̇`, `ϕ̂`, `ϕ̇` using [`linear_regular_wave!`](@ref) in-place function.

```@example 1
linear_regular_wave!(p, H, ω)
nothing # hide
```

Now, we are ready to solve a problem. We use an in-place function [`solve_problem!`](@ref) which stores the values of solution coefficients in vectors `η̂`, `η̇`, `ϕ̂`, `ϕ̇`, `ψ̂`, `ψ̇`. In our case only `η̂` will be further processed.

```@example 1
solve_problem!(p)
nothing # hide
```

We define a function that calculates free-surface elevation at specified location `x` and time instant `n` and velocity components `u` and `w` at specified location (`x`, `z`) and time instant `n` using [`water_surface`](@ref), [`water_velocity`](@ref) functions. Please note that the last argument to [`water_velocity`](@ref) specifies the axis of projection of the velocity vector using symbolic variables `:x` or  `:z`.

```@example 1
η(x, n) = water_surface(p, x, n)
u(x, z, n) = (z < η(x, n)) * water_velocity(p, x, z, n, :x)
w(x, z, n) = (z < η(x, n)) * water_velocity(p, x, z, n, :z)
v(x, z, n) = sqrt.(u(x, z, n)^2 + w(x, z, n)^2)
nothing # hide
```

We need a grid to present results.

```@example 1
x = range(start = 0, stop = ℓ, length = 21)
z = range(start = -d, stop = H, length = 11)
nothing # hide
```

Finally we are ready to plot results.

```@example 1
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
```

## Linear shoaling 

We are modelling linear and regular waves of length `L` and height `H` climbing up a slope. We apply a linear wavemaker at both sides of the domain, while we consider only a half of the domain of length `ℓ`.

```@example 2
using SpectralWaves
using CairoMakie # plotting package

L = 5.0 # wavelength (m)
H = 0.05 # wave height (m)
d = 1.0 # water depth (m)
ℓ = 120.0 # fluid domain length (m)
nothing # hide
```

Again, we need a wave period `T`.

```@example 2
k = 2π / L # wave number (rad/m)
ω = sqrt(g * k * tanh(k * d)) # angular wave frequency (rad/s)
T = 2π / ω # wave period (s)
nothing # hide
```

We define a number of numerical model parameters. In order to secure a smooth start of the wavemaker paddle, we define a number of ramped wave periods `nT₀` in addition to a total number of simulated wave periods `nT`.

```@example 2
ℐ = 120 # number of harmonics
nT = 25 # number of simulated wave periods
nT₀ = 3 # number of ramped wave periods
nΔt = 50 # number of time steps per wave period
Δt = T / nΔt # time step (s)
t₀ = 0.0 # initial time (s)
τ = nT * T # total simulation time (s)
t = range(start = t₀, stop = τ, step = Δt) # time range
nothing # hide
```

We initialize wave problem `p` with `M_b=40`.

```@example 2
p = Problem(ℓ, d, ℐ, t; M_b=40)
nothing # hide
```

We use [`linear_wavemaker!`](@ref) function to define wavemaker paddle motion.

```@example 2
linear_wavemaker!(p, H, T, L, nT₀)
nothing # hide
```

The slope of height `h` is introduced using [`bottom_slope!`](@ref) function.

```@example 2
h = 0.95d
bottom_slope!(p, h)
nothing # hide
```

And we solve the problem.

```@example 2
solve_problem!(p)
nothing # hide
```

We calculate free surface elevation and bottom surface position using [`water_surface`](@ref) and [`bottom_surface`](@ref) functions

```@example 2
η(x, n) = water_surface(p, x, n)
β(x) = bottom_surface(p, x)
nothing # hide
```

and we are ready to animate results and see how the waves shoal.


```@example 2
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
record(fig, "shoaling.mp4", lastindex(t)-nΔt+1:lastindex(t);
        framerate = 50) do n
    η₀[] = η.(x, n)
end
nothing # hide
```

```@raw html
<video width="auto" controls autoplay loop>
<source src="../shoaling.mp4" type="video/mp4">
</video>
```

## Wave transformation at a step

We are modelling linear and regular waves of length `L` and height `H` undergoing a transformation at an underwater step. We apply a linear wavemaker at both sides of the domain, while we consider only a half of the domain of length `ℓ`.

```@example 3
using SpectralWaves
using CairoMakie # plotting package

L = 5.0 # wavelength (m)
H = 0.05 # wave height (m)
d = 1.0 # water depth (m)
ℓ = 60.0 # fluid domain length (m)
nothing # hide
```

We need a wave period `T`.

```@example 3
k = 2π / L # wave number (rad/m)
ω = sqrt(g * k * tanh(k * d)) # angular wave frequency (rad/s)
T = 2π / ω # wave period (s)
nothing # hide
```

We define a number of numerical model parameters (cf. linear shoaling example).

```@example 3
ℐ = 120 # number of harmonics
nT = 12 # number of simulated wave periods
nT₀ = 3 # number of ramped wave periods
nΔt = 50 # number of time steps per wave period
Δt = T / nΔt # time step (s)
t₀ = 0.0 # initial time (s)
τ = nT * T # total simulation time (s)
t = range(start = t₀, stop = τ, step = Δt) # time range
nothing # hide
```

We initialize wave problem `p` with `M_b=40`.

```@example 3
p = Problem(ℓ, d, ℐ, t; M_b=40)
nothing # hide
```

We use [`linear_wavemaker!`](@ref) function to define wavemaker paddle motion.

```@example 3
linear_wavemaker!(p, H, T, L, nT₀)
nothing # hide
```

The step of height `h` is introduced using [`bottom_step!`](@ref) function.

```@example 3
h = 0.8d
bottom_step!(p, h)
nothing # hide
```

We solve the problem.

```@example 3
solve_problem!(p)
nothing # hide
```

We calculate free surface elevation and bottom surface position using [`water_surface`](@ref) and [`bottom_surface`](@ref) functions

```@example 3
η(x, n) = water_surface(p, x, n)
β(x) = bottom_surface(p, x)
nothing # hide
```

and we animate the results to see how the waves undergo transformation at a step.


```@example 3
x = range(start = ℓ / 8, stop = 3ℓ / 8, length = 501) # spatial range
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
record(fig, "step_transformation.mp4", lastindex(t)-nΔt+1:lastindex(t);
        framerate = 50) do n
    η₀[] = η.(x, n)
end
nothing # hide
```

```@raw html
<video width="auto" controls autoplay loop>
<source src="../step_transformation.mp4" type="video/mp4">
</video>
```
