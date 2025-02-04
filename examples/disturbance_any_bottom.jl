using SpectralWaves
using GLMakie # plotting package

d = 1.0 # water depth (m)
ℓ = 100 # fluid domain length (m)

ℐ = 100 # number of harmonics
Δt = 0.01 # time step (s)
t₀ = 0.0 # initial time (s)
τ = 12.0 # total simulation time (s)
t = range(start = t₀, stop = τ, step = Δt) # time range

p = Problem(ℓ, d, ℐ, t)
p2 = Problem(ℓ, d, ℐ, t; M_b = 30)

h = 0.5d # bump height (m)
λ = 0.02ℓ # bump length (m)
x₀ = 0.5ℓ # bump center (m)
surface_bump!(p, h, λ)
surface_bump!(p2, h, λ)
x = range(- ℓ / 2, ℓ / 2, length = 1001)
βₙ = @. - 0.3 * cos(1.5 * x) * (cos(x) - 1)
bottom_vector!(p2, x, βₙ)

solve_problem!(p)
solve_problem!(p2)

η(x, n) = water_surface(p, x, n)
η_2(x, n) = water_surface(p2, x, n)
β(x) = bottom_surface(p2, x)
η₀ = Observable(η.(x, 1))
η₂ = Observable(η_2.(x, 1))
β₀ = Observable(β.(x) .- d)

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "x (m)", ylabel = "η (m)")
lines!(ax, x, η₀, color = :blue, linewidth = 2)
lines!(ax, x, η₂, color = :black, linewidth = 2)
lines!(ax, x, β₀, color = :yellow, linewidth = 2)
lines!(ax, x, βₙ .- d, color = :black, linewidth = 1, linestyle = :dash)
limits!(ax, - ℓ / 2, ℓ / 2, -1.5d, d)
display(fig)

for n in eachindex(t)
    η₀[] = η.(x, n)
    η₂[] = η_2.(x, n)
    sleep(0.001)
end
