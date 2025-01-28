# figure 1 - Scheme of wave propagation over topography

using CairoMakie
include("../../src/plotting.jl")

# Figure size parameters:
aspect = 2 # aspect ratio (width/height)
width = 600 # figure width (px)
height = width / aspect # figure height (px)

# Input data:
ℓ = 2π # fluid domain length (m)
d = 1 # fluid depth (m)
x = range(0, ℓ, length = 100) # x-coordinates (m)
η(x) = 0.05 * sin(4 * x) * (cos(x) - 1) # free surface elevation (m)
β(x) = - 0.05 * cos(1.5 * x) * (cos(x) - 1) - d # bottom topography (m)

# Plot:
set_theme!(theme_latexfonts())
update_theme!(fontsize=16)
fig = Figure(size = (width, height))
ax = Axis(fig[1, 1],
    limits=(-0.01ℓ, 1.01ℓ, -1.1d, 0.45d),
    aspect=aspect,)
aspect2 = 1.02ℓ / 1.4d / aspect
Δ = 3 # label offset
nx = 48 # η label position
nx2 = 68 # β label position
band!(ax, x, η.(x), β.(x), color=water_bulk) # water bulk
band!(ax, x, β.(x), -1.1d, color=sand_bulk) # bottom bulk
lines!(ax, x, η.(x), color=:black, linewidth = 1) # free surface
lines!(ax, x, β.(x), color=:black, linewidth = 1) # bottom surface
lines!(ax, [0, 0, 2π, 2π, 0], [-d, 0, 0, -d, -d],
    color=:black, linestyle=:dash, linewidth = 0.7) # fluid domain rectangle
lines!(ax, [0, 0.7d], [0, 0], color = :black, linewidth = 2) # x-axis
lines!(ax, [0, 0], [0, 0.7d / aspect2], color = :black, linewidth = 2) # z-axis
scatter!(ax, [0.8d, 0], [0, 0.8d / aspect2],
    color=:black, markersize = 10,
    marker=arrow_head,
    rotation=[3π / 2, 0]) # axis arrows
scatter!(ax, [0, 0, 2π, 2π], [-d, 0, 0, -d],
    color=:white, markersize=8, strokewidth=0.5) # fluid domain rectangle corners
text!(ax, [0, 0, 2π, 2π], [-d, 0, 0, -d];
    text=[L"(0, -d)", L"(0, \, 0)", L"(ℓ, \, 0)", L"(ℓ, -d)"],
    align=[(:left, :bottom), (:left, :top), (:right, :top), (:right, :bottom)],
    offset=[(Δ, Δ), (Δ, -Δ), (-Δ, -Δ), (-Δ, Δ)],
    ) # fluid domain rectangle corners labels
text!(ax, [x[nx], x[nx2]], [η(x[nx]), β(x[nx2])];
    text=[L"\eta(x)", L"\beta(x)"],
    align=[(:left, :bottom), (:left, :bottom)],
    offset=[(Δ, Δ), (Δ, Δ)]) # water and bottom surface label
text!(ax, [0.8d, 0], [0, 0.8d / aspect2],
    text=[L"x", L"z"],
    align=[(:left, :bottom), (:left, :bottom)],
    offset=[(0, Δ), (Δ, 0)]) # axis labels
hidedecorations!(ax)
hidespines!(ax)
