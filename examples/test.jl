using GLMakie

points = Observable(Point2f[randn(2)])

fig, ax = scatter(points)
limits!(ax, -4, 4, -4, 4)
display(fig)

fps = 60
nframes = 120

for i = 1:nframes
    new_point = Point2f(randn(2))
    points[] = push!(points[], new_point)
    sleep(1/fps) # refreshes the display!
end
