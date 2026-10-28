using Gradus, Plots

# Changing definition of inner_radius to not throw and error when a is greater than 1.
function Gradus.inner_radius(m::Gradus.KerrMetric{T}) where {T}     
    disc = m.M^2 - m.a^2 # Operand of the sqrt
    disc < 0 ? 4eps(T) : m.M + √disc # "If sqrt operand is less than 0 then return as close to zero as possible but not zero
end


# Metric and initials
m = KerrMetric(M = 1.0, a = 1.968)
x = SVector(0.0, 2.0, deg2rad(70.0), 0.0)
v  = [SVector(0.0, -1.0, 0.0,  0.0), SVector(0.0,  1.0, 0.0,  0.0)]  
xs = fill(x, length(v))


# Solutions to ODEs
sols = tracegeodesics(m, xs, v, (0.0, 10.0); μ = 0.0)

p = plot(xlabel = "r", 
ylabel = "t"; 
title = "Spacetime Diagram of Kerr Metric a = $(m.a)",
ylims = (0, 15),
xlims = (0.5, 6),
fontfamily="Computer Modern",
linewidth=2.0,
framestyle=:zerolines,
label="Light",
grid=false,
legend=true
)

# Plotting each solution in turn
for s in sols.u
    r = [u[2] for u in s.u]
    t = [u[1] for u in s.u]
    
    # Plots forwards in time then fills in to the backwards in time.
    plot!(p, r, t;
        color = :black,
        label = false
    )
    # TODO: Create a method to solve for actual negative time coords
    # plot!(p, r, -t; color = :black, label = false) # Assumes the metric is time symmetric
end

if m.a <= 1.0
    vline!(p, [Gradus.inner_radius(m)]; linestyle=:dash, color=:red, label="Event Horizon")
end

r_ergo = m.M + sqrt(m.M^2 - m.a^2 * cos(x[3])^2)
vline!(p, [r_ergo]; linestyle=:dash, color=:blue, label="Ergosurface")

paths = plot_paths_3d(sols, legend=false, extent = 5.0, t_span = 100.0, title = "Cartesian Traced Path", linewidth = 3.0)
phi_circle = range(0, 2π; length=361)
plot!(paths,
    r_ergo .* cos.(phi_circle),
    r_ergo .* sin.(phi_circle);
    color=:blue,
    linestyle=:dash,
    linewidth=2,
    label="Ergosurface (equatorial section)",
)

if m.a <= 1.0
    plot_horizon!(m)
end

combined = plot(p, paths; layout=(1, 2), size=(1200, 600))
display(combined)