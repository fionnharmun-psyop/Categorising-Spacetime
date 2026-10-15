using Gradus, Plots


# Metric and initials
m = KerrMetric(M = 1.0, a = 1.968)

x = SVector(0.0, 1.5, deg2rad(60), 0.0)

v  = [SVector(0.0, -1.0, 0.0,  0.0), SVector(0.0,  1.0, 0.0,  0.0)]  

xs = fill(x, length(v))

# Changing definition of inner_radius to not throw and error when a is greater than 1.
function Gradus.inner_radius(m::Gradus.KerrMetric{T}) where {T}     
    disc = m.M^2 - m.a^2 # Operand of the sqrt
    disc < 0 ? 4eps(T) : m.M + √disc # "If sqrt operand is less than 0 then return as close to zero as possible but not zero
end

# Solutions to ODEs
sols = tracegeodesics(m, xs, v, (0.0, 10.0); μ = 0.0)

default(
    fontfamily="Computer Modern",
    linewidth=2.0,
    framestyle=:zerolines,
    label="Light",
    grid=false,
    legend=true,
)

p = plot(xlabel = "r", 
ylabel = "t"; 
title = "Minkowski Diagram of Kerr Metric a = $(m.a)",
ylims = (-15, 15),
xlims = (0, 6)
)

# Plotting each solution in turn
for s in sols.u
    r = [u[2] for u in s.u]
    t = [u[1] for u in s.u]
    
    # Plots forwards in time then fills in to the backwards in time.
    plot!(p, r, t;
        #fillrange = -t,
        #fillcolor = :grey,
        #fillalpha = 0.25,
        color = :black,
        label = false
    )
    # To Do: Create a method to solve for actual negative time coords
    plot!(p, r, -t; color = :black, label = false) # Assumes the metric is time symmetric
end

if m.a <= 1.0
    vline!(p, [Gradus.inner_radius(m)]; linestyle=:dash, color=:red, label="Event Horizon")
end

r_ergo = m.M + sqrt(m.M^2 - m.a^2 * cos(x[3])^2)
vline!(p, [r_ergo]; linestyle=:dash, color=:blue, label="Ergosurface")

display(p)