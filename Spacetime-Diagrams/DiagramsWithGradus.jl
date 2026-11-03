using Gradus, Plots

# Changing definition of inner_radius to not throw and error when a is greater than 1.
function Gradus.inner_radius(m::Gradus.KerrMetric{T}) where {T}     
    disc = m.M^2 - m.a^2 # Operand of the sqrt
    disc < 0 ? 4eps(T) : m.M + √disc # "If sqrt operand is less than 0 then return as close to zero as possible but not zero
end


function RenderTrace(spin, x, direction ; λ_max = 25.0, rmax = 6, tmax = 15, extent = 3.0)

    # Metric and initials
    m = KerrMetric(M = 1.0, a = spin)
    v  = [direction]  
    xs = fill(x, length(v))


    # ------ Solutions to ODEs ------ 
    sols = tracegeodesics(m, xs, v, (0.0, λ_max); μ = 0.0)

    # ------ Spacetime Diagram ------ 
    p = plot(xlabel = "r", 
    ylabel = "t"; 
    title = "Spacetime Diagram of Kerr Metric a = $(m.a)",
    ylims = (0, tmax),
    xlims = (-rmax, rmax),
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

    # Event Horizon Condition
    if m.a <= 1.0
        vline!(p, [Gradus.inner_radius(m)]; linestyle=:dash, color=:red, label="Event Horizon")
    end

    # Ergosurfance Condition
    ergo_disc = m.M^2 - m.a^2 * cos(x[3])^2
    if ergo_disc >= 0
        r_ergo = m.M + sqrt(ergo_disc)
        vline!(p, [r_ergo]; linestyle=:dash, color=:blue, label="Ergosurface")
    end

    # ------ Trace ------ 
    paths = plot_paths_3d(sols, 
    legend=false, 
    extent = extent, 
    t_span = 100.0, 
    n_points = 20000, 
    title = "Cartesian Traced Path", 
    linewidth = 3.0)

    if m.a <= 1.0
        plot_horizon_3d!(paths, m)
    end

    # ------ Potential ------

    angle = x[3]
    velocity = Gradus.constrain_all(m, x, v[1], 0.0)
    # Energy (measured from infinity)
    E = Gradus.E(m, x, velocity)
    # Angular Momentum (about the axis of the black hole rotation)
    L = Gradus.Lz(m, x, velocity)
    # Carter Constant (??)
    g_θθ = Gradus.metric_components(m, SVector(x[2], x[3]))[3]
    Q = (g_θθ * velocity[3])^2 + cos(angle)^2 * (L^2 / sin(angle)^2 - m.a^2 * E^2)

    # Rescaling (?)
    l, q = L/E, Q/E^2

    # Photon potential: R(r) as a lambda function
    R = r -> (r^2 + m.a^2 - m.a*l)^2 - (r^2 - 2r + m.a^2) * (q + (l - m.a)^2)
    # Fixed r's for the x-axis
    rs = range(0, 3, length = 500)

    pot_plot = plot(xlabel = "r",
    ylabel = "R(r)";
    title = "Photon Radial Potential, a = $(m.a)",
    fontfamily="Computer Modern",
    framestyle=:box,
    grid=false,
    legend=false,
    )

    #plot!(pot_plot, rs, R.(rs))
    
    combined = plot(p, paths; layout=(1, 2), size=(1800, 600), dpi = 800)
    display(combined)
    # savefig(combined, "Trace_a$(m.a)_r$(x[2])_th$(round(Int, rad2deg(x[3]))).png")
end

v  = SVector(0.0, -1.0, eps(),  0.0)
x = SVector(0.0, 5.0, deg2rad(60), 0.0)

RenderTrace(1.55, x, v, tmax = 50)