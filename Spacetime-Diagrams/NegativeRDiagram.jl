using Gradus, Plots


#function Gradus.inner_radius(m::Gradus.KerrMetric{T}) where {T}
    #disc = m.M^2 - m.a^2
    #disc < 0 ? 4eps(T) : m.M + √disc
#end

function jagged_vline!(p, x0, ymin, ymax; n = 20, A = 0.2, color = :black)
    ys = range(ymin, ymax; length = 2n + 1)
    xs = [x0 + (isodd(i) ? -A : A) for i in eachindex(ys)]
    plot!(p, xs, ys; color = color, linewidth = 2.0, label = false)
end

function RenderNegativeRTrace(spin, x, direction; λ_max = 25.0, rmin = -10.0, rmax = 10.0, tmax = 30.0)

    # Metric and initial conditions
    m = KerrMetric(M = 1.0, a = spin)
    v = [direction]
    xs = fill(x, length(v))

    # Geodesic solutions
    chart = Gradus.PolarChart(rmin, rmax)
    sols = tracegeodesics(m, xs, v, (0.0, λ_max); μ = 0.0, chart = chart)

    # Spacetime diagram
    p = plot(
        xlabel = "r",
        ylabel = "t",
        title = "Kerr Spacetime Diagram, a = $(m.a)",
        ylims = (-tmax, tmax),
        xlims = (rmin, rmax),
        fontfamily = "Computer Modern",
        linewidth = 2.0,
        framestyle = :zerolines,
        label = "Light",
        grid = false,
        legend = true,
    )

    for s in sols.u
        r = [u[2] for u in s.u]
        t = [u[1] for u in s.u]
        plot!(p, r, t; color = :black, label = false)
    end

    jagged_vline!(p, 0.0, -tmax, tmax; color = :black)

    # Radial photon potential
    angle = x[3]
    velocity = Gradus.constrain_all(m, x, v[1], 0.0)
    E = Gradus.E(m, x, velocity)
    L = Gradus.Lz(m, x, velocity)
    g_θθ = Gradus.metric_components(m, SVector(x[2], x[3]))[3]
    Q = (g_θθ * velocity[3])^2 + cos(angle)^2 * (L^2 / sin(angle)^2 - m.a^2 * E^2)

    l, q = L / E, Q / E^2
    R = r -> (r^2 + m.a^2 - m.a * l)^2 - (r^2 - 2r + m.a^2) * (q + (l - m.a)^2)
    rs = range(rmin, rmax; length = 1000)

    potential = plot(
        rs,
        R.(rs),
        xlabel = "r",
        ylabel = "R(r)",
        title = "Photon Radial Potential",
        fontfamily = "Computer Modern",
        framestyle = :box,
        grid = false,
        legend = false,
        linewidth = 2.0,
    )
    hline!(potential, [0.0]; linestyle = :dash, color = :gray)

    combined = plot(p, potential; layout = (1, 2), size = (1400, 600), dpi = 800)
    display(combined)
end

v = SVector(0.0, -1.0, eps(), 0.0)
x = SVector(0.0, 3.0, deg2rad(1), 0.0)

RenderNegativeRTrace(1.968, x, v, tmax = 15)

# TODO: Normalise potential, plot ray in real space, see if can use gradus functionality to make array of rays, condition function to check the ray endpoint: disc, -r, +r and in what order.