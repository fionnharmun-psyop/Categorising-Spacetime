    using Gradus, Plots

    # Changing definition of inner_radius to not throw and error when a is greater than 1.
    function Gradus.inner_radius(m::Gradus.KerrMetric{T}) where {T}
        disc = m.M^2 - m.a^2 # Operand of the sqrt
        disc < 0 ? 4eps(T) : m.M + √disc # "If sqrt operand is less than 0 then return as close to zero as possible but not zero
    end

    m = KerrMetric(M = 1.0, a =  1.05)
    x = SVector(0.0, 3000.0, deg2rad(60.0), 0.0)

    isco = Gradus.isco(m)
    d = ThinDisc(isco, 30.0)

    pf = ConstPointFunctions.redshift(m, x) ∘ ConstPointFunctions.filter_intersected()

    α, β, img = rendergeodesics(
        m, x, d,
        20000.0,
        αlims = (-6, 6),
        βlims = (-6, 6),
        image_width = 1000,
        image_height = 1000,
        verbose = true,
        pf = pf,
    )

    default(
        fontfamily = "Computer Modern",                 
        grid = false,
        tickdirection = :in,              
        minorticks = 5,                   
        guidefontsize = 11,               
        tickfontsize = 9,
        titlefontsize = 11,
        colorbar_titlefontsize = 11,
    )

    heatmap(α, β, img;
        color = :coolwarm,
        xlabel = "α",
        ylabel = "β",
        colorbar_title = "g",
        title = "Kerr naked singularity, a = $(m.a), θ = $(round(rad2deg(x[3]), digits=1))°",
        background_color_inside = :black,
        aspect_ratio = :equal,
      size = (500, 500),
      clims = (0.1, 1.0)
      )