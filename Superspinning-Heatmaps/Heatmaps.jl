using Gradus, Plots

function Gradus.inner_radius(m::Gradus.KerrMetric{T}) where {T}
    disc = m.M^2 - m.a^2
    disc < 0 ? 4eps(T) : m.M + √disc
end

function RenderPlot(a, size, angle; color_limits = :auto, color = :auto)

    m = KerrMetric(M = 1.0, a = a)
    x = SVector(0.0, 1000.0, angle, 0.0)
    isco = Gradus.isco(m)
    d = ThinDisc(isco, 100.0)
    pf = ConstPointFunctions.redshift(m, x) ∘ ConstPointFunctions.filter_intersected()

    αlims = (-10, 10)
    βlims = (-10, 10)
    image_height = round(Int, size * (βlims[2] - βlims[1]) / (αlims[2] - αlims[1]))

    α, β, img = rendergeodesics(
    m, x, d,
    2000.0,
    αlims = αlims,
    βlims = βlims,
    image_width = size,
    image_height = image_height,
    verbose = true,
    pf = pf
    )

    p = heatmap(α, β, img;
    title = "Redshift Map, a = $a, ϑ = $(round(rad2deg(angle), digits = 1))°",
    color = color,
    clims = color_limits,
    colorbar_title = "Redshift (g)",
    grid = false,
    size = (size / 2 , image_height / 2),
    fontfamily = "Computer Modern", 
    font_size = 11,  
    framestyle = :box,
    dpi = 800
    )

    display(p)
    savefig(p, "RedshiftMap_a$(a)_i$(round(Int, rad2deg(angle)))_zoom.png")
end



RenderPlot(1.002, 2400, deg2rad(60))
RenderPlot(1.986, 2400, deg2rad(60))
RenderPlot(1.986, 2400, deg2rad(5))
RenderPlot(1.986, 2400, deg2rad(80))
RenderPlot(1.002, 2400, deg2rad(80))