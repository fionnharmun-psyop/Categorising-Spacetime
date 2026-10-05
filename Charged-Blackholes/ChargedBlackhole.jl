using Gradus, Plots

m = KerrNewmanMetric(M = 1.0, a = 0.0, Q = 1.0)
x = SVector(0.0, 1000.0, deg2rad(80.0), 0.0)

isco = Gradus.isco(m)
d = ThinDisc(isco, 30.0)

pf = ConstPointFunctions.redshift(m, x) ∘ ConstPointFunctions.filter_intersected()

α, β, img = rendergeodesics(
    m, x, d,
    2000.0,
    αlims = (-10, 10),
    βlims = (-6, 6),
    image_width = 800,
    image_height = 800,
    verbose = true,
    pf = pf,
)

default(
    fontfamily = "Computer Modern",   
    framestyle = :box,                
    grid = true,
    tickdirection = :in,              
    minorticks = 5,                   
    guidefontsize = 11,               
    tickfontsize = 9,
    titlefontsize = 11,
    colorbar_titlefontsize = 11,
    size = (700, 420),
    dpi = 300,
)


heatmap(α, β, img;
    xlabel = "α",
    ylabel = "β",
    colorbar_title = "Redshift (g)",
    title = "Kerr-Newman Blackhole, Q = $(m.Q), θ = $(round(rad2deg(x[3]), digits=1))°"

)
