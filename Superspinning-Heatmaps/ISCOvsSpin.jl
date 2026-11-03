using Gradus, Plots

a = range(0.0, 4, length=10000)
isco_values = Float64[]

for spin in a
    m = KerrMetric(M = 1.0, a = spin)
    isco = Gradus.isco(m)
    push!(isco_values, isco)
end

default(
    fontfamily = "Computer Modern",   # serif, LaTeX-like look (GR backend ships with it)
    framestyle = :box,                # full box around the axes, not just left/bottom
    grid = false,
    tickdirection = :in,              # ticks point inward
    minorticks = 5,                   # 4 minor ticks between majors
    guidefontsize = 11,               # axis labels
    tickfontsize = 9,
    titlefontsize = 11,
    colorbar_titlefontsize = 11,
    size = (700, 420),
    dpi = 300,
)

plot(a, isco_values;
    xlabel = "Spin parameter (a)",
    ylabel = "ISCO radius (r)",
    color = :black,
    grid = true,
    title = "ISCO radius vs Spin parameter for Kerr black holes",
    legend = false
)