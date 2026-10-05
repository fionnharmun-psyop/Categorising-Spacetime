using Gradus, Plots

Q = range(0.0, 0.9, length=100)
A = range(0.0, 0.9, length=100)

isco_values_q = [Gradus.isco(KerrNewmanMetric(M = 1.0, a = 0.0, Q = q)) for q in Q]
isco_values_a = [Gradus.isco(KerrNewmanMetric(M = 1.0, a = s, Q = 0.0)) for s in A]

default(
    fontfamily = "Computer Modern",   
    framestyle = :box,                
    grid = false,
    tickdirection = :in,              
    minorticks = 5,                   
    guidefontsize = 11,               
    tickfontsize = 9,
    titlefontsize = 11,
    colorbar_titlefontsize = 11,
    size = (700, 420),
    dpi = 300,
)

p = [isco_values_q, isco_values_a]

plot(Q, p;
    xlabel = "Charge (Q)",
    ylabel = "ISCO Radius (r)",
    grid = true,
    title = "ISCO Vs Charge for Kerr-Newman Metric",
    label=["Charge" "Spin"],
    legend = true
)
