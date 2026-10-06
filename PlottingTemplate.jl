using Plots

x_data = range(0, stop=10, length=100)
y_data = sin.(x_data) .+ 0.1 .* randn(length(x_data))  

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

plot(x_data, y_data;
    xlabel = "dependent variable (a)",
    ylabel = "free variable (b)",
    color = :black,
    grid = true,
    title = "My Graph",
    label="Noise",
    legend = :topright
)

heatmap(rand(100, 100); 
    xlabel = "X-axis", 
    ylabel = "Y-axis", 
    grid = true,
    title = "Heatmap Example", 
    colorbar_title = "Intensity"
)
