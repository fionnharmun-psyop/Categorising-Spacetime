using Gradus, Plots, LaTeXStrings
include(joinpath(@__DIR__, "..", "common.jl"))

# Primary function
function DoLineProfile(spin, model, angle, resolution)

    metric = KerrMetric(M=1.0, a = spin)
    x = SVector(0.0, 10_000.0, deg2rad(angle), 0.0)
    disc = ThinDisc(0.0, Inf)

    # Resolution selector. :low for working plots, :high for final (30mins)
    if resolution == :low
        bins, plane = (collect(range(0.1, 1.7, length = 100)), PolarPlane(GeometricGrid(); Nr = 300,  Nθ = 700,  r_max = 45.0))
    elseif resolution == :med 
        bins, plane = (collect(range(0.1, 1.7, length = 100)),PolarPlane(GeometricGrid(); Nr = 800,  Nθ = 1600, r_max = 45.0))
    elseif resolution == :high 
        bins, plane = (collect(range(0.1, 1.7, length = 300)),PolarPlane(GeometricGrid(); Nr = 1600, Nθ = 3200, r_max = 45.0))
    else error("Pick a valid resolution: :low, :med: high.")
    end

    # Model Config
    @info "Computing Paths"
    sols = tracegeodesics(metric, model, 2000.0, n_samples = 10)
    model_plot = plot_paths_3d(sols, legend=false, extent = 10, t_span = 100.0)
    spin <= 1.0 && plot_horizon_3d!(model_plot, metric)
    display(model_plot)

    profile = emissivity_profile(metric, disc, model)
    rmin, rmax = Gradus.isco(metric), 30.0

    # Tabulate the emissivity once on a log-spaced grid and interpolate.
    # Passing a RingCorona profile straight to lineprofile evaluates it ray by ray,
    # which takes ~20 mins and then errors.
    r_grid = collect(exp.(range(log(rmin), log(rmax), length = 300)))
    ε_grid = emissivity_at(profile, r_grid)
    function ε(r)
        i = clamp(searchsortedlast(r_grid, r), 1, length(r_grid) - 1)
        w = (log(r) - log(r_grid[i])) / (log(r_grid[i+1]) - log(r_grid[i]))
        (1 - w) * ε_grid[i] + w * ε_grid[i+1]
    end


    # Plot style
    p = plot(;
    xlabel = L"Photon Energy $E$, keV",
    ylabel = "Arbitrarily normalised flux",
    legend = :topleft,
    xlims = (0.5, 11.0),
    ylims = (0.0, 1.05),
    size = (600, 600),
    dpi = 600,
    fontfamily = "Computer Modern",
    framestyle = :box,
    grid = true,
    tickdirection = :in,
    minorticks = 5,
    guidefontsize = 11,
    tickfontsize = 9,
    titlefontsize = 11,
    )

    # Plotting the rest frame iron-line
    E_line = 6.4
    vline!([E_line], color = :black, linestyle = :dashdot, label = "6.4KeV Emission Line")

    
    @info "Computing Profile"
    shift, flux = lineprofile(
    bins,
    ε,
    metric,
    x,
    disc;
    method = BinningMethod(),
    plane = plane,
    minrₑ = rmin,
    maxrₑ = rmax,
    verbose = true
    )

    E, F = UnitChange(shift, flux)

    plot!(
        p, E, F,  
        linewidth = 2.0,
        label = false,
    )

    # filename = "LineProfile_spin-$(join(spin, "-"))" * "_$(join(model, "-"))" * "_angle-$(angle)_resolution-$(resolution).png"
    
    display(p)
    # savefig(p, filename)
end

model = LampPostModel()
model = RingCorona()
DoLineProfile(1.998, model , 75.0 , :low)