using Gradus, Plots, LaTeXStrings

function Gradus.inner_radius(m::Gradus.KerrMetric{T}) where {T}
    disc = m.M^2 - m.a^2 # Operand of the sqrt
    disc < 0 ? 4eps(T) : m.M + √disc # "If sqrt operand is less than 0 then return as close to zero as possible but not zero
end

# Function that fixes units to conventional ones in literature
function UnitChange(g, F) 
    observed_energy = g .* 6.4
    normalised_flux = g .^2 .* F ./ maximum(g .^2 .* F)
    return observed_energy, normalised_flux
end

# Primary function
function DoLineProfile(spin, emissivity, angle, resolution)

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

    # Plot style
    p = plot(;
    xlabel = L"Photon Energy $E$, keV",
    ylabel = "Arbitrarily normalised flux",
    legend = :topleft,
    xlims = (0.5, 11.0),
    ylims = (0.0, 1.05),
    size = (800, 800),
    dpi = 800,
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

    # Pallete and styles representing different emissivity and spins
    colors = palette(:viridis, length(emissivity))
    styles = [:solid, :dash, :dot, :dashdot]

    # Enumerating plotting of each spin for each emissivity
    for (a_index, a) in enumerate(spin)
        plot!(p, [NaN], [NaN]; color = :black, linestyle = styles[a_index], label = L"a = %$a")
        for (y_index, y) in enumerate(emissivity)
    
            metric = KerrMetric(M=1.0, a = a)
            @info "Computing" a y
            shift, flux = lineprofile(
                bins,
                r -> r^(-y),
                metric, 
                x, 
                disc;
                verbose = true,
                method = BinningMethod(),
                plane = plane,
                minrₑ = Gradus.isco(metric),
                maxrₑ = 30.0
            )

            E, F = UnitChange(shift, flux)

            plot!(
                p, E, F, 
                color = colors[y_index], 
                linestyle = styles[a_index], linewidth = 1.5,
                label = a_index == 1 ? L"\gamma = %$y" : "",
            )

        end
    end

     filename = "LineProfile_spin-$(join(spin, "-"))" * "_emissivity-$(join(emissivity, "-"))" * "_angle-$(angle)_resolution-$(resolution).png"
    
     display(p)
     savefig(p, filename)
end

# TODO: Make sure Units are right and congruent with Andy's update

DoLineProfile([0.998, 1.002] , [3,4,5,6,7] , 75.0 , :low)