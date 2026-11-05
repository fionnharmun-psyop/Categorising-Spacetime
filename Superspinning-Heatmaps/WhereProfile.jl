using Gradus, Plots, LaTeXStrings
include("common.jl")

const E_LINE = 6.4 # Rest-frame iron line energy (keV)


# r_out must match maxr in LineProfile so the render and the profile see the same disc.
function RenderDisc(spin, position::SVector{4,T}, emissivity_index, size; r_out = 30.0,
    αlims = (-33, 33), βlims = (-22, 22)) where {T}
    m = KerrMetric(M = 1.0, a = spin)
    d = ThinDisc(Gradus.isco(m), r_out)

    
    image_height = round(Int, size * (βlims[2] - βlims[1]) / (αlims[2] - αlims[1]))
    # Using prerender since we can do it once then apply pointfunctions seperately
    α, β, img = prerendergeodesics(
    m, position, d,
    2 * position[2], # max integration time: must exceed the observer distance or rays stop before reaching the disc
    αlims = αlims,
    βlims = βlims,
    image_width = size,
    image_height = image_height, # Dimensions of image depend on the ratio of impact parameters.
    verbose = true,
    )

    # Finding the redshift and radius endpoint of each ray since that is what emissivity and flux depend on.
    # r * abs(sin(θ)) is the distance from the spin axis
    redshift = apply(ConstPointFunctions.redshift(m, position) ∘ ConstPointFunctions.filter_intersected(), img)
    radius = apply(PointFunction((m, gp, t) -> Gradus._equatorial_project(gp.x)) ∘ ConstPointFunctions.filter_intersected(), img)

  
    # Cutting at the ISCO and filtering out rays that don't intersect the disc
    valid = isfinite.(radius) .& (radius .> Gradus.isco(m)) .& isfinite.(redshift)
    flux = fill(NaN, Base.size(radius))
    
    # Flux calculation: F = r^(-emissivity_index) * g^3.
    flux[valid] .= radius[valid] .^ (-emissivity_index) .* redshift[valid] .^ 3

    if !any(isfinite, flux)
        error("Flux contains no finite values; check the ray-intersection results.")
    end

    redshift[.!valid] .= NaN
    return (; α, β, redshift, radius, flux, image_height)
end

# Redshift map with log-spaced contours of constant per-pixel flux.
function FluxHeatmap(spin, position::SVector{4,T}, emissivity_index, size; colour_min = 0.6, colour_max = 1.2, contour_levels = 7) where {T}
    (; α, β, redshift, flux, image_height) = RenderDisc(spin, position, emissivity_index, size)

    p = heatmap(α, β, redshift;
    title = "Redshift Map, a = $spin, ϑ = $(round(rad2deg(position[3]), digits = 1))°",
    color = :auto,
    clims = (colour_min, colour_max),
    colorbar_title = "Redshift (g)",
    grid = false,
    size = (size / 2, image_height / 2),
    fontfamily = "Computer Modern",
    framestyle = :box,
    dpi = 600
    )

    # Due to contour!() weirdness using the same colour range as the heatmap we need to normalise contour values to the colour range.
    # Log spacing the levels.
    log_flux = log10.(flux)
    z_min, z_max = extrema(filter(isfinite, log_flux))
    if z_min == z_max
        error("Flux is constant across all finite pixels; cannot generate contour levels.")
    end
    normalised_levels(z) = colour_min + (z - z_min) / (z_max - z_min) * (colour_max - colour_min) # fixing contours to colour range

    levels = normalised_levels.(range(z_min, z_max, length = contour_levels)[2:end-1])
    contour!(p, α, β, normalised_levels.(log_flux);
    levels = levels,
    c = cgrad([:grey, :grey]),
    clims = (colour_min, colour_max), # same as the heatmap
    colorbar_entry = false,           
    linewidth = 1.0,
    label = false,
    )

    display(p)
    return p
end

# line profile calculation. largely same as iteratinglineprofiles.jl
function CalculateLineProfile(spin, position::SVector{4,T}, emissivity_index; resolution = :low, r_out = 30.0) where {T}
    m = KerrMetric(M = 1.0, a = spin)

    # Resolution selector
    if resolution == :low
        bins, plane = (collect(range(0.1, 1.7, length = 100)), PolarPlane(GeometricGrid(); Nr = 300,  Nθ = 700,  r_max = 45.0))
    elseif resolution == :med
        bins, plane = (collect(range(0.1, 1.7, length = 100)), PolarPlane(GeometricGrid(); Nr = 800,  Nθ = 1600, r_max = 45.0))
    elseif resolution == :high
        bins, plane = (collect(range(0.1, 1.7, length = 200)), PolarPlane(GeometricGrid(); Nr = 2000, Nθ = 3600, r_max = 45.0))
    else error("Pick a valid resolution: :low, :med, :high.")
    end

    shift, flux = lineprofile(
        bins,
        r -> r^(-emissivity_index),
        m,
        position,
        ThinDisc(0.0, Inf);
        verbose = true,
        method = BinningMethod(),
        plane = plane,
        minrₑ = Gradus.isco(m),
        maxrₑ = r_out,
    )
    return UnitChange(shift, flux)
end

"""
Disc render coloured by observed energy with contours at chosen energies, plus the line
profile with vertical lines at the same energies, so each contour shows where on the disc
that part of the profile comes from.
"""
function WhereProfile(spin, position::SVector{4,T}, emissivity_index, size;
   num_energies = 5, # number of energy contours to show
   spec_energies = nothing, # specific energies to show
   resolution = :low,
   αlims = (-5, 5),
   βlims = (-5, 5),
   oversample = 2 # output pixels per traced ray; >1 doesn't add disc detail but keeps contour lines and text crisp
   ) where {T}


    (; α, β, redshift, image_height) = RenderDisc(spin, position, emissivity_index, size; 
    αlims = αlims, βlims = βlims
    )

    # Plotting observed energy on disc render and profile
    E_img = E_LINE .* redshift
    E_min, E_max = extrema(filter(isfinite, E_img))
    
    if spec_energies !== nothing
        energies = spec_energies
    else
        energies = round.(range(E_min, E_max, length = num_energies + 2)[2:end-1], digits = 1)
    end 

    line_colours = palette(:cool, max(length(energies), 2))
    clims = (E_min, E_max) # Shared colour scale for the disc and the profile

    style = (
        fontfamily = "Computer Modern", 
        framestyle = :box, 
        grid = false,
        guidefontsize = 13, 
        tickfontsize = 11, 
        titlefontsize = 13
    )

    ϑ = round(rad2deg(position[3]), digits = 1)

    # Figure sizing: the canvas is fixed in points (so fonts, which are in points, keep the same
    # proportions), and dpi sets the pixel count. Output pixels = canvas × dpi / 100, so this dpi
    # gives `oversample` output pixels per traced ray across the disc. Text and lines scale with it.
    plot_w = 600                                          # points across the disc plot area
    canvas_disc = (plot_w + 150, plot_w * image_height / size + 100) # + room for colourbar, title, labels
    canvas_prof = (plot_w + 150, 0.75 * (plot_w + 150))
    dpi = round(Int, 100 * oversample * size / plot_w)

    # ---- Disc ----
    p_disc = heatmap(α, β, E_img;
        c = :inferno, # TODO: Colourmap of observed energy should be different to that we normally use for redshift; it's confusing.
        clims = clims,
        colorbar_title = "Observed energy (keV)",
        xlabel = L"\alpha", 
        ylabel = L"\beta",
        title = "a = $spin, ϑ = $(ϑ)°, γ = $emissivity_index",
        aspect_ratio = 1, 
        xlims = extrema(α), 
        ylims = extrema(β),
        size = canvas_disc,
        dpi = dpi,
        style...
    )

    for (i, E) in enumerate(energies)
        contour!(p_disc, α, β, E_img;
            levels = [E], 
            c = cgrad([line_colours[i], line_colours[i]]), # No gradient across the line since const
            colorbar_entry = false, 
            linewidth = 1.5, 
            label = false
        )
        plot!(p_disc, [NaN], [NaN]; color = line_colours[i], label = "$E keV")
    end

    

    # ---- Line profile ----
    E_prof, F_prof = CalculateLineProfile(spin, position, emissivity_index; resolution = resolution)
    p_prof = plot(E_prof, F_prof;
        xlabel = L"Photon Energy $E$, keV", 
        ylabel = L"Normalised Flux ($E^2\,\mathrm{d}N_{\gamma}/\mathrm{d}E$)",
        title = "Line profile",
        legend = :bottomright,
        xlims = (0.5, 11.0), ylims = (0.0, 1.05),
        tickdirection = :in, minorticks = 5,
        color = :black,
        linewidth = 1.5,
        label = false,
        size = canvas_prof,
        dpi = dpi,
        style...
    )

    vline!(p_prof, [E_LINE]; color = :black, linestyle = :dashdot, linewidth = 2.0, label = "6.4 keV Emission Line")

    for (i, E) in enumerate(energies)
        vline!(p_prof, [E]; color = line_colours[i], linewidth = 2.0, label = "$E keV")
    end

    name = "WhereProfile_a$(spin)_i$(round(Int, ϑ))_γ$(emissivity_index)"
    display(p_disc)
    savefig(p_disc, name * "_disc.png")
    display(p_prof)
    savefig(p_prof, name * "_profile.png")
end

# size = 3600 is around 30mins. 
x1 = SVector(0.0, 10_000.0, deg2rad(60), 0.0)
x2 = SVector(0.0, 10000.0, deg2rad(80), 0.0)
x3 = SVector(0.0, 10000.0, deg2rad(5), 0.0)

WhereProfile(1.968, x2, 3, 1200; resolution = :med, num_energies = 6, αlims = (-5, 5), βlims = (-5, 5))
WhereProfile(1.968, x1, 3, 1200; resolution = :med, num_energies = 6, αlims = (-5, 5), βlims = (-5, 5))
WhereProfile(1.968, x3, 3, 1200; resolution = :med, num_energies = 6, αlims = (-5, 5), βlims = (-5, 5))

WhereProfile(1.05, x2, 3, 1200; resolution = :med, num_energies = 6, αlims = (-5, 5), βlims = (-5, 5))

WhereProfile(0.0, x3, 3, 1200; resolution = :med, num_energies = 6, αlims = (-33, 33), βlims = (-22, 22))
WhereProfile(0.0, x2, 3, 1200; resolution = :med, num_energies = 6, αlims = (-33, 33), βlims = (-22, 22))

WhereProfile(0.998, x1, 4, 1200; resolution = :med, num_energies = 6, αlims = (-5, 5), βlims = (-5, 5))
WhereProfile(1.05, x1, 4, 1200; resolution = :med, num_energies = 6, αlims = (-5, 5), βlims = (-5, 5))

