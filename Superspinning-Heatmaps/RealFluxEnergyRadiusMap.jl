using Gradus, Plots, LaTeXStrings, Random
include("common.jl")

# NOT COMMENTED ENOUGH BECAUSE I WROTE THIS SO LATE AT NIGHT AAHH
# TODO: Comment FluxRadiusEnergy.jl

const E_LINE = 6.4 # Rest-frame iron line energy (keV)

"Function returns the radius, energy, redshift, and flux image data of rays incident on the disc.

Flux is given in raw units and must be normalised first.
"
function RenderDisc(spin, pos, emissivity, imgsize; rmax = 30.0,
    αlims = (-5, 5), βlims = (-5, 5), looping = true, flux_compression = 1.0)
# TODO: make it so there is a rmin kwarg that allows rendering below ISCO.

    m = KerrMetric(M = 1.0, a = spin)
    rmin = Gradus.isco(m)
    d = looping ? ThinDisc(rmin, rmax) : ThinDisc(0.0, Inf)

    image_height = round(Int, imgsize * (βlims[2] - βlims[1]) / (αlims[2] - αlims[1]))

    α, β, imgdata = prerendergeodesics(
        m, pos, d,
        2000,
        αlims = αlims, βlims = βlims,
        image_width = imgsize, image_height = image_height,
        verbose = true
    )

    redshift = apply(ConstPointFunctions.redshift(m, pos) ∘ ConstPointFunctions.filter_intersected(), imgdata)
    radius = apply(PointFunction((m, gp, t) -> Gradus._equatorial_project(gp.x)) ∘ ConstPointFunctions.filter_intersected(), imgdata)

    valid = isfinite.(radius) .& (radius .> rmin) .& (radius .< rmax) .& isfinite.(redshift)
    redshift[.!valid] .= NaN
    radius[.!valid] .= NaN
    energy = E_LINE .* redshift

    flux = fill(NaN, Base.size(radius)) # Base. because the script defines a global `size`
    if emissivity isa Real
        println("Power-law flux")
        # Flux calculation: F = r^(-q) * g^3
        flux[valid] .= radius[valid] .^ (-emissivity) .* redshift[valid] .^ 3
    elseif emissivity isa Gradus.AbstractCoronaModel
        println("Coronal model flux: $(nameof(typeof(emissivity)))")
        # Trace corona → disc to get ε(r), then F = ε(r) * g^3
        ϵ = emissivity_profile(m, d, emissivity; n_samples = 2000)
        flux[valid] .= emissivity_at(ϵ, radius[valid]) .* redshift[valid] .^ 3
    else
        error("emissivity must be a power-law index, a coronal model; got $(typeof(emissivity))")
    end

    if !any(isfinite, flux)
        error("Flux contains no finite values; check the ray-intersection results.")
    end

    norm_flux = (flux ./ maximum(filter(isfinite, flux))) .^ flux_compression
    norm_flux[.!valid] .= NaN 

    return (; α, β, radius, energy, redshift, flux, norm_flux, image_height)
end

"Probabilistically sampling flux pixels based on the flux.

Requires a normalised flux data set. `colourimg` is the image whose value is recorded for each dot.
"
function FluxSample(α, β, fluximg, colourimg; numdots = 60_000, seed = 1)
    Random.seed!(seed)
    weighting = fluximg
    dot_α, dot_β, dot_img = Float64[], Float64[], Float64[]
    pixel_w, pixel_h = α[2] - α[1], β[2] - β[1]

    # Rejection sampling.
    while length(dot_α) < numdots
        ind = rand(CartesianIndices(fluximg))

        w = weighting[ind]

        if isfinite(w) && rand() < w
            push!(dot_α, α[ind[2]] + (rand() - 0.5) * pixel_w) # jitter inside the pixel so dots don't sit on a lattice
            push!(dot_β, β[ind[1]] + (rand() - 0.5) * pixel_h)
            push!(dot_img, colourimg[ind])
        end
    end

    return dot_α, dot_β, dot_img
end

function Heatmapping(α, β, img, style;
    colour = cgrad([cgrad(:inferno)[z] for z in range(0.35, 1.0, length = 256)], rev = true)
    ,markersize = 1.0)

    p = scatter(α, β;
    zcolor = img,
    c = colour,
    markersize = markersize,
    background_color_inside = :black,
    label = false,
    style...
    )
    return p
end


function Contouring(α, β, img, plot, levels; colour = :cool, width = 1.5)
    line_colours = palette(colour, max(length(levels), 2))
    for (i, E) in enumerate(levels)
        contour!(plot, α, β, img;
        levels = [E],
        c = cgrad([line_colours[i], line_colours[i]]), # No gradient across the line since const
        colorbar_entry = false,
        linewidth = width,
        label = false
    )
    plot!(plot, [NaN], [NaN]; color = line_colours[i], label = "$E keV")
    end
    return plot
end

function FluxRadiusEnergy(spin, position, emissivity, size, alims, blims; 
    rmax = 30.0, looping = true, flux_compression = 1.0 )

    render = RenderDisc(spin, position, emissivity, size; 
    rmax = rmax, 
    αlims = alims, βlims = blims, 
    looping = looping, flux_compression = flux_compression)
    
    (; α, β, radius, energy, norm_flux, image_height) = render

    dot_α, dot_β, dot_logr = FluxSample(α, β, norm_flux, log10.(radius); seed = 3435, numdots = 100_000)

    ϑ = round(rad2deg(position[3]), digits = 1)
    emiss_label = emissivity isa Real ? "y = $emissivity" : string(nameof(typeof(emissivity)))
    # No spaces for the filename; coronal height added so different lamp posts don't overwrite each other
    emiss_tag = emissivity isa Real ? "y$emissivity" :
        string(nameof(typeof(emissivity))) * (hasproperty(emissivity, :h) ? "_h$(emissivity.h)" : "")
    plot_w = 600

    style = (
        colorbar_title = L"\log_{10} r",
        xlabel = L"\alpha",
        ylabel = L"\beta",
        title = "a = $spin, ϑ = $(ϑ)°, $emiss_label",
        aspect_ratio = 1,
        xlims = extrema(α),
        ylims = extrema(β),
        markerstrokewidth = 0,
        fontfamily = "Computer Modern",
        framestyle = :box,
        grid = false,
        size = (plot_w + 150, plot_w * image_height / size + 100), # + room for colourbar, title, labels
    )

    p = Heatmapping(dot_α, dot_β, dot_logr, style; markersize = 0.5)

    # Energy contours evenly spaced between the extremes
    E_min, E_max = extrema(filter(isfinite, energy))
    levels = round.(range(E_min, E_max, length = 7)[2:end-1], digits = 1)
    Contouring(α, β, energy, p, levels)

    name = "FluxDensity_EnergyContours_RadiusColour_a$(spin)_i$(round(Int, ϑ))_$(emiss_tag)_c$(flux_compression)_looping$(looping)"
    savefig(p, name * ".png")
    return p
end



size = 400
spins = [1.2]
angles = [60.0, 80.0, 5.0]
winding = [true, false]
models = [3, 6, LampPostModel()]

for ϵ in models
    for a in spins
        for theta in angles
            for yn in winding
              x = SVector(0.0, 10_000.0, deg2rad(theta), 0.0)
             FluxRadiusEnergy(a, x, ϵ, 400, (-5, 5), (-5, 5); flux_compression = 0.5, looping = yn)
            end
        end
    end
end