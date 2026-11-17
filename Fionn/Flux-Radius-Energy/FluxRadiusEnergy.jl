using Gradus, Plots, LaTeXStrings, Random
include("common.jl")

const E_LINE = 6.4 # Rest-frame iron line energy (keV)


# Same render as WhereProfile.jl. Returns the redshift, landing radius and flux of every pixel.
# gap_images = true keeps the hole inside the ISCO so rays can pass through it (what we have been doing).
# gap_images = false stops every ray at its first plane crossing (what M&I do).
function RenderDisc(spin, position::SVector{4,T}, emissivity_index, size; r_out = 30.0,
    αlims = (-5, 5), βlims = (-5, 5), gap_images = true) where {T}
    m = KerrMetric(M = 1.0, a = spin)
    d = gap_images ? ThinDisc(Gradus.isco(m), r_out) : ThinDisc(0.0, Inf)

    image_height = round(Int, size * (βlims[2] - βlims[1]) / (αlims[2] - αlims[1]))
    # Using prerender since we can do it once then apply pointfunctions seperately
    α, β, img = prerendergeodesics(
    m, position, d,
    2 * position[2], # max integration time: must exceed the observer distance or rays stop before reaching the disc
    αlims = αlims,
    βlims = βlims,
    image_width = size,
    image_height = image_height,
    verbose = true,
    )

    redshift = apply(ConstPointFunctions.redshift(m, position) ∘ ConstPointFunctions.filter_intersected(), img)
    radius = apply(PointFunction((m, gp, t) -> Gradus._equatorial_project(gp.x)) ∘ ConstPointFunctions.filter_intersected(), img)

    # Cutting at the ISCO and the outer edge, and filtering out rays that don't intersect the disc
    valid = isfinite.(radius) .& (radius .> Gradus.isco(m)) .& (radius .< r_out) .& isfinite.(redshift)
    flux = fill(NaN, Base.size(radius))

    # Flux calculation: F = r^(-emissivity_index) * g^3.
    flux[valid] .= radius[valid] .^ (-emissivity_index) .* redshift[valid] .^ 3

    if !any(isfinite, flux)
        error("Flux contains no finite values; check the ray-intersection results.")
    end

    redshift[.!valid] .= NaN
    radius[.!valid] .= NaN
    return (; α, β, redshift, radius, flux, image_height)
end


# Cuts the dark end off a colour map so no dot disappears into the black background.
# rev = true puts the bright end at small radius, where the flux is highest.
function TruncatedMap(name; low = 0.35, high = 1.0, rev = false)
    full = cgrad(name)
    return cgrad([full[z] for z in range(low, high, length = 256)], rev = rev)
end


# Luminet-style dots: random pixels are kept with probability equal to their (compressed) flux,
# so the dot density follows the flux. Returns the dot positions and the log radius each one landed at.
function SampleFluxDots(α, β, flux, radius; n_dots = 60_000, compress = 0.5)
    # Raw flux spans orders of magnitude, so compress or every dot ends up at the inner edge
    weight = (flux ./ maximum(filter(isfinite, flux))) .^ compress
    pixel_w, pixel_h = α[2] - α[1], β[2] - β[1]

    dot_α, dot_β, dot_radius = Float64[], Float64[], Float64[]

    # Rejection sampling. Slow if compress is near 1, since most pixels get rejected
    while length(dot_α) < n_dots
        index = rand(CartesianIndices(flux))
        w = weight[index]
        if isfinite(w) && rand() < w
            push!(dot_α, α[index[2]] + (rand() - 0.5) * pixel_w) # jitter inside the pixel so dots don't sit on a lattice
            push!(dot_β, β[index[1]] + (rand() - 0.5) * pixel_h)
            push!(dot_radius, log10(radius[index]))
        end
    end
    return dot_α, dot_β, dot_radius
end


# Regular grid of dots whose area follows the (compressed) flux. Goes on top of a radius heatmap.
function PlotFluxGrid!(p, α, β, flux; spacing = 30, max_size = 7, compress = 0.3)
    # Taking every `spacing`-th pixel. Rows of the image are β, columns are α
    rows = 1:spacing:length(β)
    cols = 1:spacing:length(α)

    dot_α = [α[j] for i in rows, j in cols]
    dot_β = [β[i] for i in rows, j in cols]

    weight = (flux[rows, cols] ./ maximum(filter(isfinite, flux))) .^ compress
    valid = isfinite.(weight)

    scatter!(p, dot_α[valid], dot_β[valid];
        markersize = max_size .* sqrt.(weight[valid]), # sqrt so dot area, not diameter, tracks the weight
        color = :white,
        markerstrokecolor = :black,
        markerstrokewidth = 0.5,
        label = false
    )
    return p
end


"""
Disc image showing three things at once: flux (dots), landing radius (colour) and observed
energy (contours).

dots = :density gives Luminet-style dots on black, with density following flux and colour the radius.
dots = :grid gives a radius heatmap with a regular grid of dots sized by flux.

energy_colour = :shift colours each contour by its shift from 6.4 keV (red below, blue above).
Any palette name (e.g. :cool) colours them in order of energy instead.
"""
function FluxRadiusEnergy(spin, position::SVector{4,T}, emissivity_index, size;
   dots = :density,
   num_energies = 5, # number of energy contours to show
   spec_energies = nothing, # specific energies to show
   αlims = (-5, 5),
   βlims = (-5, 5),
   gap_images = true,
   render = nothing, # pass the output of RenderDisc to skip the ray tracing, e.g. when only changing compress or colours
   n_dots = 60_000, # :density only
   markersize = 1.2, # :density only. Tune together with n_dots: too many big dots saturates the inner edge
   compress = 0.5, # flux^compress sets the dot density / size. 1 is true flux, smaller shows more of the outer disc
   spacing = 30, # :grid only, pixels between dots
   radius_colour = TruncatedMap(:inferno; rev = true),
   energy_colour = :cool,
   tag = "", # added to the filename so different settings don't overwrite each other
   seed = 1 # fixed so the same dots come out each run
   ) where {T}

    if render === nothing
        render = RenderDisc(spin, position, emissivity_index, size;
        αlims = αlims, βlims = βlims, gap_images = gap_images
        )
    end
    (; α, β, redshift, radius, flux, image_height) = render

    # Energies for the contours. Same as WhereProfile.jl
    E_img = E_LINE .* redshift
    E_min, E_max = extrema(filter(isfinite, E_img))

    if spec_energies !== nothing
        energies = spec_energies
    else
        energies = round.(range(E_min, E_max, length = num_energies + 2)[2:end-1], digits = 1)
    end

    if energy_colour == :shift
        # Red for redshifted, blue for blueshifted, white at the rest energy
        shift_map = cgrad(:RdBu)
        span = max(E_LINE - E_min, E_max - E_LINE)
        line_colours = [shift_map[0.5 + 0.5 * (E - E_LINE) / span] for E in energies]
    else
        line_colours = palette(energy_colour, max(length(energies), 2))
    end

    style = (
        fontfamily = "Computer Modern",
        framestyle = :box,
        grid = false,
        guidefontsize = 13,
        tickfontsize = 11,
        titlefontsize = 13
    )

    ϑ = round(rad2deg(position[3]), digits = 1)
    plot_w = 600
    canvas = (plot_w + 150, plot_w * image_height / size + 100) # + room for colourbar, title, labels

    shared = (
        colorbar_title = L"\log_{10} r", # kept short so it doesn't run into the tick labels
        xlabel = L"\alpha",
        ylabel = L"\beta",
        title = "a = $spin, ϑ = $(ϑ)°, γ = $emissivity_index",
        aspect_ratio = 1,
        xlims = extrema(α),
        ylims = extrema(β),
        size = canvas,
        dpi = 400
    )

    # ---- Flux and radius ----
    if dots == :density
        Random.seed!(seed)
        dot_α, dot_β, dot_radius = SampleFluxDots(α, β, flux, radius; n_dots = n_dots, compress = compress)

        p = scatter(dot_α, dot_β;
            zcolor = dot_radius,
            c = radius_colour,
            markersize = markersize,
            markerstrokewidth = 0,
            background_color_inside = :black,
            label = false,
            shared...,
            style...
        )
    elseif dots == :grid
        p = heatmap(α, β, log10.(radius);
            c = radius_colour,
            shared...,
            style...
        )
        PlotFluxGrid!(p, α, β, flux; spacing = spacing, compress = compress)
    else error("Pick a valid dot style: :density, :grid.")
    end

    # ---- Energy contours ----
    for (i, E) in enumerate(energies)
        contour!(p, α, β, E_img;
            levels = [E],
            c = cgrad([line_colours[i], line_colours[i]]), # No gradient across the line since const
            colorbar_entry = false,
            linewidth = 1.5,
            label = false
        )
        plot!(p, [NaN], [NaN]; color = line_colours[i], label = "$E keV")
    end

    name = "FluxRadEnergy_a$(spin)_i$(round(Int, ϑ))_γ$(emissivity_index)_c$(compress)_pastISCO$(gap_images)" * tag
    display(p)
    savefig(p, name * ".png")
    return p
end

observer(angle) = SVector(0.0, 10_000.0, deg2rad(angle), 0.0)

# Put it in a try loop so one bad render doesn't kill the whole batch
function TryPlot(args...; kwargs...)
    try
        FluxRadiusEnergy(args...; kwargs...)
    catch err
        @error "Plot failed; skipping it" args exception = err
    end
end



size = 800
spins = [1.2]
angles = [60.0, 80.0]
winding = [true, false]

for a in spins
    for theta in angles
        for winds in winding
            TryPlot(a, observer(theta), 4, size; 
            n_dots = 120_000, markersize = 0.5, 
            seed = 545354, 
            αlims = (-4, 4), βlims = (-4, 4), 
            compress = 0.8,
            gap_images = winds
            )
        end
    end
end


# TODO: Pick out the best plots and run them on a compression of zero and higher res. 
# TODO: Check flux calculations.
# TODO: Change flux equation to a disc emissivity profile (lampost, corona)