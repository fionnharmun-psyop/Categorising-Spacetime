using Gradus, Plots, Random
include("common.jl")

# Recreating the style of Luminet (1979): white dots on black, dot density proportional to the
# bolometric flux of a thin disc. No colour, contours or axes.


# Page & Thorne flux for Schwarzschild (M = 1), up to a constant. Zero at the ISCO (r = 6), peaks near r ≈ 9.5.
# This is what Luminet used. Only kept as a check on ThinDiscFlux below.
function PageThorneFlux(r)
    x = sqrt(r)
    log_term = log(((x + √3) * (√6 - √3)) / ((x - √3) * (√6 + √3)))
    return (x - √6 + (√3 / 2) * log_term) / (r^2.5 * (r - 3))
end

# Newtonian thin disc with zero torque at the inner edge. Approximation, works for any spin
function ZeroTorqueFlux(r, isco)
    return r^(-3) * (1 - sqrt(isco / r))
end


# Circular equatorial orbits in Kerr (M = 1, prograde): angular velocity, energy and angular momentum.
# The square root stays real for a > 1, so the same formulas are used for superspinning.
OrbitΩ(r, a) = 1 / (r^1.5 + a)
OrbitE(r, a) = (r^1.5 - 2 * sqrt(r) + a) / (r^0.75 * sqrt(r^1.5 - 3 * sqrt(r) + 2a))
OrbitL(r, a) = (r^2 - 2a * sqrt(r) + a^2) / (r^0.75 * sqrt(r^1.5 - 3 * sqrt(r) + 2a))

# Thin disc flux (Page & Thorne / Novikov & Thorne) for any spin, done numerically, up to a constant:
#   F(r) = (-dΩ/dr) / (E - ΩL)^2 / r * ∫ (E - ΩL) dL/dr dr, integrated from the ISCO to r.
# Zero torque at the ISCO is assumed, including for a > 1 where that is an extension, not a derived result.
# Returns a function of radius. For a = 0 it should match PageThorneFlux (see CheckThinDiscFlux).
function ThinDiscFlux(spin, isco; r_max = 1e4, n_radii = 4000)
    radii = exp.(range(log(isco), log(r_max), length = n_radii)) # log spaced: fine near the ISCO, coarse far out

    Ω = OrbitΩ.(radii, spin)
    E = OrbitE.(radii, spin)
    L = OrbitL.(radii, spin)

    # dL/dr by finite differences on the grid
    dL = similar(L)
    dL[1] = (L[2] - L[1]) / (radii[2] - radii[1])
    dL[end] = (L[end] - L[end-1]) / (radii[end] - radii[end-1])
    for i in 2:n_radii-1
        dL[i] = (L[i+1] - L[i-1]) / (radii[i+1] - radii[i-1])
    end

    # Running trapezium integral from the ISCO outwards
    integrand = (E .- Ω .* L) .* dL
    integral = zeros(n_radii)
    for i in 2:n_radii
        integral[i] = integral[i-1] + 0.5 * (integrand[i] + integrand[i-1]) * (radii[i] - radii[i-1])
    end

    dΩ = -1.5 .* sqrt.(radii) ./ (radii .^ 1.5 .+ spin) .^ 2
    table = (-dΩ) ./ (E .- Ω .* L) .^ 2 ./ radii .* integral

    # Linear interpolation in the table. Radii past the ends take the end value
    function flux_at(r)
        r <= radii[1] && return 0.0
        r >= radii[end] && return table[end]
        i = searchsortedlast(radii, r)
        fraction = (r - radii[i]) / (radii[i+1] - radii[i])
        return table[i] + fraction * (table[i+1] - table[i])
    end
    return flux_at
end

# The numerical flux divided by Luminet's formula should print the same number at every radius.
function CheckThinDiscFlux()
    numerical = ThinDiscFlux(0.0, 6.0)
    for r in [6.5, 8.0, 9.5, 15.0, 30.0, 100.0]
        println("r = $r   numerical / Page & Thorne = ", numerical(r) / PageThorneFlux(r))
    end
end


# Same render as FluxRadiusEnergy.jl, but with a bolometric flux: emitted flux * g^4.
# flux_law = :thin_disc is the Page & Thorne flux for any spin (Luminet's own law at a = 0).
# flux_law = :power is r^(-emissivity_index), brightest at the inner edge (not Luminet-like, here for comparison).
function RenderDisc(spin, position::SVector{4,T}, size; r_out = Inf,
    αlims = (-47, 47), βlims = (-23, 26), gap_images = true,
    flux_law = :thin_disc, emissivity_index = 3) where {T}
    m = KerrMetric(M = 1.0, a = spin)
    isco = Gradus.isco(m)
    # gap_images = true keeps the hole inside the ISCO, which gives the ring inside the shadow
    d = gap_images ? ThinDisc(isco, r_out) : ThinDisc(0.0, Inf)

    image_height = round(Int, size * (βlims[2] - βlims[1]) / (αlims[2] - αlims[1]))
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
    valid = isfinite.(radius) .& (radius .> isco) .& (radius .< r_out) .& isfinite.(redshift)

    # Emitted flux at the radius each ray lands on
    if flux_law == :thin_disc
        emitted = ThinDiscFlux(spin, isco).(radius[valid])
    elseif flux_law == :zero_torque
        emitted = ZeroTorqueFlux.(radius[valid], isco)
    elseif flux_law == :power
        emitted = radius[valid] .^ (-emissivity_index)
    else error("Pick a valid flux law: :thin_disc, :zero_torque, :power.")
    end

    # Bolometric flux: g^4, not the g^3 we use for the line
    flux = fill(NaN, Base.size(radius))
    flux[valid] .= emitted .* redshift[valid] .^ 4

    if !any(isfinite, flux)
        error("Flux contains no finite values; check the ray-intersection results.")
    end

    return (; α, β, flux, image_height)
end


# Places dots with density proportional to flux^compress.
# compress = 1.0 is the true flux, which is what Luminet plotted.
# floor adds a faint even sprinkle over the whole disc, like the sparse outer dots in Luminet's drawing.
# It is cosmetic (a fraction of the peak flux), so leave it at 0 for anything quantitative.
function SampleFluxDots(α, β, flux; n_dots = 120_000, compress = 1.0, floor = 0.0)
    disc_pixels = findall(isfinite, flux)
    weight = (flux[disc_pixels] ./ maximum(flux[disc_pixels])) .^ compress
    weight = max.(weight, floor)

    # Drawing straight from the cumulative weights. Unlike rejection sampling this never stalls,
    # even when nearly all the flux is in a few pixels (which happens for superspinning)
    cumulative = cumsum(weight)
    pixel_w, pixel_h = α[2] - α[1], β[2] - β[1]

    dot_α, dot_β = zeros(n_dots), zeros(n_dots)
    for i in 1:n_dots
        pick = min(searchsortedfirst(cumulative, rand() * cumulative[end]), length(disc_pixels))
        index = disc_pixels[pick]
        dot_α[i] = α[index[2]] + (rand() - 0.5) * pixel_w # jitter inside the pixel so dots don't sit on a lattice
        dot_β[i] = β[index[1]] + (rand() - 0.5) * pixel_h
    end
    return dot_α, dot_β
end


"""
Luminet-style image: dots only, density proportional to bolometric flux.

Pass `render` (the output of RenderDisc) to redraw with different dots without ray tracing again.
"""
function LuminetPlot(spin, position::SVector{4,T}, size;
   αlims = (-47, 47), # estimated from Luminet's drawing: the frame is about 6 shadow-widths across
   βlims = (-23, 26), # and the hole sits just below the middle
   r_out = Inf, # his disc has no visible outer edge
   gap_images = true,
   flux_law = :thin_disc,
   emissivity_index = 3, # only used by flux_law = :power
   render = nothing,
   n_dots = 120_000,
   markersize = 0.5, # tune together with n_dots: many small dots look closest to the original
   compress = 1.0,
   floor = 0.0,
   ink = :white,
   paper = :black, # swap ink and paper for a printable negative
   show_title = false,
   tag = "", # added to the filename so different settings don't overwrite each other
   seed = 1 # fixed so the same dots come out each run
   ) where {T}

    if render === nothing
        render = RenderDisc(spin, position, size;
        r_out = r_out, αlims = αlims, βlims = βlims, gap_images = gap_images,
        flux_law = flux_law, emissivity_index = emissivity_index
        )
    end
    (; α, β, flux) = render

    Random.seed!(seed)
    dot_α, dot_β = SampleFluxDots(α, β, flux; n_dots = n_dots, compress = compress, floor = floor)

    ϑ = round(rad2deg(position[3]), digits = 1)
    canvas_w = 900

    # No ticks, colourbar or legend: the dots are the whole picture
    p = scatter(dot_α, dot_β;
        color = ink,
        markersize = markersize,
        markerstrokewidth = 0,
        background_color = paper,
        title = show_title ? "a = $spin, ϑ = $(ϑ)°" : "",
        titlefontcolor = ink,
        fontfamily = "Computer Modern",
        framestyle =:none,
        ticks = false,
        grid = false,
        legend = false,
        colorbar = false,
        aspect_ratio = 1,
        xlims = extrema(α),
        ylims = extrema(β),
        size = (canvas_w, canvas_w * length(β) / length(α)), # shape taken from the render itself, so a stored render can't mismatch
        dpi = 400
    )

    name = "Luminet_a$(spin)_i$(round(Int, ϑ))_$(flux_law)_improved" * tag
    display(p)
    savefig(p, name * ".png")
    return p
end


observer(angle) = SVector(0.0, 10_000.0, deg2rad(angle), 0.0)

size = 2400 # 800 for working plots, 2400 for final

luminet_paper = RGB(0.07, 0.09, 0.09) # his print is a dark grey-green, not pure black
luminet_look = (n_dots = 120_000, markersize = 0.5, paper = luminet_paper)

# Run this once first: every line should print the same ratio
CheckThinDiscFlux()

LuminetPlot(1.968, observer(80), 800;
    flux_law = :zero_torque,
    r_out = 30.0,
    αlims = (-3.1, 4.9),
    βlims = (-2.4, 5.6),
    n_dots = 150_000,
    markersize = 0.8,
    compress = 1.0,
    show_title = true,
    tag = "_zoom4_centred"
)