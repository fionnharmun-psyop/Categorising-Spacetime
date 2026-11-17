using Gradus, Plots, Random
include("common.jl") # inner_radius override so a > 1 doesn't throw


# Same render as Heatmaps.jl, but with prerendergeodesics so the endpoint of every pixel's ray is kept.
# That tells us where each pixel's ray ended up, and lets us check the 3D traces against the heatmap.
function RenderHeatmap(m, x, d, size; αlims = (-10, 10), βlims = (-10, 10), λ_max = 2000.0)
    image_height = round(Int, size * (βlims[2] - βlims[1]) / (αlims[2] - αlims[1]))

    α, β, cache = prerendergeodesics(
    m, x, d,
    λ_max,
    αlims = αlims,
    βlims = βlims,
    image_width = size,
    image_height = image_height,
    verbose = true,
    ensemble = Gradus.EnsembleEndpointThreads() # what rendergeodesics uses, so this is the Heatmaps.jl image exactly
    )

    redshift = apply(ConstPointFunctions.redshift(m, x) ∘ ConstPointFunctions.filter_intersected(), cache)

    # points[i, j] is the endpoint of the ray with impact parameters (α[j], β[i])
    return (; α, β, redshift, points = cache.points, image_height)
end


# Smallest window containing every pixel whose ray missed the disc, plus some padding.
# The islands sit inside the hole in the image, so sampling here doesn't waste rays on the outer disc.
function HoleWindow(α, β, points; padding = 0.05)
    missed = [p.status != StatusCodes.IntersectedWithGeometry for p in points]
    if !any(missed)
        return extrema(α), extrema(β)
    end
    hole_α = extrema(α[j] for j in eachindex(α) if any(view(missed, :, j)))
    hole_β = extrema(β[i] for i in eachindex(β) if any(view(missed, i, :)))
    pad_α = padding * (hole_α[2] - hole_α[1])
    pad_β = padding * (hole_β[2] - hole_β[1])
    return (hole_α[1] - pad_α, hole_α[2] + pad_α), (hole_β[1] - pad_β, hole_β[2] + pad_β)
end


# Picks n random pixels of the heatmap, without repeats, inside the sampling window.
# Pixels whose ray lands on the disc further out than r_max are left out: plot_paths_3d draws each
# path around its endpoint, so those would sit outside the 3D box.
function SamplePixels(α, β, points, n; sample_αlims, sample_βlims, r_max)
    inside(value, lims) = lims[1] <= value <= lims[2]
    lands_in_box(p) = p.status != StatusCodes.IntersectedWithGeometry || Gradus._equatorial_project(p.x) <= r_max

    # Rows of the image are β, columns are α
    candidates = [CartesianIndex(i, j) for i in eachindex(β), j in eachindex(α)
        if inside(α[j], sample_αlims) && inside(β[i], sample_βlims) && lands_in_box(points[i, j])]

    if isempty(candidates)
        error("No pixels to sample; check sample_αlims, sample_βlims and extent.")
    end
    return shuffle(candidates)[1:min(n, length(candidates))]
end


# Number of times a traced path passes through the equatorial plane (z = r cos θ changes sign).
# Landing on the disc doesn't count, since the trace stops just before the plane.
function PlaneCrossings(sol)
    height = [u[2] * cos(u[3]) for u in sol.u]
    return count(k -> height[k] * height[k+1] < 0, 1:length(height)-1)
end


# Listed in drawing order, so the rays that make the islands end up on top.
const CATEGORIES = [
    (label = "Hits disc directly",                     colour = :grey60,     linewidth = 0.6),
    (label = "Escapes to infinity",                    colour = :limegreen,  linewidth = 1.0),
    (label = "Reaches the inner boundary",             colour = :magenta,    linewidth = 1.0),
    (label = "Still orbiting at the end of the trace", colour = :red,        linewidth = 1.0),
    (label = "Hits disc after 1 plane crossing",       colour = :dodgerblue, linewidth = 1.6),
    (label = "Hits disc after 2+ plane crossings",     colour = :darkorange, linewidth = 1.6),
]

# Sorts a traced ray by where it ended up and how it got there. Returns an index into CATEGORIES.
function Categorise(sol, r_escape)
    status = Gradus.get_status_code(sol.prob.p)
    if status == StatusCodes.IntersectedWithGeometry
        crossings = PlaneCrossings(sol)
        return crossings == 0 ? 1 : crossings == 1 ? 5 : 6
    elseif status == StatusCodes.WithinInnerBoundary
        return 3
    elseif sol.u[end][2] > r_escape
        return 2 # ran out of integration time far from the hole, on its way out
    else
        return 4
    end
end


# plot_horizon_3d needs a horizon radius, which a > 1 doesn't have. This draws the surface the
# integrator stops rays at instead, in the same coordinates as plot_paths_3d. That is the inner radius
# of the integration chart (inner_radius(m) plus 1%); rays reaching it get the status WithinInnerBoundary.
function PlotInnerBoundary3D!(p, m; n_points = 32)
    R = Gradus.chart_for_metric(m).inner_radius

    if m.a > m.M
        # No horizon, so inner_radius is r ≈ 0: the ring singularity and the disc it spans.
        # plot_paths_3d treats (r, θ, ϕ) as spherical coordinates, where all of r = 0 is one point.
        scatter3d!(p, [0.0], [0.0], [0.0];
            color = :black,
            markersize = 5,
            label = "Ring singularity (r = 0)"
        )
    else
        u = range(0, 2π, length = n_points)
        v = range(0, π, length = n_points)
        surface!(p, R .* cos.(u) .* sin.(v)', R .* sin.(u) .* sin.(v)', R .* ones(n_points) .* cos.(v)';
            color = cgrad([:black, :black]),
            colorbar = false,
            label = false
        )
        plot!(p, [NaN], [NaN], [NaN]; color = :black, linewidth = 4, label = "Event horizon")
    end
    return p
end


# Disc drawn as rings in the equatorial plane. The inner edge is bold since that is the gap rays pass through.
function PlotDisc3D!(p, r_in, r_out, extent; rings = 6)
    ϕ = range(0, 2π, length = 200)
    r_edge = min(r_out, extent)

    for r in range(r_in, r_edge, length = rings)[2:end]
        plot!(p, r .* cos.(ϕ), r .* sin.(ϕ), zero(ϕ); color = :grey80, linewidth = 0.5, label = false)
    end
    plot!(p, r_in .* cos.(ϕ), r_in .* sin.(ϕ), zero(ϕ); color = :black, linewidth = 1.5, label = "Disc inner edge")
    return p
end


# An axis pointing straight at the camera has all of its ticks drawn on one spot, so it is left blank.
function EndOnAxis(camera)
    azimuth, elevation = camera
    if abs(elevation) == 90
        return (zformatter = _ -> "", zlabel = "")
    elseif elevation == 0 && mod(azimuth, 180) == 0
        return (yformatter = _ -> "", ylabel = "")
    elseif elevation == 0 && mod(azimuth, 180) == 90
        return (xformatter = _ -> "", xlabel = "")
    end
    return (;)
end


# One 3D view of the traced paths with the disc, the inner boundary and the direction to the observer.
# category[k] is the index into CATEGORIES of paths[k].
function PlotPaths3D(paths, category, m, r_in, r_out, angle; extent, t_span, n_points, camera, legend, kwargs...)
    colours = [CATEGORIES[c].colour for c in category]
    linewidths = [CATEGORIES[c].linewidth for c in category]

    p = plot_paths_3d(paths;
        extent = extent,
        t_span = t_span,
        n_points = n_points,
        color = permutedims(colours), # row vectors give one value per path
        linewidth = permutedims(linewidths),
        label = false,
        xlabel = "x",
        ylabel = "y",
        zlabel = "z",
        camera = camera,
        legend = legend,
        EndOnAxis(camera)...,
        kwargs...
    )

    # GR doesn't clip 3D lines to the axes, so the parts of each path outside the box are blanked
    for series in p.series_list
        outside = @. abs(series[:x]) > extent || abs(series[:y]) > extent || abs(series[:z]) > extent
        series[:x][outside] .= NaN
        series[:y][outside] .= NaN
        series[:z][outside] .= NaN
    end

    PlotDisc3D!(p, r_in, r_out, extent)

    # Direction the rays arrive from
    to_observer = 0.95 * extent .* (sin(angle), 0.0, cos(angle))
    plot!(p, [0.0, to_observer[1]], [0.0, to_observer[2]], [0.0, to_observer[3]];
        color = :black, linestyle = :dash, linewidth = 1.0, label = "To observer"
    )

    # Where each ray landed on the disc
    landed = [c == 1 || c == 5 || c == 6 for c in category]
    ends = [Gradus.to_cartesian(sol.u[end]) for sol in paths[landed]]
    scatter3d!(p, [e[1] for e in ends], [e[2] for e in ends], [e[3] for e in ends];
        color = colours[landed],
        markersize = 3,
        markerstrokewidth = 0.5,
        label = false
    )

    PlotInnerBoundary3D!(p, m)
    return p
end


"""
Redshift heatmap of a thin disc (as Heatmaps.jl) next to 3D traces of a random sample of the same rays.

Each sampled pixel is marked on the heatmap and its ray is drawn in 3D in the same colour, sorted by
where the ray ended up and how many times it passed through the equatorial plane on the way there.
Rays that reach the disc after passing through the plane are the ones that make the inner islands.

The islands are a small part of the image, so a plain random sample mostly misses them. Instead
n_candidates random pixels are traced and rays_per_kind of each kind of ray are drawn.
"""
function IslandTrace(a, size, angle;
    rays_per_kind = 12, # rays of each kind drawn in 3D
    n_candidates = 1500, # random pixels traced to find them
    r_in = nothing, # disc inner edge. nothing puts it at the ISCO
    r_out = 500.0, # disc outer edge
    αlims = (-10, 10),
    βlims = (-10, 10),
    sample_αlims = nothing, # window the rays are sampled from. nothing fits it to the hole in the image
    sample_βlims = nothing,
    color_limits = (0.5, 1.2),
    color = :inferno,
    extent = 8, # half-width of the 3D box
    t_span = 4 * extent, # affine length of path drawn up to each endpoint. Long enough to cover the whole way through the box
    n_points = 4000, # points drawn per path
    cameras = [(40, 25), (0, 0)], # (azimuth, elevation) of each 3D view. (0, 0) is side on, with the observer off to the right
    seed = 1 # fixed so the same rays come out each run
    )

    m = KerrMetric(M = 1.0, a = a)
    x = SVector(0.0, 1000.0, angle, 0.0)
    r_in = r_in === nothing ? Gradus.isco(m) : r_in
    d = ThinDisc(r_in, r_out)
    λ_max = 2 * x[2] # max integration time: must exceed the observer distance or rays stop before reaching the disc

    (; α, β, redshift, points, image_height) = RenderHeatmap(m, x, d, size; αlims = αlims, βlims = βlims, λ_max = λ_max)

    # ---- Sampling ----
    if sample_αlims === nothing || sample_βlims === nothing
        hole_α, hole_β = HoleWindow(α, β, points)
        sample_αlims = sample_αlims === nothing ? hole_α : sample_αlims
        sample_βlims = sample_βlims === nothing ? hole_β : sample_βlims
    end

    Random.seed!(seed)
    pixels = SamplePixels(α, β, points, n_candidates; sample_αlims = sample_αlims, sample_βlims = sample_βlims, r_max = extent)

    # ---- 3D trace ----
    # The function the render itself uses to turn a pixel number into a velocity, so these are the
    # same impact parameters as the heatmap, including the small offset Gradus adds to avoid α = 0.
    velocity = Gradus._render_velocity_function(m, x, size, image_height, αlims, βlims)
    pixel_number = LinearIndices(points)

    vs = [velocity(pixel_number[pixel]) for pixel in pixels]
    xs = fill(x, length(vs))
    sols = tracegeodesics(m, xs, vs, d, λ_max; save_on = true) # save_on keeps the whole path, not just the endpoint

    # Each traced ray should end exactly where its pixel in the heatmap did
    mismatched = count(k -> Gradus.get_status_code(sols.u[k].prob.p) != points[pixels[k]].status, eachindex(pixels))
    if mismatched > 0
        @warn "$mismatched of $(length(pixels)) traced rays ended differently to their heatmap pixel."
    end

    category = [Categorise(sol, x[2] / 2) for sol in sols.u]
    for c in eachindex(CATEGORIES)
        println(rpad(CATEGORIES[c].label, 42), count(==(c), category), " of $(length(pixels)) sampled rays")
    end

    # Keeping the first few rays of each kind. The pixels are in random order, so this is a random pick.
    # Going through the kinds in turn also leaves everything in drawing order.
    keep = reduce(vcat, [first(findall(==(c), category), rays_per_kind) for c in eachindex(CATEGORIES)])
    pixels, category, paths = pixels[keep], category[keep], sols.u[keep]

    style = (
        fontfamily = "Computer Modern",
        guidefontsize = 13,
        tickfontsize = 11,
        titlefontsize = 13,
        legendfontsize = 10
    )

    ϑ = round(rad2deg(angle), digits = 1)

    # ---- Heatmap ----
    p_map = heatmap(α, β, redshift;
        title = "Redshift Map, a = $a, ϑ = $(ϑ)°",
        color = color,
        clims = color_limits,
        colorbar_title = "Redshift (g)",
        xlabel = "α",
        ylabel = "β",
        aspect_ratio = 1,
        xlims = extrema(α),
        ylims = extrema(β),
        grid = false,
        framestyle = :box,
        legend = :bottomleft,
        style...
    )

    plot!(p_map,
        [sample_αlims[1], sample_αlims[2], sample_αlims[2], sample_αlims[1], sample_αlims[1]],
        [sample_βlims[1], sample_βlims[1], sample_βlims[2], sample_βlims[2], sample_βlims[1]];
        color = :white, linestyle = :dash, linewidth = 1.0, label = false
    )

    # One set of markers per kind of ray, only for the kinds that came up
    for c in eachindex(CATEGORIES)
        chosen = pixels[category .== c]
        if !isempty(chosen)
            scatter!(p_map, [α[pixel[2]] for pixel in chosen], [β[pixel[1]] for pixel in chosen];
                color = CATEGORIES[c].colour,
                markersize = 4,
                markerstrokecolor = :black,
                markerstrokewidth = 0.8,
                label = CATEGORIES[c].label
            )
        end
    end

    # ---- 3D paths ----
    p_3d = [PlotPaths3D(paths, category, m, r_in, r_out, angle;
        extent = extent,
        t_span = t_span,
        n_points = n_points,
        camera = camera,
        legend = i == 1 ? :topleft : false, # same in every view, so only shown once
        title = "Traced Paths, view = $camera",
        style...
        ) for (i, camera) in enumerate(cameras)]

    # The heatmap panel is wider than the 3D ones to leave room for its colourbar
    widths = [1020; fill(765, length(cameras))]
    combined = plot(p_map, p_3d...;
        layout = grid(1, length(widths), widths = widths ./ sum(widths)),
        size = (sum(widths), 850),
        dpi = 200,
        left_margin = 16Plots.mm,
        bottom_margin = 8Plots.mm
    )
    display(combined)
    savefig(combined, "IslandTrace_a$(a)_i$(round(Int, ϑ)).png")
    return combined
end


# size = 600 takes around 15 seconds each
IslandTrace(1.986, 600, deg2rad(80))
IslandTrace(1.986, 600, deg2rad(60))

#TODO: Make plotting clearer and focus sampling.
#TODO: Make each plot of the sample a seperate figure. 

