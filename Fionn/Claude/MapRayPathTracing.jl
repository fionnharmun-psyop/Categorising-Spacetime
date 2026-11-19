using Gradus, Plots, LaTeXStrings, Random
include("common.jl")

# Render -> Sampler -> Heatmap -> Tracer

function RenderDisc(spin, pos, imgsize; rmax = 30.0,
    αlims = (-5, 5), βlims = (-5, 5), λmax = 2000.0)

    m = KerrMetric(M = 1.0, a = spin)
    d = ThinDisc(Gradus.isco(m), rmax)

    image_height = round(Int, imgsize * (βlims[2] - βlims[1]) / (αlims[2] - αlims[1]))

    α, β, cache = prerendergeodesics(
        m, pos, d,
        λmax,
        αlims = αlims, βlims = βlims,
        image_width = imgsize, image_height = image_height,
        verbose = true
    )

    redshift = apply(ConstPointFunctions.redshift(m, pos) ∘ ConstPointFunctions.filter_intersected(), cache)
    radius = apply(PointFunction((m, gp, t) -> Gradus._equatorial_project(gp.x)) ∘ ConstPointFunctions.filter_intersected(), cache)
    return m, d, α, β, radius, redshift
end


function SampleRays(α, β, img; samples = nothing, num_samples = 20, rng = Random.default_rng())

    if isnothing(samples)
        pool = [CartesianIndex(i, j) for i in eachindex(β), j in eachindex(α)]
        
        isempty(pool) && error("No pixels to sample")

        chosen = shuffle(rng, pool)[1:min(num_samples, length(pool))]
        return [(α[I[2]], β[I[1]]) for I in chosen]
    else
        return [(Float64(a), Float64(b)) for (a, b) in samples]
    end
end

style = (
    fontfamily = "Computer Modern",
    framestyle = :box,
    grid = false,
    guidefontsize = 13,
    tickfontsize = 11,
    titlefontsize = 13
)

function HeatmappingSamples(α, β, img, samples; kwargs...)
    plot_w = 600                                                   
    canvas = (plot_w + 150, plot_w * length(β) / length(α) + 100)  
    dpi = 600

    p = heatmap(α, β, img;
        xlabel = L"\alpha",
        ylabel = L"\beta",
        aspect_ratio = 1,
        xlims = extrema(α),
        ylims = extrema(β),
        size = canvas,
        dpi = dpi,
        style...,
        kwargs...
    )

    scatter!(p, first.(samples), last.(samples);
        label = "Traced rays", markercolor = :white, markerstrokecolor = :black, markersize = 4)
    return p
end

function RenderTrace(m, pos, samples, disc; λmax = 2000.0, extent = 20, kwargs...)
    vs = [map_impact_parameters(m, pos, a, b) for (a, b) in samples]
    xs = fill(pos, length(vs))

    sols = tracegeodesics(m, xs, vs, disc, λmax)

    p = plot_paths_3d(sols; extent = extent, legend = false, kwargs...)
    return p, sols
end


spin = 1.998
θ = deg2rad(60)
pos = SVector(0.0, 1000.0, θ, 0.0)

m, d, α, β, radius, redshift = RenderDisc(spin, pos, 200)

samples = SampleRays(α, β, redshift; num_samples = 50)

p_map = HeatmappingSamples(α, β, redshift, samples;
    title = "a = $spin, ϑ = $(round(rad2deg(θ), digits = 1))°",
    colorbar_title = L"Redshift ($g$)"
)
display(p_map)

p_paths, sols = RenderTrace(m, pos, samples, d)
display(p_paths)
