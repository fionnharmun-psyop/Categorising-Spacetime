using Gradus, Plots
include("common.jl")

function RenderPlot(a, size, angle; color_limits = (0.3, 1.3), color = :inferno, αlims = (-5, 5), βlims = (-5, 5))

    m = KerrMetric(M = 1.0, a = a)
    x = SVector(0.0, 1000.0, angle, 0.0)
    isco = Gradus.isco(m)
    d = ThinDisc(isco, 30.0)
    pf = ConstPointFunctions.redshift(m, x) ∘ ConstPointFunctions.filter_intersected()

    image_height = round(Int, size * (βlims[2] - βlims[1]) / (αlims[2] - αlims[1]))

    α, β, img = rendergeodesics(
    m, x, d,
    2000.0,
    αlims = αlims,
    βlims = βlims,
    image_width = size,
    image_height = image_height,
    verbose = true,
    pf = pf
    )

    # Canvas fixed in points with zero margins so the image fills the frame; dpi sets the pixel count
    # (output pixels = canvas × dpi / 100), here 4 output pixels per traced ray.
    canvas_w = 600
    p = heatmap(α, β, img;
    color = color,
    clims = color_limits,
    colorbar = false,
    background_colour = :ivory,
    grid = false,
    framestyle = :none,
    ticks = false,
    margin = 0Plots.mm,
    aspect_ratio = 1,
    xlims = extrema(α),
    ylims = extrema(β),
    size = (canvas_w, canvas_w * image_height / size),
    dpi = round(Int, 100 * 4 * size / canvas_w)
    )

    any(isfinite, img) || @warn "No rays hit the disc at ϑ = $(round(rad2deg(angle), digits = 1))°; frame will be blank."
    return p
end


fps = 24
runtime = 4
angles = range(1.0, 89.0, fps * runtime)

n_frames = length(angles)
anim = Animation()
for (i, theta) in enumerate(angles)
    @info "Rendering frame $i of $n_frames" theta
    # Put it in a try loop so it doesn't crash while i'm away
    try
        frame(anim, RenderPlot(1.002, 800, deg2rad(theta)))
    catch err
        @error "Frame $i failed; skipping it" exception = err # a bad frame shouldn't kill an unattended run
    end
end
mp4(anim, "a1.002-inclination.mp4", fps = fps)  

spins = range(0.0, 1.9, fps * runtime)
anim2 = Animation()
for (i, a) in enumerate(spins)
    @info "Rendering frame $i of $n_frames" 
    try
        frame(anim2, RenderPlot(a, 800, deg2rad(80); αlims = (-12, 12), βlims = (-12, 12)))
    catch err
        @error "Frame $i failed; skipping it" exception = err
    end
end
mp4(anim2, "theta80-spins.mp4", fps = fps)
