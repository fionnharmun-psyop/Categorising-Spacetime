using Gradus, Plots, Statistics

# a² > M² has no horizon, so return a tiny inner radius instead of NaN
function Gradus.inner_radius(m::Gradus.KerrMetric{T}) where {T}
    disc = m.M^2 - m.a^2
    disc < 0 ? 4eps(T) : m.M + √disc
end

function render_disc(a;
        θ = 80.0,
        emis_exp = 2.0,            # emission ∝ r^-emis_exp (was 3.0)
        r_floor = 6.0,             # radii below this are treated as r_floor in the emission profile
        gamma = 0.45,              # was 0.5
        black_point = 0.03,
        norm_ref = nothing,        # reuse another render's norm for identical brightness scaling
        αlims = (-34, 34),
        βlims = (-12, 16),
        width = 2400,              # pixels along α; height follows from the limits
        outfile = "render_a$(a)")

    bg = RGB(0.12, 0.12, 0.12)
    height = round(Int, width * (βlims[2] - βlims[1]) / (αlims[2] - αlims[1]))

    m  = KerrMetric(M = 1.0, a = a)
    xo = SVector(0.0, 1000.0, deg2rad(θ), 0.0)
    isco = Gradus.isco(m)
    @show a isco width height
    isfinite(isco) || error("ISCO is not finite for a = $a")
    d  = ThinDisc(isco, 30.0)

    kw = (; αlims, βlims, image_width = width, image_height = height, verbose = true)

    # g image
    pf_g = ConstPointFunctions.redshift(m, xo) ∘ ConstPointFunctions.filter_intersected()
    α, β, g_img = rendergeodesics(m, xo, d, 2000.0; pf = pf_g, kw...)

    # radius image
    pf_r = PointFunction((m, gp, max_time; kwargs...) -> gp.x[2]) ∘
           ConstPointFunctions.filter_intersected()
    _, _, r_img = rendergeodesics(m, xo, d, 2000.0; pf = pf_r, kw...)

    # diagnostics
    @show extrema(filter(!isnan, g_img))
    rr = filter(!isnan, r_img)
    @show extrema(rr) quantile(rr, [0.01, 0.1, 0.5])

    # tone mapping from this call's own images
    g = replace(g_img, NaN => 0.0)
    r = replace(r_img, NaN => 1.0)

    I = g .^ 3 .* max.(r, r_floor) .^ -emis_exp     # g³ (was g⁴): gentler on the receding side
    nrm = isnothing(norm_ref) ? quantile(I[I .> 0], 0.99) : norm_ref
    I = clamp.(I ./ nrm, 0, 1)
    I = clamp.((I .- black_point) ./ (1 - black_point), 0, 1) .^ gamma

    p = heatmap(α, β, I;
        c = cgrad([bg, RGB(0.40, 0.35, 0.30), RGB(1.0, 0.96, 0.88)]),
        clims = (0, 1),
        colorbar = false,
        background_color = bg,
        titlefontcolor = :white,
        title = "Kerr Blackhole, a = $a, ϑ = $(round(θ, digits = 1))°",
        aspect_ratio = 1,
        fontfamily = "Computer Modern",
        framestyle = :none,
        grid = false,
        size = (width ÷ 2, height ÷ 2),
        dpi = 200,
    )

    savefig(p, outfile * ".png")
    return (; p, g_img, r_img, norm = nrm, α, β)
end

# ---------------- run ----------------
common = (; θ = 80.0, αlims = (-34, 34), βlims = (-12, 16), width = 2400)   # raise to 2400 for the final render

a0 = render_disc(0.0; common..., outfile = "full_a0.0")
display(a0.p)

nks = render_disc(1.968; common..., norm_ref = a0.norm, outfile = "full_a1.968")
display(nks.p)