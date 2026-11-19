using Gradus, Plots, LaTeXStrings
include("common.jl")

# Which photons pass through the r = 0 disc of a superspinning Kerr metric into r < 0, which of
# those turn round and come back, and where each one ends up.
#
# All ray tracing and all constants of motion come from Gradus. The only physics added here is the
# radial and polar potentials and what follows from them. Sources, by number, in the comments below:
#   [1] Carter (1968), Phys. Rev. 174, 1559. Constants of motion, the r < 0 extension, causality violation.
#   [2] Bardeen, Press & Teukolsky (1972), ApJ 178, 347. R(r) and Θ(θ) in the form used here.
#   [3] Bardeen (1973), in "Black Holes" (Les Houches), eds. DeWitt & DeWitt. Impact parameters α, β.
#   [4] Chandrasekhar (1983), "The Mathematical Theory of Black Holes", ch. 7. Null geodesics with negative Carter constant.
#   [5] Visser (2007), arXiv:0706.0622. Kerr-Schild coordinates and the r = 0 disc.
# The references are given from memory: check them (and find equation numbers) before citing.
# Each derivation is written out in full so it can be checked without them. M = 1 throughout.


# ---------------- Physics ----------------

# Constants of motion of a photon, from its Gradus position and velocity.
#   E, L : energy and axial angular momentum (Gradus.E, Gradus.Lz)
#   Q    : Carter constant, Q = p_θ² + cos²θ (L² / sin²θ - a² E²) for a photon [1, 2]
#   l, q : L / E and Q / E². A photon's path only depends on these two ratios, not on its energy
#   k    : q + (l - a)², the other common form of the Carter constant. Never negative (see PolarPotential)
function PhotonConstants(m, x, v)
    velocity = Gradus.constrain_all(m, x, v, 0.0)
    E = Gradus.E(m, x, velocity)
    L = Gradus.Lz(m, x, velocity)
    g_θθ = Gradus.metric_components(m, SVector(x[2], x[3]))[3]
    Q = (g_θθ * velocity[3])^2 + cos(x[3])^2 * (L^2 / sin(x[3])^2 - m.a^2 * E^2)

    l, q = L / E, Q / E^2
    k = q + (l - m.a)^2
    return (; E, L, Q, l, q, k)
end

# Radial potential of a photon [2]. The radial motion obeys Σ² (dr/dλ)² = E² R(r), so the photon can
# only be where R ≥ 0 and it turns round where R = 0.
#   R(r) = (r² + a² - a l)² - (r² - 2r + a²) (q + (l - a)²)
# Multiplying out gives
#   R(r) = r⁴ + (a² - l² - q) r² + 2 k r - a² q
# Three things follow from this:
#   1. R(0) = -a² q. A photon can only reach r = 0 if q ≤ 0. With q > 0 it always turns at some r > 0.
#   2. If q < 0 every term is positive for r > 0 (a² - l² > 0 is needed for the photon to exist, see
#      PhotonAllowed). So there is no turning point outside: every q < 0 photon goes through the disc.
#   3. For r < 0 the term 2 k r is negative, and it is the only negative term. If k is large enough
#      it pulls R down to zero and the photon turns round inside r < 0. If k is small it never does,
#      and the photon carries on to r → -∞. k = 0 gives R = (r² + a² - a l)², which has no turning point at all.
function RadialPotential(r, a, l, q)
    return (r^2 + a^2 - a * l)^2 - (r^2 - 2r + a^2) * (q + (l - a)^2)
end

# Polar potential of a photon [2]. The polar motion obeys Σ² (dθ/dλ)² = E² Θ(θ), so Θ ≥ 0 is needed.
#   Θ(θ) = q + cos²θ (a² - l² / sin²θ)
# On the equator Θ = q, so a photon that ever touches the equatorial plane has q ≥ 0. A q < 0 photon
# stays in one hemisphere, which is why nothing emitted by a thin equatorial disc can enter r < 0 [4].
# Also Θ + (a sinθ - l / sinθ)² = q + (l - a)² = k, so k ≥ 0 always.
function PolarPotential(θ, a, l, q)
    return q + cos(θ)^2 * (a^2 - l^2 / sin(θ)^2)
end

# Whether a photon with these constants can exist at all, i.e. whether Θ ≥ 0 somewhere.
# For q ≥ 0 it always can. For q < 0, maximising Θ over θ (the maximum is at sin²θ = |l| / a)
# gives Θ_max = q + (a - |l|)², so the condition is |l| < a and q ≥ -(a - |l|)².
function PhotonAllowed(a, l, q)
    q >= 0 && return true
    return abs(l) < a && q >= -(a - abs(l))^2
end

# Bisection for a root of f between lo and hi, where f changes sign.
function Bisect(f, lo, hi; n_steps = 60)
    f_lo = f(lo)
    for i in 1:n_steps
        mid = 0.5 * (lo + hi)
        if (f(mid) > 0) == (f_lo > 0)
            lo = mid
        else
            hi = mid
        end
    end
    return 0.5 * (lo + hi)
end

# Radius where the photon turns round, found from the first place R(r) drops to zero along its path.
#   q < 0 : the photon comes through the disc, so search from r = 0 into negative r. NaN means no turning point.
#   q ≥ 0 : the photon comes in from far away, so search inwards from large positive r.
function TurningPoint(a, l, q; r_search = 100.0, n_scan = 4000)
    # Distances from r = 0, log spaced since a rebound can be very shallow
    distances = exp.(range(log(1e-6), log(r_search), length = n_scan))
    radii = q < 0 ? -distances : reverse(distances)

    for i in 1:n_scan-1
        R_here = RadialPotential(radii[i], a, l, q)
        R_next = RadialPotential(radii[i+1], a, l, q)
        if R_here > 0 && R_next <= 0
            return Bisect(r -> RadialPotential(r, a, l, q), radii[i], radii[i+1])
        end
    end
    return NaN
end

# What the potentials say a photon coming in from far away will do.
#   :forbidden : no photon has these constants
#   :outside   : q ≥ 0, turns round at r > 0 and never reaches the disc
#   :rebound   : q < 0, passes through the disc, turns round in r < 0 and comes back out
#   :through   : q < 0, passes through the disc and carries on to r → -∞
function PredictedFate(a, l, q)
    PhotonAllowed(a, l, q) || return :forbidden
    q >= 0 && return :outside
    return isnan(TurningPoint(a, l, q)) ? :through : :rebound
end


# ---------------- Plot helpers ----------------

FATES = [:forbidden, :outside, :rebound, :through, :disc]
FATE_LABELS = ["no such photon", "turns at r > 0", "rebounds in r < 0", "through to r → -∞", "hits the disc"]
FATE_COLOURS = [:gray85, :lightsteelblue, :darkorange, :black, :gold]

FateCode(fate) = findfirst(==(fate), FATES)

style = (
    fontfamily = "Computer Modern",
    framestyle = :box,
    grid = false,
    guidefontsize = 13,
    tickfontsize = 11,
    titlefontsize = 13
)

# Heatmap of fate codes with a legend instead of a colourbar. `shown` picks which fates go in the legend
function FateHeatmap(xs, ys, codes; shown = FATES, kwargs...)
    p = heatmap(xs, ys, codes;
        c = cgrad(FATE_COLOURS, categorical = true),
        clims = (1, length(FATES)),
        colorbar = false,
        style...,
        kwargs...
    )
    for fate in shown
        i = FateCode(fate)
        scatter!(p, [NaN], [NaN]; markershape = :square, markersize = 6, color = FATE_COLOURS[i], label = FATE_LABELS[i])
    end
    return p
end

function jagged_vline!(p, x0, ymin, ymax; n = 20, A = 0.2, color = :black)
    ys = range(ymin, ymax; length = 2n + 1)
    xs = [x0 + (isodd(i) ? -A : A) for i in eachindex(ys)]
    plot!(p, xs, ys; color = color, linewidth = 2.0, label = false)
end


# ---------------- 1. Which constants give which fate ----------------

# Map of the (l, q) plane, from the potentials alone. No ray tracing.
function PlotConstantsMap(spin; n = 300)
    ls = range(-1.3 * spin, 1.3 * spin, length = n)
    qs = range(-1.1 * spin^2, 0.5 * spin^2, length = n)

    # Rows are q, columns are l
    codes = [FateCode(PredictedFate(spin, l, q)) for q in qs, l in ls]

    p = FateHeatmap(ls, qs, codes;
        shown = [:forbidden, :outside, :rebound, :through],
        xlabel = L"l = L / E",
        ylabel = L"q = Q / E^2",
        title = "Fate of an incoming photon, a = $spin",
        legend = :bottomright,
        size = (800, 650),
        dpi = 300
    )
    hline!(p, [0.0]; color = :white, linestyle = :dash, label = false)

    display(p)
    savefig(p, "NegativeR_constants_a$(spin).png")
    return p
end


# ---------------- 2. Where on the image plane those photons are ----------------

# Traces a grid of rays back from the observer with Gradus and records where each one ends.
# The chart is extended to negative r so rays can pass through the disc instead of stopping at r = 0.
function TraceEndpoints(spin, angle; αlims = (-3, 3), βlims = (-3, 3), n = 150,
    r_observer = 1000.0, r_escape = 50.0, with_disc = true, r_out = 30.0)
    m = KerrMetric(M = 1.0, a = spin)
    x = SVector(0.0, r_observer, deg2rad(angle), 0.0)
    chart = Gradus.PolarChart(-r_escape, 1.1 * r_observer)

    α = range(αlims[1], αlims[2], length = n)
    β = range(βlims[1], βlims[2], length = n)

    # Rows are β, columns are α, to match heatmap
    vs = vec([map_impact_parameters(m, x, a, b) for b in β, a in α])
    xs = fill(x, length(vs))

    if with_disc
        d = ThinDisc(Gradus.isco(m), r_out)
        sols = tracegeodesics(m, xs, vs, d, 3 * r_observer; chart = chart, save_on = false)
    else
        sols = tracegeodesics(m, xs, vs, 3 * r_observer; chart = chart, save_on = false)
    end
    points = reshape(unpack_solution.(sols.u), (n, n))
    constants = reshape([PhotonConstants(m, x, v) for v in vs], (n, n))

    return (; m, α, β, points, constants, angle)
end

# Image plane coloured by where each ray ends, next to how far into r < 0 the rebounding rays get.
# The traced fate comes from the Gradus status code of each ray:
#   IntersectedWithGeometry : hit the disc
#   WithinInnerBoundary     : left the chart at r = -r_escape, so went through to negative r
#   OutOfDomain             : left the chart at large positive r. With q < 0 it must have been through the
#                             disc first (point 2 under RadialPotential), so that is a rebound.
# The rebound / turns-outside split of the OutOfDomain rays therefore uses q, not the path itself.
function PlotEndpointMap(spin, angle; kwargs...)
    (; m, α, β, points, constants) = TraceEndpoints(spin, angle; kwargs...)
    a = m.a

    traced = map(points, constants) do gp, c
        if gp.status == StatusCodes.IntersectedWithGeometry
            FateCode(:disc)
        elseif gp.status == StatusCodes.WithinInnerBoundary
            FateCode(:through)
        elseif gp.status == StatusCodes.OutOfDomain
            c.q < 0 ? FateCode(:rebound) : FateCode(:outside)
        else
            NaN # still travelling when the integration stopped: raise the max time
        end
    end

    # What the potentials predict for the same rays, from the Gradus constants of each one
    predicted = [FateCode(PredictedFate(a, c.l, c.q)) for c in constants]
    turning = [c.q < 0 ? TurningPoint(a, c.l, c.q) : NaN for c in constants]

    # Checking the trace against the prediction, leaving out rays the disc caught
    compared = (traced .!= FateCode(:disc)) .& isfinite.(traced)
    agree = count(traced[compared] .== predicted[compared])
    println("a = $spin, ϑ = $(angle)°: traced fate matches the potential for $agree of $(count(compared)) rays")

    p_fate = FateHeatmap(α, β, traced;
        shown = [:outside, :rebound, :through, :disc],
        xlabel = L"\alpha",
        ylabel = L"\beta",
        title = "Where each ray ends, a = $spin, ϑ = $(angle)°",
        aspect_ratio = 1,
        xlims = extrema(α),
        ylims = extrema(β),
        legend = :topright
    )

    # Edge of the q < 0 region for an observer far away [3]:
    #   q = β² + cos²ϑ (α² - a²), so q = 0 is the ellipse α² + β² / cos²ϑ = a².
    # Only rays inside it can reach the disc. It does not depend on the sign convention of α.
    φ = range(0, 2π, length = 400)
    plot!(p_fate, a .* cos.(φ), a .* cos(deg2rad(angle)) .* sin.(φ);
        color = :red,
        linestyle = :dash,
        linewidth = 1.5,
        label = "q = 0"
    )

    p_turn = heatmap(α, β, turning;
        c = :viridis,
        colorbar_title = "Turning radius",
        xlabel = L"\alpha",
        ylabel = L"\beta",
        title = "How far into r < 0 the rebounding rays get",
        aspect_ratio = 1,
        xlims = extrema(α),
        ylims = extrema(β),
        style...
    )

    combined = plot(p_fate, p_turn; layout = (1, 2), size = (1500, 650), dpi = 300)
    display(combined)
    savefig(combined, "NegativeR_endpoints_a$(spin)_i$(angle).png")
    return combined
end


# ---------------- 3. One ray in detail ----------------

# Traces one photon and shows its spacetime diagram, its radial potential and its path in space,
# with the turning point predicted from the potential marked on each for comparison with the trace.
function RenderNegativeRTrace(spin, x::SVector, v::SVector; λ_max = 3 * x[2] + 100.0,
    r_escape = 30.0, r_window = 6.0, name = "ray")

    # Metric and initial conditions
    m = KerrMetric(M = 1.0, a = spin)
    a = m.a
    v_list = [v]
    xs = fill(x, length(v_list))

    # Geodesic solution, with the chart extended to negative r
    chart = Gradus.PolarChart(-r_escape, max(1.1 * x[2], 2 * r_window))
    sols = tracegeodesics(m, xs, v_list, (0.0, λ_max); μ = 0.0, chart = chart)
    path = sols.u[1]

    t = [u[1] for u in path.u]
    r = [u[2] for u in path.u]
    θ = [u[3] for u in path.u]

    # Constants of motion and what the potential predicts
    (; l, q, k) = PhotonConstants(m, x, v)
    fate = PredictedFate(a, l, q)
    r_turn = TurningPoint(a, l, q)

    println("---- $name ----")
    println("l = $(round(l, digits = 4)), q = $(round(q, digits = 4)), k = $(round(k, digits = 4))")
    println("Predicted: $(FATE_LABELS[FateCode(fate)]), turning radius = $(round(r_turn, digits = 4))")
    println("Traced:    smallest r reached = $(round(minimum(r), digits = 4)), final r = $(round(r[end], digits = 2))")
    # Near the ring on the r < 0 side is where the Kerr metric allows closed timelike curves [1],
    # so this is the place to look for coordinate time running backwards along a ray.
    println("Coordinate time decreases somewhere along the path: $(any(diff(t) .< 0))")

    # Only the part of the path near the hole is plotted, with t = 0 where the ray enters the window
    start = findfirst(<(r_window), abs.(r))
    start === nothing && (start = 1)
    t = t .- t[start]
    near = abs.(r) .< r_window

    # ---- Spacetime diagram ----
    t_lo, t_hi = extrema(t[near])
    p_time = plot(r[near], t[near];
        xlabel = "r",
        ylabel = "t",
        title = "Spacetime Diagram, a = $a",
        xlims = (-r_window, r_window),
        color = :black,
        linewidth = 2.0,
        label = "Light",
        legend = :topright,
        style...,
        framestyle = :zerolines
    )
    jagged_vline!(p_time, 0.0, t_lo, t_hi; color = :black)
    if isfinite(r_turn)
        vline!(p_time, [r_turn]; color = :darkorange, linestyle = :dash, label = "Predicted turning point")
    end

    # ---- Radial potential ----
    # Divided by (r² + a²)² so it stays of order one: it tends to 1 far away and equals -q / a² at r = 0.
    # The sign and the zeros are unchanged, and those are all that matter.
    rs = range(-r_window, r_window; length = 2000)
    p_potential = plot(rs, RadialPotential.(rs, a, l, q) ./ (rs .^ 2 .+ a^2) .^ 2;
        xlabel = "r",
        ylabel = L"R(r) / (r^2 + a^2)^2",
        title = "Photon Radial Potential, l = $(round(l, digits = 3)), q = $(round(q, digits = 3))",
        color = :black,
        linewidth = 2.0,
        label = false,
        legend = :topright,
        style...
    )
    hline!(p_potential, [0.0]; linestyle = :dash, color = :gray, label = false)
    vline!(p_potential, [0.0]; linestyle = :dot, color = :gray, label = false)
    if isfinite(r_turn)
        scatter!(p_potential, [r_turn], [0.0]; color = :darkorange, markersize = 6, label = "Turning point")
    end

    # ---- Path in space ----
    # Kerr-Schild coordinates [5]: distance from the spin axis ρ = sqrt(r² + a²) sinθ, height z = r cosθ.
    # r = 0 is the flat disc z = 0, ρ < a, and the singularity is only its rim ρ = a.
    # The r < 0 sheet covers the same (ρ, z) points again, so that part of the path is drawn dashed.
    ρ = sqrt.(r .^ 2 .+ a^2) .* sin.(θ)
    z = r .* cos.(θ)
    outside = near .& (r .>= 0)
    inside = near .& (r .< 0)

    p_space = plot([0.0, a], [0.0, 0.0];
        color = :gray,
        linewidth = 4.0,
        label = "r = 0 disc",
        xlabel = L"\rho",
        ylabel = L"z",
        title = "Path in Space",
        aspect_ratio = 1,
        legend = :topright,
        style...
    )
    scatter!(p_space, [a], [0.0]; color = :red, markersize = 6, label = "Ring singularity")
    plot!(p_space, ifelse.(outside, ρ, NaN), ifelse.(outside, z, NaN); color = :black, linewidth = 2.0, label = "r > 0")
    plot!(p_space, ifelse.(inside, ρ, NaN), ifelse.(inside, z, NaN); color = :darkorange, linewidth = 2.0, linestyle = :dash, label = "r < 0")

    combined = plot(p_time, p_potential, p_space; layout = (1, 3), size = (1900, 600), dpi = 300)
    display(combined)
    savefig(combined, "NegativeR_trace_a$(spin)_$(name).png")
    return combined
end

# Same, but choosing the ray by where it lands on a distant observer's image plane.
function RenderNegativeRTrace(spin, angle::Real, α::Real, β::Real; r_observer = 1000.0, kwargs...)
    m = KerrMetric(M = 1.0, a = spin)
    x = SVector(0.0, r_observer, deg2rad(angle), 0.0)
    v = map_impact_parameters(m, x, α, β)
    RenderNegativeRTrace(spin, x, v; name = "i$(angle)_alpha$(α)_beta$(β)", kwargs...)
end


# ---------------- Run ----------------

# 1. Which (l, q) rebound, which go through. Instant
PlotConstantsMap(1.968)

# 2. Where those rays sit on the image plane. n = 150 is 22,500 rays
PlotEndpointMap(1.968, 60; n = 150)
PlotEndpointMap(1.968, 5; n = 150)

# 3. Single rays at ϑ = 60°. The expected fates are from q = β² + cos²ϑ (α² - a²) worked by hand,
#    so go by what each call prints
RenderNegativeRTrace(1.968, 60, 0.0, 0.75) # q just below zero: should rebound
RenderNegativeRTrace(1.968, 60, 0.0, 0.3)  # q well below zero: should go through
RenderNegativeRTrace(1.968, 60, 0.0, 1.2)  # q above zero: should turn at r > 0

# The ray from the old version of this script: straight down the spin axis.
# Here l ≈ 0 and q ≈ -a², so k ≈ 0 and R ≈ (r² + a²)², which is why it showed no turning point
v = SVector(0.0, -1.0, eps(), 0.0)
x = SVector(0.0, 3.0, deg2rad(1), 0.0)
RenderNegativeRTrace(1.968, x, v; name = "on_axis")

# TODO: put the endpoint condition as a pointfunction on rendergeodesics. 