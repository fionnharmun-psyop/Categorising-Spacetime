using Gradus, Plots

# --- patch: let inner_radius degrade gracefully instead of throwing when
#     there's no real horizon (a > M). KerrMetric has no constructor guard,
#     so this is the only change needed to build and trace it.
function Gradus.inner_radius(m::Gradus.KerrMetric{T}) where {T}
    disc = m.M^2 - m.a^2
    disc < 0 ? 4eps(T) : m.M + √disc
end

m = KerrMetric(M = 1.0, a = 1.5)

x = SVector(0.0, 1000.0, deg2rad(80.0), 0.0)

# accretion disc -- inner edge pinned at the ISCO, computed via the same
# closed-form formula the redshift pf itself uses, so the plunging-region
# formulas are never evaluated arbitrarily close to the true r=0 singularity
isco = Gradus.isco(m)
d = ThinDisc(isco, 30.0)

pf = ConstPointFunctions.redshift(m, x) ∘ ConstPointFunctions.filter_intersected()

α, β, img = rendergeodesics(
    m, x, d,
    20000.0,
    αlims = (-10, 10),
    βlims = (-5, 5),
    image_width = 800,
    image_height = 800,
    verbose = true,
    pf = pf,
)

heatmap(α, β, img;
    xlabel = "α",
    ylabel = "β",
    colorbar_title = "g",
    title = "Kerr naked singularity, a = $(m.a)M, θ_obs = $(round(rad2deg(x[3]), digits=1))°"

)