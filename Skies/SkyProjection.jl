using Gradus, Plots

# This is all largely copies straight from the tutorial apart from the patch at the bottom

# Changing definition of inner_radius to not throw and error when a is greater than 1.
function Gradus.inner_radius(m::Gradus.KerrMetric{T}) where {T}     
    disc = m.M^2 - m.a^2 # Operand of the sqrt
    disc < 0 ? 4eps(T) : m.M + √disc # "If sqrt operand is less than 0 then return as close to zero as possible but not zero
end

m = KerrMetric(1.0, 1.0)
# the ring-corona can also be used to represent any off-axis point
corona = RingCorona(; r = 3.0, h = 3.0)
d = ThinDisc(0.0, Inf)

sky = AdaptiveSky(m, corona, d)

parameter_label(m, corona) =
    "Spin a = $(m.a) | Corona height h = $(corona.h) | Corona radius r = $(corona.r)"

# trace an initial grid of points
trace_initial!(sky)

function plot_sky(sky; kwargs...)
    phi, cos_theta, values = unpack_sky(sky)
    theta = acos.(cos_theta)
    # get the radius that each geodesic hit
    Rs = [v.r for v in values]

    fig, ax, _ = scatter(phi, theta; color = log10.(Rs), kwargs...)
    ax.xlabel = "ϕ"
    ax.ylabel = "θ"
    ax.title = parameter_label(m, corona)
    fig
end

function refine_edges(sky::AdaptiveSky, i1::Int, i2::Int)
    v1 = sky.values[i1]
    v2 = sky.values[i2]
    if isnan(v1.r) && isnan(v2.r)
        false # do not refine
    else
        isnan(v1.r) || isnan(v2.r)
    end
end

trace_step!(sky; check_refine = refine_edges, verbose = true)

for i in 1:2
    trace_step!(sky; check_refine = refine_edges, verbose = true)
end

function refine_distant(sky::AdaptiveSky, i1::Int, i2::Int)
    v1 = sky.values[i1]
    v2 = sky.values[i2]
    if isnan(v1.r) && isnan(v2.r)
        false # do not refine
    else
        # refine if they are more than 10% different
        !isapprox(v1.r, v2.r, rtol = 0.1)
    end
end

# create a fresh sky
sky = AdaptiveSky(m, corona, d)
trace_initial!(sky)
for i in 1:3
    trace_step!(sky; check_refine = refine_distant, verbose = true)
end

# Patch
# Work around Gradus constructing a six-element vector from only five fields.
Gradus.vector_average(
    weights::Gradus.SVector{N,T},
    values::Gradus.SVector{N,V},
) where {N,T<:Number,V<:Gradus.CoronaGridValues} = V(
    sum(weights[i] * values[i].t for i in eachindex(weights)),
    sum(weights[i] * values[i].r for i in eachindex(weights)),
    sum(weights[i] * values[i].ϕ for i in eachindex(weights)),
    sum(weights[i] * values[i].g for i in eachindex(weights)),
    sum(weights[i] * values[i].J for i in eachindex(weights)),
    sum(weights[i] * values[i].status for i in eachindex(weights)),
)

phi, theta, grid = fill_sky_values(sky, 1080)
Rs = [v.r for v in grid]

p  = heatmap(phi, theta, log10.(Rs)')

display(p)