using Gradus, Plots, OrdinaryDiffEq

# METRIC
m = KerrMetric(M = 1.0, a = 0.9)

# Changing definition of inner_radius to not throw and error when a is greater than 1.
function Gradus.inner_radius(m::Gradus.KerrMetric{T}) where {T}
    disc = m.M^2 - m.a^2 # Operand of the sqrt
    disc < 0 ? 4eps(T) : m.M + √disc # "If sqrt operand is less than 0 then return as close to zero as possible but not zero
end

# INITIAL CONDITIONS
θ = π/2  
r0 = 4.0

# LIGHT RAY ODE
function lightray(u, p, t)
    metric, θ, direction = p
    r = u
    gtt, grr = Gradus.metric_components(metric, SVector(r, θ))
    return direction * sqrt(-gtt / grr) # Condition for null-rays where the ST interval = 0
end

# SOLVING
function solve_ray(direction, tspan)
    solve(ODEProblem(lightray, r0, tspan, (m, θ, direction)), Tsit5(); reltol = 1e-10, abstol = 1e-10)
end

# PLOTTING
default(
    fontfamily="Computer Modern",
    linewidth=2.0,
    framestyle=:zerolines,
    label=nothing,
    grid=false,
    legend=false,
)

Tmax = 10.0
n = 400

future_times = collect(range(0.0, Tmax; length=n))
past_times = collect(range(0.0, -Tmax; length=n))

future_out = solve_ray(1, (0.0, Tmax))
future_in  = solve_ray(-1, (0.0, Tmax))
past_out   = solve_ray(1, (0.0, -Tmax))
past_in    = solve_ray(-1, (0.0, -Tmax))

p = plot(xlabel="r", ylabel="t", title = "Minkowski Diagram of a Schwarzchild Blackhole")

# AI did the shading of the cones so I only vaguely get how it works.
# Shade the regions between the two future/past light-ray boundaries.
plot!(p, Shape(
    vcat([future_out(t) for t in future_times],
         reverse([future_in(t) for t in future_times])),
    vcat(future_times, reverse(future_times)),
); color=:deepskyblue, fillalpha=0.25, linealpha=0, label=false)

plot!(p, Shape(
    vcat([past_out(t) for t in past_times],
         reverse([past_in(t) for t in past_times])),
    vcat(past_times, reverse(past_times)),
); color=:orange, fillalpha=0.25, linealpha=0, label=false)

# Draw the light-ray boundaries on top of the shading.
for ray in (future_out, future_in, past_out, past_in)
    plot!(p, ray.u, ray.t; color=:black, label=false)
end

println("All Working")
vline!(p, [Gradus.inner_radius(m)]; linestyle=:dash, color=:red, label=false)
display(p)