using Gradus, Plots, OrdinaryDiffEq

function Gradus.inner_radius(m::Gradus.KerrMetric{T}) where {T}
    disc = m.M^2 - m.a^2 # Operand of the sqrt
    disc < 0 ? 4eps(T) : m.M + √disc # "If sqrt operand is less than 0 then return as close to zero as possible but not zero
end

m = KerrMetric(M=1.0, a=1.2)

# INITIAL CONDITIONS
θ = π / 2
r0 = 3.0
t0 = 0.0
Tmax = 10.0
E = 1.0
Lz = 0.0 
horizon = Gradus.inner_radius(m)

# SET UP ODE FUNCTION FOR NULL RAYS
function lightray(r, p, t)
    metric, θ, energy, angular_momentum, radial_direction = p
    gtt, grr, _, gphiphi, gtphi =
        Gradus.metric_components(metric, SVector(r, θ))

    determinant = gtt * gphiphi - gtphi^2
    gtt_inverse = gphiphi / determinant
    gtphi_inverse = -gtphi / determinant
    gphiphi_inverse = gtt / determinant

    # For p_t = -E, p_phi = Lz, and equatorial p_theta = 0,
    # solve g^mu_nu p_mu p_nu = 0 for the radial momentum.
    radial_momentum_squared =
        -grr * (gtt_inverse * energy^2 -
                2gtphi_inverse * energy * angular_momentum +
                gphiphi_inverse * angular_momentum^2)
    radial_momentum_squared >= 0 ||
        throw(DomainError(radial_momentum_squared,
            "The selected null geodesic has reached a radial turning point."))

    radial_momentum = radial_direction * sqrt(radial_momentum_squared)
    dt_dλ = -gtt_inverse * energy + gtphi_inverse * angular_momentum
    dr_dλ = radial_momentum / grr
    return dr_dλ / dt_dλ
end

# SOLVING THE ODE FOR EACH TIME DIRECTION
function solve_ray(direction, tspan)
    parameters = (m, θ, E, Lz, direction)
    problem = ODEProblem(lightray, r0, tspan, parameters)
    return solve(problem, Tsit5(); reltol=1e-10, abstol=1e-10)
end

default(
    fontfamily="Computer Modern",
    linewidth=2,
    framestyle=:zerolines,
    grid=false,
    legend=false,
)

future_times = collect(range(t0, Tmax; length=400))
past_times = collect(range(t0, -Tmax; length=400))

# SOLVING FOR DIRECTION IN BOTH TIME AND SPACE
future_out = solve_ray(1.0, (t0, Tmax))
future_in = solve_ray(-1.0, (t0, Tmax))
past_out = solve_ray(1.0, (t0, -Tmax))
past_in = solve_ray(-1.0, (t0, -Tmax))

diagram = plot(
    xlabel="r",
    ylabel="t (Boyer-Lindquist coordinate time)",
    title="Minkowski Diagram of a Kerr Blackhole (a = $(a))",
)

# Shade the causal past and future in the equatorial, fixed-Lz radial slice.
plot!(diagram, Shape(
    vcat([future_out(t) for t in future_times],
         reverse([future_in(t) for t in future_times])),
    vcat(future_times, reverse(future_times)),
); color=:deepskyblue, fillalpha=0.25, linealpha=0, label=false)

plot!(diagram, Shape(
    vcat([past_out(t) for t in past_times],
         reverse([past_in(t) for t in past_times])),
    vcat(past_times, reverse(past_times)),
); color=:orange, fillalpha=0.25, linealpha=0, label=false)

for ray in (future_out, future_in, past_out, past_in)
    plot!(diagram, ray.u, ray.t; color=:black, label=false)
end

vline!(diagram, [horizon]; linestyle=:dash, color=:red, label="event horizon")
scatter!(diagram, [r0], [t0]; color=:black, markersize=4, label="observer")
display(diagram)
