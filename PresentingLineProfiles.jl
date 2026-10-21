using Gradus, Plots, LaTeXStrings, Serialization

function Gradus.inner_radius(m::Gradus.KerrMetric{T}) where {T}
    disc = m.M^2 - m.a^2
    disc < 0 ? 4eps(T) : m.M + √disc
end

# ---------- SETTINGS ----------
θ_obs  = 60.0                  # the paper's figure uses 75
spins  = [0.998, 1.05]         # the paper's dashed curve is a = 1.002
γs     = 3:7
E_line = 6.4
RES    = :high                  # :low, :medium or :high

x = SVector(0.0, 10_000.0, deg2rad(θ_obs), 0.0)

bins, plane =
    RES == :low    ? (collect(range(0.1, 1.7, length = 100)),
                      PolarPlane(GeometricGrid(); Nr = 300,  Nθ = 700,  r_max = 45.0)) :
    RES == :medium ? (collect(range(0.1, 1.7, length = 200)),
                      PolarPlane(GeometricGrid(); Nr = 800,  Nθ = 1600, r_max = 45.0)) :
                     (collect(range(0.1, 1.7, length = 300)),
                      PolarPlane(GeometricGrid(); Nr = 1600, Nθ = 3200, r_max = 45.0))

# ---------- LINE PROFILE ----------
function CalculateLineProfile(metric, pos, disc; y = 4.0, kwargs...)
    return lineprofile(
        bins,
        r -> r^(-y),
        metric, pos, disc;
        verbose = true,
        method = BinningMethod(),
        plane = plane,
        maxrₑ = 30.0,
        kwargs...,
    )
end

function UnitChange(g, F)
    observed_energy = g .* E_line
    normalised_flux = g .* F ./ maximum(g .* F)
    return observed_energy, normalised_flux
end

# ---------- SWEEP ----------
results = Dict{Tuple{Float64,Int},Any}()

for a in spins
    m = KerrMetric(M = 1.0, a = a)
    d = ThinDisc(0.0, Inf)
    for γ in γs
        @info "Computing" a γ
        g, F = CalculateLineProfile(m, x, d; y = Float64(γ))
        results[(a, γ)] = UnitChange(g, F)
    end
end

serialize("lineprofiles_$(RES).jls", results)   # reload later with deserialize() to replot without recomputing

# ---------- PLOTTING ----------
default(
    fontfamily = "Computer Modern",
    framestyle = :box,
    grid = true,
    tickdirection = :in,
    minorticks = 5,
    guidefontsize = 11,
    tickfontsize = 9,
    titlefontsize = 11,
)

colors = Dict(3 => :blue, 4 => :deepskyblue, 5 => :lightgreen, 6 => :orange, 7 => :darkred)
styles = Dict(spins[1] => :solid, spins[2] => :dash)

p = plot(;
    xlabel = L"Photon Energy $E$, keV",
    ylabel = "Arbitrarily normalised flux",
    legend = :topleft,
    xlims = (0.5, 11.0),
    ylims = (0.0, 1.05),
)

vline!(p, [E_line]; color = :black, linestyle = :dashdot, label = "")

for a in spins, γ in γs
    e, f = results[(a, γ)]
    plot!(p, e, f;
        color = colors[γ],
        linestyle = styles[a],
        linewidth = 1.5,
        label = a == spins[1] ? L"\gamma = %$γ" : "",   # one legend entry per γ
    )
end

# legend entries for line style
for a in spins
    plot!(p, [NaN], [NaN]; color = :black, linestyle = styles[a], label = L"a = %$a")
end

display(p)
savefig(p, "lineprofiles_emissivity.pdf")