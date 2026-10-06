using Gradus, Plots

function Gradus.inner_radius(m::Gradus.KerrMetric{T}) where {T}
    disc = m.M^2 - m.a^2 # Operand of the sqrt
    disc < 0 ? 4eps(T) : m.M + √disc # "If sqrt operand is less than 0 then return as close to zero as possible but not zero
end

# DEFINE METRICS
Sch = KerrMetric(M = 1.0, a = 0.0)
Kerr = KerrMetric(M=1.0, a = 0.998)
SuperKerr = KerrMetric(M=1.0, a = 1.5)
x = SVector(0.0, 10_000.0, deg2rad(60), 0.0)


# LINE PROFILE FUNCTION (binning method)
bins = collect(range(0.1, 1.5, length=360))
plane = PolarPlane(GeometricGrid(); Nr=450, Nθ=1300, r_max=250.0)

# For a higher-resolution final plot, try 720 bins and Nr=900, Nθ=2600.
# Doubling both plane dimensions gives about 4x as many rays, so runtime may
# rise from about 20 seconds to roughly 80 seconds per metric; check convergence.
# bins = collect(range(0.1, 1.5, length=720))
# plane = PolarPlane(GeometricGrid(); Nr=900, Nθ=2600, r_max=250.0)
function CalculateLineProfile(metric, pos, disc; kwargs...)
    bins, flux = lineprofile(metric, 
    pos, 
    disc; 
    verbose = true, 
    method = BinningMethod(),
    plane = plane, 
    kwargs...)
end

# DEFINE THE DISKS
sch_d = ThinDisc(0.0, Inf)
kerr_d = ThinDisc(0.0, Inf)
Skerr_d = ThinDisc(0.0, Inf)

# CALCULATE THE LINE PROFILES
sch_bins, sch_line = CalculateLineProfile(Sch, x, sch_d; bins = bins)
kerr_bins, kerr_line = CalculateLineProfile(Kerr, x, kerr_d; bins = bins)
skerr_bins, skerr_line = CalculateLineProfile(SuperKerr, x, Skerr_d; bins=bins)

# PLOTTING
default(
    fontfamily="Computer Modern",
    linewidth=2.0,
    framestyle=:zerolines,
    label=nothing,
    grid=false,
    legend=false,
)

p = plot(; 
title = "Line Profile of Different Metrics",
xlabel = "Redshift (g, Eobs/Eemit)",
ylabel = "Flux",
legend=:topleft,
)

plot!(sch_bins, sch_line, label = "Schwarzschild")
plot!(kerr_bins, kerr_line, label = "Kerr")
plot!(skerr_bins, skerr_line, label = "Super Spinning")
display(p)