using Gradus, Plots, LaTeXStrings

function Gradus.inner_radius(m::Gradus.KerrMetric{T}) where {T}
    disc = m.M^2 - m.a^2 # Operand of the sqrt
    disc < 0 ? 4eps(T) : m.M + √disc # "If sqrt operand is less than 0 then return as close to zero as possible but not zero
end

Sch = KerrMetric(M = 1.0, a = 0.0)
Kerr = KerrMetric(M=1.0, a = 0.8)
SuperKerr = KerrMetric(M=1.0, a = 1.5)
x = SVector(0.0, 1000.0, deg2rad(80.0), 0.0)

# define custom bins for g
bins = collect(range(0.1, 1.4, 200))
# define the plane to perform the binning over
plane = PolarPlane(GeometricGrid(); Nr = 1000, Nθ = 1000, r_max = 50.0)

plot(
    PolarPlane(GeometricGrid(); Nr = 10, Nθ = 20, r_max = 50.0)
)

function calculate_line_profile(m, x, d, bins, plane)
    _, f = lineprofile(
        m,
        x,
        d,
        method = BinningMethod(),
        # no false images
        callback = domain_upper_hemisphere(),
        verbose = true,
        bins = bins,
        plane = plane,
    )
    return f
end

sch_d = ThinDisc(Gradus.isco(Sch), 200.0)
kerr_d = ThinDisc(Gradus.isco(Kerr), 200.0)
#skerr_d = ThinDisc(Gradus.isco(SuperKerr), 200.0)

sch_line = calculate_line_profile(m, x, sch_d, bins, plane)
kerr_line = calculate_line_profile(m, x, kerr_d, bins, plane)
#skerr_line = calculate_line_profile(m, x, skerr_d, bins, plane)

