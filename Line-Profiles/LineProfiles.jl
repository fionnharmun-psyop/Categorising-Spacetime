using Gradus, Plots, LaTeXStrings

function Gradus.inner_radius(m::Gradus.KerrMetric{T}) where {T}
    disc = m.M^2 - m.a^2 # Operand of the sqrt
    disc < 0 ? 4eps(T) : m.M + √disc # "If sqrt operand is less than 0 then return as close to zero as possible but not zero
end

# DEFINE METRICS
Sch = KerrMetric(M = 1.0, a = 0.0)
Kerr = KerrMetric(M=1.0, a = 0.998)
SuperKerr = KerrMetric(M=1.0, a = 1.05)
x = SVector(0.0, 10_000.0, deg2rad(60.0), 0.0)


# LINE PROFILE FUNCTION (binning method)
# Medium Res 
#bins  = collect(range(0.1, 1.7, length = 200))       
#plane = PolarPlane(GeometricGrid(); Nr = 800, Nθ = 1600, r_max = 45.0)
# High Res (3 mins per plot)
bins  = collect(range(0.1, 1.7, length = 300))                     
plane = PolarPlane(GeometricGrid(); Nr = 1600, Nθ = 3200, r_max = 45.0)
# Low Res (7 secs per plot)
#bins  = collect(range(0.1, 1.7, length = 100))                     
# plane = PolarPlane(GeometricGrid(); Nr = 300, Nθ = 700, r_max 
function CalculateLineProfile(metric, pos, disc; bins = bins, y = 4.0 , kwargs...)
    return lineprofile(
        bins,                    
        r -> r^(-y),             # emissivity function (can be replaced with a disc emissivity model)
        metric,
        pos,
        disc;
        verbose = true,
        method = BinningMethod(),
        plane = plane,
        maxrₑ = 30.0,
        kwargs...,
    )
end

# DEFINE THE DISKS
sch_d = ThinDisc(0.0, Inf)
kerr_d = ThinDisc(0.0, Inf)
Skerr_d = ThinDisc(0.0, Inf)

# CALCULATE THE LINE PROFILES
#sch_bins, sch_line = CalculateLineProfile(Sch, x, sch_d; bins = bins)
kerr_bins, kerr_line = CalculateLineProfile(Kerr, x, kerr_d; y = 4.0)
skerr_bins, skerr_line = CalculateLineProfile(SuperKerr, x, Skerr_d; y = 4.0)

# ADJUST UNITS OF LINE PROFILE TO MATCH M&I
# Changes gravitational redshift and flux to observed energy and conventional normalised flux
# The flux should match that in Mummery & Ingram (2024) Figure 9, i.e., E^2 dN_gamma/dE
function UnitChange(g, F) 
    observed_energy = g .* 6.4
    normalised_flux = g.^2 .* F ./ maximum(g.^2 .* F)
    return observed_energy, normalised_flux
end

#sch_e , sch_f = UnitChange(sch_bins, sch_line)
kerr_e , kerr_f = UnitChange(kerr_bins, kerr_line)
skerr_e, skerr_f = UnitChange(skerr_bins, skerr_line)


# PLOTTING
default(
    fontfamily = "Computer Modern",   
    framestyle = :box,                
    grid = true,
    tickdirection = :in,              
    minorticks = 5,                   
    guidefontsize = 11,               
    tickfontsize = 9,
    titlefontsize = 11,
    colorbar_titlefontsize = 11,
)

p = plot(; 
title = "Line Profile of Different Metrics",
xlabel = "Photon Energy E (keV)",
ylabel = L"Normalised Flux ($E^2\,\mathrm{d}N_{\gamma}/\mathrm{d}E$)",
legend=:topright,
ylims = (0.0, 1.2),
xlims = (2.0, 11.0)
)

E_line = 6.4
vline!([E_line], color = :black, linestyle = :dashdot, label = "6.4 keV")
#plot!(sch_e, sch_f, label = "Schwarzschild")
plot!(kerr_e, kerr_f, label = L"a = 0.998")
plot!(skerr_e, skerr_f, label = L"a = 1.05")
display(p)