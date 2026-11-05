using Gradus, Plots

# NOT SURE ABOUT THIS COMMON FILE IDEA.

function Gradus.inner_radius(m::Gradus.KerrMetric{T}) where {T}
    disc = m.M^2 - m.a^2
    disc < 0 ? 4eps(T) : m.M + √disc
end

function UnitChange(g, F) 
    observed_energy = g .* 6.4
    normalised_flux = g .^2 .* F ./ maximum(g .^2 .* F)
    return observed_energy, normalised_flux
end
