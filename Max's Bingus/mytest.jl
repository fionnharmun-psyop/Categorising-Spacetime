using Gradus
#defining a spactime
struct Schwarzschild{T} <: AbstractStaticAxisSymmetric{T} #Defining a struct that parameterises spacetime as a subtype of
    M::T                                                  #AbstractStaticAxisSymmetric, since we know we want our spacetime to be static (no time dependence) and axisymmetric (no dependence on the azimuthal angle phi).
end                                             #T is the number type of this metric.

function Gradus.metric_components(m::Schwarzschild, x)     #implementation of metric. x only contains r and theta because it is static and axisymmetric, so we don't need to include t or phi in the metric components.
    r, theta = x
    M = m.M

    dt2 = -(1 - (2M/r))
    dr2 = -inv(dt2)
    dtheta2 = r^2
    dphi2 = r^2 * sin(theta)^2
    dtdphi = zero(r)

    SVector(dt2, dr2, dtheta2, dphi2, dtdphi)
end

Gradus.inner_radius(m::Schwarzschild) = 2 * m.M   #specifying inner radius for integrattion to avoid numerical errors. Here it is set to the Scwarzschild radius

#photon trajectories
#The constrain function is automatically applied by the tracegeodesics function
# we place ourselves 1000r_g away from the singularity, which is at the origin.
# setup our spactime with choice of mass of M = 1.0 in standard GR units
# initial velocity arbitrarily set to v^r = -1 towards the singularity and v^phi very small so it grazes past the singularity.

m = Schwarzschild(1.0)
x = SVector(0.0, 1000.0, pi/2, 0.0) #initial position of the photon in (t,r,theta,phi) coordinates
v = SVector(1.0, -1.0, 0.0, -8e-6)

lambda_max = 2000.0 #maximum affine parameter to integrate to
sol = tracegeodesics(m, x, v, lambda_max) #trace the geodesics of the photon in the Scwarzschild spacetime

#plot the thing
using Plots

plot_paths(sol)
plot_horizon!(m)


