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
x = SVector(0.0, 1000.0, deg2rad(70), 0.0) #rotating the camera, or the x plane
lambda_max = 2000.0 #maximum affine parameter to integrate to

#IMPACT PARAMETERS

# grid of impact parameters in horizontal direction
# keeping β fixed at 0
α = range(-10.0, 10.0, 30)
vs = map_impact_parameters(m, x, α, 0.0)

# need a position for each velocity vector
xs = fill(x, size(vs))

sols = tracegeodesics(m, xs, vs, lambda_max)

using Plots
p = plot_paths(sols, legend=false)
plot_horizon!(m, color = :black)

#time based imaging
α = range(-10.0, 10.0,100) #100x100 matrix of velocity vectors
β = range(-10.0, 10.0,100)
#flatten object using vec
vs = vec([map_impact_parameters(m, x, a, b) for a in α, b in β])
xs = fill(x, size(vs))

#trace geodesics
sols = tracegeodesics(m, xs, vs, lambda_max, save_on = false)

points = unpack_solution.(sols.u)
points =  reshape(points ,(100,100))

times = map(points) do gp
    if gp.status == StatusCodes.WithinInnerBoundary
        gp.x[1]
    else
        NaN
    end
end

heatmap(α, β, times, xlabel = "α", ylabel = "β", colorbar_title = "Time to reach horizon", title = "Time to reach horizon for different impact parameters", aspect_ratio = 1)

time_coord = PointFunction((m, gp, lambda) -> gp.x[1])
filter_event_horizon = FilterStatusCode(StatusCodes.WithinInnerBoundary)
pf_geometry = time_coord ∘ ConstPointFunctions.filter_intersected()
times = pf_geometry.(m, points, lambda_max)

#adding an accretion disk to the image
#define cross section of accretion disk
function cross_section(x)
    center = 8
    radius = 3

    if (x < center - radius) || (x > center + radius)
        zero(x)
    else
        r = x - center
        sqrt(radius^2 - r^2) + (0.5sin(3x))
    end
end

#lessavalook
sample = collect(range(0.0,20.0,300))
y = cross_section.(sample)

#plot(sample, y, xlabel = "r", ylabel = "height", aspect_ratio = 1)

d = ThickDisc(r -> cross_section(r))


# this function returns the impact parameter axes

x = SVector(0.0, 1000.0, deg2rad(40), 0.0) #rotating the camera, or the x plane
α, β, image = rendergeodesics(
    m,
    x,
    d,
    # no longer need to specify the velocities
    # these are automatically calculated
    lambda_max,
    pf = pf_geometry,
    # image parameters
    image_width = 1200,
    image_height = 800,
    # the "zoom" -- use the impact parameter axes
    αlims = (-20, 20),
    βlims = (-15, 15),
    verbose = true
)

heatmap(α, β, image, aspect_ratio = 1)
