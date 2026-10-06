using Gradus, Plots 

m = KerrMetric(M=1.0, a=0.9)
x = SVector(0.0, 1000.0, deg2rad(80), 0.0) # Camera Position in Boyer-Lindquist coords
max_int_time = 2000.0

# Accretion Disk
function cross_section(x)
    # Defining two concentric rings
    centre1 = 8
    centre2 = 19
    radius1 = 1
    radius2 = 0.5
    # Condition for first ring
    ring1 = if (x < centre1 - radius1) || (radius1 + centre1 < x)
        zero(x)
    else
        r = x - centre1
        sqrt(radius1^2 - r^2) + (0.5 * sin(3x))
    end
    # Condition for second ring
    ring2 = if (x < centre2 - radius2) || (radius2 + centre2 < x)
       zero(x)
    else
        r = x - centre2
        sqrt(radius2^2 - r^2) + (0.5 * sin(3x))
    end

    ring1 + ring2
end

#Plotting Cross Section
sample = collect(range(0.0, 20.0, 300))
y = cross_section.(sample)
plot(sample, y, xlabel = "r", ylabel = "height", aspect_ratio = 1)

#Wrapping as ThickDisk type
d = ThickDisc(r -> cross_section(r))

# Defining the point function for redshift and filtering intersected rays
pf = ConstPointFunctions.redshift(m, x) ∘ ConstPointFunctions.filter_intersected()

# Rendering
alpha, beta, image = rendergeodesics(
    m,
    x,
    d,
    max_int_time,
    αlims = (-20, 20),
    βlims = (-15, 15),
    image_width = 1600,
    image_height = 800,
    verbose = true,
    pf = pf,
)
heatmap(alpha, beta, image, aspect_ratio = 1)