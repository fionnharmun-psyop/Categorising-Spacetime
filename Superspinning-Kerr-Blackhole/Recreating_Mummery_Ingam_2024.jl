using Gradus, Plots

# metric and metric parameters
m = KerrNewmanMetric(M=1.0, a=0.8, Q=0.0)
# observer position
x = SVector(0.0, 1000.0, deg2rad(80.0), 0.0)
# accretion disc
d = ThinDisc(1.0, 30.0)

# define point function which filters geodesics that intersected the accretion disc
# and use those to calculate redshift
pf = ConstPointFunctions.redshift(m, x) ∘ ConstPointFunctions.filter_intersected()


# For wider view set alims to pm35 and blimes to pm15
# For Narrower view set alims to pm10 and blimes to pm5
α, β, img = rendergeodesics(
    m,
    x,
    d,
    # Longer integration time for more accurate results
    20000.0,
    αlims = (-35, 35),
    βlims = (-15, 15),
    image_width = 2560,
    image_height = 1440,
    verbose = true,
    pf = pf,
)

heatmap(α, β, img)