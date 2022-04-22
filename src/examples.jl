using MechEq

# circle
p3 = x3, y3 = circle(r=100mm, N=6)
pc3 = bolt_centroid(p3)
plot_bolt_pattern(p3)
plot_bolt_pattern(p3, udf_pivot = [-30mm, 75mm])


# rectangle
p4 = rectangle(x_dist = 250mm, y_dist = 125mm, Nx = 3, Ny=4)
pc4 = bolt_centroid(p4)
plot_bolt_pattern(p4)
plot_bolt_pattern(p4, udf_pivot = (-48, 35))

# points
x5 = [-35, -30, -25, 27, 29, 45]mm
y5 = [-20, 12, 30, 27, -20, -50]mm
p5 = x5, y5
pc5 = bolt_centroid(p5)
plot_bolt_pattern(p5)
plot_bolt_pattern(p5, udf_pivot = (-11, -20))

# custom area for points p5
A5 = [1,1,1,1,1,20]mm^2
pc5_custom = bolt_centroid(p5, A = A5)
plot_bolt_pattern(p5, A = A5)