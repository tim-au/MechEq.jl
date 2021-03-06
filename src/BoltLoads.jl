#####################################
#####################################

module BoltLoads

export bolt_loads, plot_bolt_loads

####################################
###################################
include("BoltPattern.jl")

using .BoltPattern: bolt_centroid

using LinearAlgebra

using Plots

using DataFrames

using Unitful
using Unitful: mm, cm, m, inch, N, lbf, kN

########################################
########################################
"""
    bolt_loads(points; Fc=[0,0,0]N, Mc = [0,0,0]N*m, A=1mm^2, udf_pivot=false, length_format = "mm", load_format = "kN", return_plot_data=false)

Compute the bolt loads for the forces and moments acting about the centroid of the bolt pattern.

`points` is a tuple of x & y coordinates ([x], [y]).  The `points` may be generated by functions [`circle`](@ref) and [`rectangle`](@ref).
`Fc` is the force acting at the centroid of the bolt pattern entered as a vector [Fxc, Fyc, Fzc]. 
`Mc` is the moment acting at the centroid of the bolt pattern entered as a vector [Mxc, Myc, Mzc].
`udf_pivot` is false when it is to be calculated, or enetered as [x,y] when a specific pivot point is desired to be manually entered.
`A` is the bolt stress area, either entered as a single value for all bolts or different values for each bolt (i.e. vector).
`length_format` may be 'mm', `cm`, 'm` or `inch`, and sets the display output returned.
`load_format` may by 'N', 'kN', or `lbf`, and sets the display output returned.
Units conforming to the `Unitful` package should be applied.

# Examples
```julia-repl
julia> using Unitful: mm, m, N
julia> p = x, y = circle(r=100mm, N=7);
julia> Fc = [7000, 7000, 5000]N;
julia> Mc = [30_000, 20_000, 15_000]N*m;
julia> bolt_loads(p, Fc=Fc, Mc=Mc)
7×4 DataFrame
 Row │ x               y            Paxial       Pshear     
     │ Quantity…       Quantity…    Quantity…    Quantity…  
─────┼──────────────────────────────────────────────────────
   1 │ 6.12323e-15 mm     100.0 mm   86.4286 kN   20.453 kN
   2 │     78.1831 mm    62.349 mm   9.48018 kN  21.6326 kN
   3 │     97.4928 mm  -22.2521 mm  -74.0691 kN  22.6385 kN
   4 │     43.3884 mm  -90.0969 mm  -101.305 kN  22.7682 kN
   5 │    -43.3884 mm  -90.0969 mm  -51.7183 kN  21.9363 kN
   6 │    -97.4928 mm  -22.2521 mm   37.3512 kN  20.7108 kN
   7 │    -78.1831 mm    62.349 mm   98.8324 kN  20.0239 kN
```
"""
function bolt_loads(points; Fc=[0,0,0]N, Mc = [0,0,0]N*m, A=1mm^2, udf_pivot=false,
                    length_format = "mm", load_format = "kN", return_plot_data=false)

    x,y = points

    if length(A) == 1
        A = fill(A, length(x))
    end

    Fcx, Fcy, Fcz = Fc
    Mcx, Mcy, Mcz = Mc

    Mcz = [0N*m,0N*m, Mcz]

    xc, yc, rcx, rcy, Icx, Icy, Icp, rcxy = bolt_centroid(points, A=A, return_all_data = true, udf_pivot = udf_pivot)
    centroid = xc, yc

    rc = [ [rcx[i], rcy[i], 0mm] for i in eachindex(rcx) ]

    # axial load in bolt
    sumA = sum(A)
    Pz_Fz = @. Fcz * A / sumA
    Pz_Mx = @. Mcx * rcy * A / Icx
    Pz_My = @. -Mcy * rcx  * A / Icy

    Paxial = Pz_Fz + Pz_Mx + Pz_My

    # shear load
    Px_Fx = @. Fcx * A / sumA
    Py_Fy = @. Fcy * A / sumA

    P_FxFy = [ [Px_Fx[i], Py_Fy[i], 0kN] for i in eachindex(Px_Fx) ]

    P_Mz = map( (rc, A) -> cross(Mcz, rc) * A / Icp, rc, A )

    Pshear = P_FxFy + P_Mz

    PshearMag = norm.(Pshear)

    if length_format == "mm"
        x = x .|> mm
        y = y .|> mm
    elseif length_format == "cm"
        x = x .|> cm
        y = y .|> cm
    elseif length_format == "m"
        x = x .|> m
        y = y .|> m
    elseif length_format == "inch"
        x = x .|> inch
        y = y .|> inch
    end

    if load_format == "N"
        Paxial = Paxial .|> N
        PshearMag = PshearMag .|> N
    elseif load_format == "kN"
        Paxial = Paxial .|> kN
        PshearMag = PshearMag .|> kN
    elseif load_format == "lbf"
        Paxial = Paxial .|> lbf
        PshearMag = PshearMag .|> lbf
    end



    # dataframe to return
    df = DataFrame(x = x, y = y, Paxial = Paxial, Pshear = PshearMag)
    
    if return_plot_data
        return x, y, Paxial, Pshear, PshearMag
    else
        return df
    end
end

############################################
############################################
"""
    plot_bolt_loads(points; Fc=[0,0,0]N, Mc = [0,0,0]N*m, A=1mm^2, udf_pivot=false, length_format = "mm", load_format = "kN")

Plot the bolt loads for the forces and moments acting about the centroid of the bolt pattern.

    `points` is a tuple of x & y coordinates ([x], [y]).  The `points` may be generated by functions [`circle`](@ref) and [`rectangle`](@ref).
    `Fc` is the force acting at the centroid of the bolt pattern entered as a vector [Fxc, Fyc, Fzc]. 
    `Mc` is the moment acting at the centroid of the bolt pattern entered as a vector [Mxc, Myc, Mzc].
    `udf_pivot` is false when it is to be calculated, or enetered as [x,y] when a specific pivot point is desired to be manually entered.
    `A` is the bolt stress area, either entered as a single value for all bolts or different values for each bolt (i.e. vector).
    `length_format` may be 'mm', `cm`, 'm` or `inch`, and sets the display output returned.
    `load_format` may by 'N', 'kN', or `lbf`, and sets the display output returned.
    Units conforming to the `Unitful` package should be applied.
"""
function plot_bolt_loads(points; Fc=[0,0,0]N, Mc = [0,0,0]N*m, A=1mm^2, udf_pivot=false, length_format = "mm", load_format = "kN")
    
    x, y, Paxial, Pshear, PshearMag = bolt_loads(points; Fc=Fc, Mc=Mc, A=A, udf_pivot=udf_pivot,
                                                    length_format = length_format, load_format = load_format,
                                                    return_plot_data = true)

    xc, yc = bolt_centroid(points; A=A, udf_pivot=udf_pivot)

    # cover Pshear to Vx and Vy vector
    Vx = map(a -> ustrip(a[1]), Pshear)
    Vy = map(a -> ustrip(a[2]), Pshear)
    Vxy = sqrt.(Vx.^2 + Vy.^2)

    # covert to desired length format
    if length_format == "mm"
        xc = xc .|> mm
        yc = yc .|> mm
    elseif length_format == "cm"
        xc = xc .|> cm
        yc = yc .|> cm
    elseif length_format == "m"
        xc = xc .|> m
        yc = yc .|> m
    elseif length_format == "inch"
        xc = xc .|> inch
        yc = yc .|> inch
    end


    # remove all units now that they are in the required format 
    x_plot = ustrip(x)
    y_plot = ustrip(y)
    xc_plot = ustrip(xc)
    yc_plot = ustrip(yc)
    Paxial_plot = ustrip(Paxial)
    PshearMag_plot = ustrip(PshearMag)

    # convert Paxial to color map
    check(x) = x>=0 ? x : 0

    color_map = map(check, Paxial_plot)
    
    # get range of x,y points
    xrange = (maximum(x_plot) - minimum(x_plot))
    yrange = (maximum(y_plot) - minimum(y_plot))
    xyrange = maximum([xrange, yrange])

    # determine xlims and ylims for the plot
    scale = 0.30
    xs = minimum(x_plot) - scale * xrange
    xf = maximum(x_plot) + scale * xrange
    ys = minimum(y_plot) - scale * yrange
    yf = maximum(y_plot) + scale * yrange
    
    # hover text for plot
    x_hover = map(x -> "$(round.(x, digits=1)) $(length_format)", x_plot)
    y_hover = map(x -> "$(round.(x, digits=1)) $(length_format)", y_plot)
    Paxial_hover = map(x -> "$(round.(x, digits=1)) $(load_format)", Paxial_plot)
    PshearMag_hover = map(x -> "$(round.(x, digits=1)) $(load_format)", PshearMag_plot)
    xc_hover = map(x -> "$(round.(x, digits=1)) $(length_format)", xc_plot)
    yc_hover = map(x -> "$(round.(x, digits=1)) $(length_format)", yc_plot)
    id = collect(1:length(x_hover))


    hovertext = map((a, b, c, d, e) -> "ID:     $(a)<br>x:      $(b)<br>y:      $(c)<br>Axial: $(d)<br>Shear:  $(e)",
                    id,x_hover, y_hover, Paxial_hover, PshearMag_hover)

    # scale arrow length
    max_arrow_length = 0.2 * xyrange
    arrow_scale = max_arrow_length / maximum(Vxy)

    # vector for arrow body
    arrow_body_x = Vx * arrow_scale
    arrow_body_y = Vy * arrow_scale

    # vector for arrow head 1 and 2
    β = 180-30
    head_scale = 0.1

    head1_x = head_scale * (cosd(β)*arrow_body_x - sind(β)*arrow_body_y)
    head1_y = head_scale * (sind(β)*arrow_body_x + cosd(β)*arrow_body_y)
    head2_x = head_scale * (cosd(-β)*arrow_body_x - sind(-β)*arrow_body_y)
    head2_y = head_scale * (sind(-β)*arrow_body_x + cosd(-β)*arrow_body_y)

    # coordinates of arrow vector
    arrow_tip_x = x_plot + arrow_body_x
    arrow_tip_y = y_plot + arrow_body_y
    h1x = arrow_tip_x + head1_x
    h1y = arrow_tip_y + head1_y
    h2x = arrow_tip_x + head2_x
    h2y = arrow_tip_y + head2_y

    
    p = plot();
    # Plot Paxial data on scatter plot
    scatter!(x_plot, y_plot,
            hover = hovertext,
            xlims = [xs, xf],
            ylims = [ys, yf],
            aspect_ratio = 1,
            markersize = 12,
            marker = :hexagon,
            marker_z = color_map, 
            markercolor = :lighttest,
            markerstrokecolor = :grey,
            markerstrokewidth = 1,
            markeralpha = 1,
            colorbar=true,
            colorbar_title = "Paxial",
            legend = false);
     
    # Plot centroid on same plot
    scatter!([xc_plot], [yc_plot], markercolor = :grey, markersize = 8, markerstrokewidth = 1, 
                markeralpha = 0.7, hover = "Centroid<br>xc: $(xc_hover)<br>yc: $(yc_hover)");
    vline!([ustrip(xc_plot)],
            color = :grey);
    hline!([ustrip(yc_plot)],
            color = :grey);

    # plot Vxy
    for (xs, xf, ys, yf, h1x, h1y, h2x, h2y) in zip(x_plot, arrow_tip_x, y_plot, arrow_tip_y,
                                                    h1x, h1y, h2x, h2y)
        # plot body of arrow
        plot!([xs, xf], [ys, yf], linecolor = :pink, linewidth = 1);
        # plot head1 of arrow
        plot!([xf, h1x], [yf, h1y], linecolor = :pink, linewidth = 1);
        #plot head2 of arrow
        plot!([xf, h2x], [yf, h2y], linecolor = :pink, linewidth = 1);
    end


    # display plot "p"
    p


end




# end of module
end