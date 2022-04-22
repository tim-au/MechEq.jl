module BoltPattern

using Plots

using Unitful


export bolt_centroid, circle, rectangle, plot_bolt_pattern

function bolt_centroid(points; A=1, return_all_data = false, udf_pivot = false)
    x, y = points
    
    if length(A) == 1
        A = fill(A, length(x))
    end  

    if udf_pivot == false
        xc = sum( x .* A) ./ sum(A)
        yc = sum( y .* A) ./ sum(A)
    else
        xc, yc = udf_pivot
    end

    # distance of each bolt from the pattern bolt_centroid
    rcx = x .- xc
    rcy = y .- yc

    # centroidal moment of inertia
    Icx = sum(rcy.^2 .* A)
    Icy = sum(rcx.^2 .* A)
    Icp = Icx + Icy

    # shortest distance between bolt and centroidal
    rcxy = @. √(rcx^2 + rcy^2)

    if return_all_data == false
        return xc, yc
    elseif return_all_data == true
        return xc, yc, rcx, rcy, Icx, Icy, Icp, rcxy
    else
        println("error: return_all_data not correctly specified")
    end
end





function circle(;r, N=4, θ_start=0)
    β = deg2rad(θ_start)
    θ = LinRange(π/2, -3π/2 + 2π/N, N)
    points = @. r * cos(θ - β) , r * sin(θ - β) 

    return points
end



function rectangle(;x_dist, y_dist, Nx=2, Ny=2)
    x = LinRange(-x_dist / 2, x_dist / 2, Nx)
    y = LinRange(y_dist / 2, -y_dist / 2, Ny)

    y_outer = [y[1], y[end]]
    x_out = [ repeat([x[1]], inner = Ny) ; repeat(x[2:end-1], inner = 2) ; repeat([x[end]], inner = Ny) ]
    y_out = [ y ; repeat(y_outer, Nx - 2) ; y ]
    points = x_out, y_out

    return points
end



function plot_bolt_pattern(points; A = 1, udf_pivot = false)

    x,y = points

    if udf_pivot != false
        xc, yc = udf_pivot
    else
        xc, yc = bolt_centroid(points, A=A)
    end

    unit_type = string(unit(x[1]))

    # get range of x,y points
    x_plot = ustrip(x)
    y_plot = ustrip(y)
    xrange = (maximum(x_plot) - minimum(x_plot))
    yrange = (maximum(y_plot) - minimum(y_plot))

    # determine xlims and ylims for the plot
    scale = 0.1
    xs = minimum(x_plot) - scale * xrange
    xf = maximum(x_plot) + scale * xrange
    ys = minimum(y_plot) - scale * yrange
    yf = maximum(y_plot) + scale * yrange
  
    scatter(x_plot, y_plot,
            xlims = [xs, xf],
            ylims = [ys, yf], 
            marker = :hexagon, 
            markersize = 10,
            markercolor = :grey,
            markerstrokecolor = :black,
            aspect_ratio = 1,
            alpha = 1,
            legend = false,
            xlabel = "x coordinate [$unit_type]",
            ylabel = "y coordinate [$unit_type]")
    
    # plot centroid
    xc_hover = map(x -> "$(round.(ustrip(x), digits=1)) $(unit(x))", xc)
    yc_hover = map(x -> "$(round.(ustrip(x), digits=1)) $(unit(x))", yc)

    scatter!([ustrip(xc)], [ustrip(yc)], markercolor = :grey, markersize = 8, markerstrokewidth = 1, 
                markeralpha = 0.7, hover = "Centroid<br>xc: $(xc_hover)<br>yc: $(yc_hover)")
    vline!([ustrip(xc)],
            color = :grey)
    hline!([ustrip(yc)],
            color = :grey)
end



end