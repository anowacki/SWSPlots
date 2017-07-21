"""
## SWSPlots

The SWSPlots module provides plotting functions to show maps with shear wave splitting
results.

"""
module SWSPlots

using Geodesy, Plots

using SWSs
using CircPlot

export
    plot_splits_bars,
    plot_splits_bars!,
    plot_splits_hist,
    plot_splits_hist!

const SWSArr = Union{SWS, Array{SWS}}
const PlotSubplot = Union{Plots.Plot,Plots.Subplot}

"""
    plot_splits_bars!(p::Plots.Plot, splits::SWSArr) -> p::Plots.Plot
"""
function plot_splits_bars!(p::PlotSubplot, splits::SWSArr;
        km=false, scale=km?10:1e4, location=:station, kwargs...)
    s = splits
    x, y = if location == :station
        s[:sx], s[:sy]
    elseif location == :midpoint
        (s[:sx] + s[:ex])./2, (s[:sy] + s[:ey])./2
    elseif location == :event
        s[:ex], s[:ey]
    else
        throw(ArgumentError("`location` ($location) must be one of:" *
            " $([:station, :location, :midpoint])"))
    end
    if km
        x, y = x/1e3, y/1e3
    end
    ϕ = deg2rad.(s[:phi])
    dt = s[:dt]
    u, v = scale.*dt.*sin.(ϕ), scale.*dt.*cos.(ϕ)
    for (U,V) in zip((u,-u), (v,-v))
        quiver!(p, x, y, quiver=(U,V), arrow=(0,0), line=(stroke(4,:blue)); kwargs...)
    end
    p
end
plot_splits_bars(splits::SWSArr; kwargs...) = plot_splits_bars!(plot(), splits; kwargs...)

function plot_splits_hist!(p::PlotSubplot, splits::SWSArr;
        km=false, binwidth=20, maxsize=km?5:5e3, minobs=3, kwargs...)
    s = splits
    stas, sx, sy = sta_coords(s)
    km && (sx /= 1e3; sy /= 1e3)
    hists = Array{Array{Plots.Shape}}(0)
    colours = [:red, :green, :blue, :orange, :white]
    for (i,sta) in enumerate(sort(stas))
        plot!(p, chistogram(s[s[:stnm] .== sta][:phi], binwidth, true, sx[i], sy[i],
            axial=true, azimuth=true, maxr=maxsize, minobs=minobs);
            c=colours[mod1(i, length(colours))], label=sta, kwargs...)
    end
    p
end
plot_splits_hist(splits::SWSArr; kwargs...) = plot_splits_hist!(plot(), splits; kwargs...)

"""
    sta_coords(s) -> stas, x, y

Return arrays of station names `stas`, x-coordinates, `x` and y-coordinates `y`.
"""
function sta_coords(s::SWSArr)
    stas = unique(s[:stnm])
    x, y = Float64[], Float64[]
    for sta in stas
        push!(x, s[s[:stnm] .== sta][1].sx)
        push!(y, s[s[:stnm] .== sta][1].sy)
    end
    stas, x, y
end

"""
    bar(ϕ, length, x=0, y=0, width=length/20) -> ::Shape

Return a `Plots.Shape` representing a bar oriented parallel to azimuth `ϕ` °
and `length` plot units long.  `width` specifies the bar width in plot units.
The bar is centred at `(x, y)`.
"""
function bar(ϕ, length, x=0., y=0., width=length/20.)
    l2 = length/2
    w2 = width/2
    ϕ = deg2rad(ϕ)
    sp, cp = sin(ϕ), cos(ϕ)
    p1x, p1y = l2*sp, l2*cp
    shape = Array{Tuple{Float64,Float64}}(4)
    i = 0
    for (lpm ,wpm) in zip((1, -1, -1, 1), (1, 1, -1, -1))
        i += 1
        shape[i] = (x + lpm*p1x + wpm*w2*sin(ϕ-π/2), y + lpm*p1y + wpm*w2*cos(ϕ-π/2))
    end
    Plots.Shape(shape)
end

end # module
