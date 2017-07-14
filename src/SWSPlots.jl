"""
## SWSPlots

The SWSPlots module provides plotting functions to show maps with shear wave splitting
results

"""
module SWSPlots

using Geodesy, Plots
using SWSs

export
    plot_splits_bars,
    plot_splits_bars!

const SWSArr = Union{SWS, Array{SWS}}

"""
    plot_splits_bars!(p::Plots.Plot, splits::SWSArr) -> p::Plots.Plot
"""
function plot_splits_bars!(p::Plots.Plot, splits::SWSArr;
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

end # module
