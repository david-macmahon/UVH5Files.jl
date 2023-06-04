module UVH5Files

using HDF5
using DataFrames

export UVH5File, corrcoeff, loadbaseline

struct UVH5File
    h5::HDF5.File
    all::DataFrame
    autos::GroupedDataFrame
    crosses::GroupedDataFrame

    function UVH5File(h5::HDF5.File)
        h5h = h5["Header"]
        # Get all baseline-time tuples
        a1s = h5h["ant_1_array"][]
        a2s = h5h["ant_2_array"][]
        jds = h5h["time_array"][]

        all = [(; a1, a2, jd, idx) for (idx, (a1, a2, jd)) in enumerate(zip(a1s, a2s, jds))] |> DataFrame;
        autos   = groupby(subset(all, [:a1,:a2]=>ByRow((a1,a2)->a1==a2), view=true), :jd)
        crosses = groupby(subset(all, [:a1,:a2]=>ByRow((a1,a2)->a1!=a2), view=true), :jd)

        new(h5, all, autos, crosses)
    end
end

function UVH5File(h5obj::Union{HDF5.Attribute, HDF5.Dataset, HDF5.Datatype, HDF5.Group})
    UVH5File(HDF5.file(h5obj))
end

function UVH5File(h5name::AbstractString)
    UVH5File(h5open(h5name))
end

import Base.close
function close(uv::UVH5File)
    close(uv.h5)
end

function ants(uv::UVH5File)
    sort(uv.all.a1 ∪ uv.all.a2)
end

import Base.==
function ==(a::UVH5File, b::UVH5File)
    a.all == b.all
end

function Base.show(io::IO, uv::UVH5File)
    nants = length(ants(uv))
    nautos = nrow(parent(uv.autos))
    ncross = nrow(parent(uv.crosses))
    ntautos = length(uv.autos)
    ntcross = length(uv.crosses)
    print(io, "UVH5File($nants ants: $nautos autos/$ntautos times, $ncross crosses/$ntcross times)")
end

# uv[Integer, Integer, Integer]

function Base.haskey(uv::UVH5File, a1::Integer, a2::Integer, tidx::Integer)
    gdf  = (a1 == a2) ? uv.autos : uv.crosses
    tidx in axes(gdf, 1) || return false
    any(gdf[tidx].a1 .== a1 .&& gdf[tidx].a2 .== a2)
end

function Base.getindex(uv::UVH5File, a1::Integer, a2::Integer, tidx::Integer)
    sdf = (a1 == a2) ? uv.autos[tidx] : uv.crosses[tidx]
    sdf.idx[sdf.a1.==a1 .&& sdf.a2.==a2]
end

# uv[Integer, Integer, AbstractFloat]

function Base.haskey(uv::UVH5File, a1::Integer, a2::Integer, jd::AbstractFloat)
    gdf  = (a1 == a2) ? uv.autos : uv.crosses
    haskey(gdf, (jd,)) || return false
    any(gdf[(jd,)].a1 .== a1 .&& gdf[(jd,)].a2 .== a2)
end

function Base.getindex(uv::UVH5File, a1::Integer, a2::Integer, jd::AbstractFloat)
    sdf = (a1 == a2) ? uv.autos[(jd,)] : uv.crosses[(jd,)]
    sdf.idx[sdf.a1.==a1 .&& sdf.a2.==a2]
end

# uv[Integer, Integer]
# uv[Integer, Integer, :]

function Base.haskey(uv::UVH5File, a1::Integer, a2::Integer, _::Colon=:)
    df = uv.all
    any(df.a1 .== a1 .&& df.a2 .== a2)
end

function Base.getindex(uv::UVH5File, a1::Integer, a2::Integer, _::Colon=:)
    df = uv.all
    df.idx[df.a1.==a1 .&& df.a2.==a2]
end

# uv[:, :, Integer]

function Base.haskey(uv::UVH5File, _::Colon, _::Colon, tidx::Integer)
    tidx in axes(unique(uv.all.jds), 1)
end

function Base.getindex(uv::UVH5File, _::Colon, _::Colon, tidx::Integer)
    # It is possible for autos and crosses to have a different number of time
    # samples, so we have to bounds check tidx against both.
    gdfa = uv.autos
    idxa = tidx in axes(gdfa, 1) ? gdfa.idx : []
    gdfc = uv.crosses
    idxc = tidx in axes(gdfc, 1) ? gdfc.idx : []
    sort(vcat(idxa, idxc))
end

# uv[:, :, AbstractFloat]

function Base.haskey(uv::UVH5File, _::Colon, _::Colon, jd::AbstractFloat)
    jd in unique(uv.all.jds)
end

function Base.getindex(uv::UVH5File, _::Colon, _::Colon, jd::AbstractFloat)
    df = uv.all
    df.idx[df.jd .== jd]
end

# uv[Integer] - All baseline-time indices involving given ant as a1 or a2
# uv[Integer, :]

function Base.haskey(uv::UVH5File, a::Integer, _::Colon=:)
    df = uv.all
    any(df.a1 .== a .|| df.a2 .== a)
end

function Base.getindex(uv::UVH5File, a::Integer, _::Colon=:)
    df = uv.all
    df.idx[df.a1.==a .|| df.a2.==a]
end

# uv[:, Integer] - All auto-correlation baseline-time indices for given time index

function Base.haskey(uv::UVH5File, _::Colon, tidx::Integer)
    haskey(uv.autos, tidx)
end

function Base.getindex(uv::UVH5File, _::Colon, tidx::Integer)
    uv.autos[tidx].idx
end

# Load baseline method

function loadbaselines(vis::HDF5.Dataset, idxs::AbstractArray{<:Integer})
    mapreduce(i->vis[:,:,i], (a,b)->cat(a,b,dims=3), idxs)
end

function loadbaseline(vis::HDF5.Dataset, uv, a1, a2)
    mapreduce(i->vis[:,:,i], (a,b)->cat(a,b,dims=3), uv[a1, a2, :])
end

function loadbaseline(uv, a1, a2)
    vis = uv.h5["Data/visdata"]
    try
        return loadbaseline(vis, uv, a1, a2)
    finally
        close(vis)
    end
end

# Corrcoeff methods

function corrcoeff(vis::AbstractArray{T,3}, uv::UVH5File, a1, a2) where T
    @assert a1 != a2 "a1 must not equal a2"
    v11 = @view vis[:, :, uv[a1,a1]]
    v22 = @view vis[:, :, uv[a2,a2]]
    v12 = haskey(uv, a1, a2) ? @view(vis[:, :, uv[a1,a2]]) :
          haskey(uv, a2, a1) ? @view(vis[:, :, uv[a2,a1]]) :
          error("no (a1,a2) or (a2,a1) baselines")
    v12 ./ sqrt.(abs.(v11) .* abs.(v22))
end

function corrcoeff(vis::AbstractArray{T,2}, uv::UVH5File, a1, a2) where T
    @assert a1 != a2 "a1 must not equal a2"
    v11 = @view vis[:, uv[a1,a1]]
    v22 = @view vis[:, uv[a2,a2]]
    v12 = haskey(uv, a1, a2) ? @view(vis[:, uv[a1,a2]]) :
          haskey(uv, a2, a1) ? @view(vis[:, uv[a2,a1]]) :
          error("no (a1,a2) or (a2,a1) baselines")
    v12 ./ sqrt.(abs.(v11) .* abs.(v22))
end

function corrcoeff(vis::HDF5.Dataset, uv::UVH5File, a1, a2)
    @assert a1 != a2 "a1 must not equal a2"
    v11 = loadbaseline(vis, uv, a1, a1)
    v22 = loadbaseline(vis, uv, a2, a2)
    v12 = haskey(uv, a1, a2) ? loadbaseline(vis, uv, a1, a2) :
          haskey(uv, a2, a1) ? loadbaseline(vis, uv, a2, a1) :
          error("no (a1,a2) or (a2,a1) baselines")
    @assert size(v11) == size(v12) "a1 auto size != cross size"
    @assert size(v22) == size(v12) "a2 auto size != cross size"
    v12 ./ sqrt.(abs.(v11) .* abs.(v22))
end

function corrcoeff(h5::HDF5.File, uv::UVH5File, a1, a2)
    corrcoeff(h5["Data/visdata"], uv, a1, a2)
end

"""
Load corrcoefs for all baselines in `vis` involving `refant`.
Returns a Dict{<:Integer,<:AbstractArray}.
"""
function corrcoeff(vis::Union{<:AbstractArray{T,2},<:AbstractArray{T,3},HDF5.Dataset}, uv::UVH5File, refant::Integer) where T
    allants = ants(uv)
    refant in allants || error("refant is not present")
    xants = filter(!=(refant), allants)
    Dict(xants .=> corrcoeff.(Ref(vis), Ref(uv), refant, xants))
end

function corrcoeff(h5::HDF5.File, uv::UVH5File, refant::Integer)
    corrcoeff(h5["Data/visdata"], uv, refant)
end

end # module UVH5Files
