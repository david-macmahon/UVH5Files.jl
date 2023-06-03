module UVH5Files

using HDF5
using DataFrames

export BLTIndex, corrcoeff, loadbaseline

struct BLTIndex
    all::DataFrame
    autos::GroupedDataFrame
    crosses::GroupedDataFrame

    function BLTIndex(h5h::HDF5.Group)
        # Get all baseline-time tuples
        a1s = h5h["ant_1_array"][]
        a2s = h5h["ant_2_array"][]
        jds = h5h["time_array"][]

        all = [(; a1, a2, jd, idx) for (idx, (a1, a2, jd)) in enumerate(zip(a1s, a2s, jds))] |> DataFrame;
        autos   = groupby(subset(all, [:a1,:a2]=>ByRow((a1,a2)->a1==a2), view=true), :jd)
        crosses = groupby(subset(all, [:a1,:a2]=>ByRow((a1,a2)->a1!=a2), view=true), :jd)

        new(all, autos, crosses)
    end
end

function BLTIndex(h5::HDF5.File)
    BLTIndex(h5["Header"])
end

function ants(blt::BLTIndex)
    sort(blt.all.a1 ∪ blt.all.a2)
end

import Base.==
function ==(a::BLTIndex, b::BLTIndex)
    a.all == b.all
end

function Base.show(io::IO, blt::BLTIndex)
    nants = length(ants(blt))
    nautos = nrow(parent(blt.autos))
    ncross = nrow(parent(blt.crosses))
    ntautos = length(blt.autos)
    ntcross = length(blt.crosses)
    print(io, "BLTIndex($nants ants: $nautos autos/$ntautos times, $ncross crosses/$ntcross times)")
end

# blt[Integer, Integer, Integer]

function Base.haskey(blt::BLTIndex, a1::Integer, a2::Integer, tidx::Integer)
    gdf  = (a1 == a2) ? blt.autos : blt.crosses
    tidx in axes(gdf, 1) || return false
    any(gdf[tidx].a1 .== a1 .&& gdf[tidx].a2 .== a2)
end

function Base.getindex(blt::BLTIndex, a1::Integer, a2::Integer, tidx::Integer)
    sdf = (a1 == a2) ? blt.autos[tidx] : blt.crosses[tidx]
    sdf.idx[sdf.a1.==a1 .&& sdf.a2.==a2]
end

# blt[Integer, Integer, AbstractFloat]

function Base.haskey(blt::BLTIndex, a1::Integer, a2::Integer, jd::AbstractFloat)
    gdf  = (a1 == a2) ? blt.autos : blt.crosses
    haskey(gdf, (jd,)) || return false
    any(gdf[(jd,)].a1 .== a1 .&& gdf[(jd,)].a2 .== a2)
end

function Base.getindex(blt::BLTIndex, a1::Integer, a2::Integer, jd::AbstractFloat)
    sdf = (a1 == a2) ? blt.autos[(jd,)] : blt.crosses[(jd,)]
    sdf.idx[sdf.a1.==a1 .&& sdf.a2.==a2]
end

# blt[Integer, Integer]
# blt[Integer, Integer, :]

function Base.haskey(blt::BLTIndex, a1::Integer, a2::Integer, _::Colon=:)
    df = blt.all
    any(df.a1 .== a1 .&& df.a2 .== a2)
end

function Base.getindex(blt::BLTIndex, a1::Integer, a2::Integer, _::Colon=:)
    df = blt.all
    df.idx[df.a1.==a1 .&& df.a2.==a2]
end

# blt[:, :, Integer]

function Base.haskey(blt::BLTIndex, _::Colon, _::Colon, tidx::Integer)
    tidx in axes(unique(blt.all.jds), 1)
end

function Base.getindex(blt::BLTIndex, _::Colon, _::Colon, tidx::Integer)
    # It is possible for autos and crosses to have a different number of time
    # samples, so we have to bounds check tidx against both.
    gdfa = blt.autos
    idxa = tidx in axes(gdfa, 1) ? gdfa.idx : []
    gdfc = blt.crosses
    idxc = tidx in axes(gdfc, 1) ? gdfc.idx : []
    sort(vcat(idxa, idxc))
end

# blt[:, :, AbstractFloat]

function Base.haskey(blt::BLTIndex, _::Colon, _::Colon, jd::AbstractFloat)
    jd in unique(blt.all.jds)
end

function Base.getindex(blt::BLTIndex, _::Colon, _::Colon, jd::AbstractFloat)
    df = blt.all
    df.idx[df.jd .== jd]
end

# blt[Integer] - All baseline-time indices involving given ant as a1 or a2
# blt[Integer, :]

function Base.haskey(blt::BLTIndex, a::Integer, _::Colon=:)
    df = blt.all
    any(df.a1 .== a .|| df.a2 .== a)
end

function Base.getindex(blt::BLTIndex, a::Integer, _::Colon=:)
    df = blt.all
    df.idx[df.a1.==a .|| df.a2.==a]
end

# Load baseline method

function loadbaseline(vis::HDF5.Dataset, blt, a1, a2)
    mapreduce(i->vis[:,:,i], (a,b)->cat(a,b,dims=3), blt[a1, a2, :])
end

# Corrcoeff methods

function corrcoeff(vis::AbstractArray{T,3}, blt::BLTIndex, a1, a2) where T
    @assert a1 != a2 "a1 must not equal a2"
    v11 = @view vis[:, :, blt[a1,a1]]
    v22 = @view vis[:, :, blt[a2,a2]]
    v12 = haskey(blt, a1, a2) ? @view(vis[:, :, blt[a1,a2]]) :
          haskey(blt, a2, a1) ? @view(vis[:, :, blt[a2,a1]]) :
          error("no (a1,a2) or (a2,a1) baselines")
    v12 ./ sqrt.(abs.(v11) .* abs.(v22))
end

function corrcoeff(vis::AbstractArray{T,2}, blt::BLTIndex, a1, a2) where T
    @assert a1 != a2 "a1 must not equal a2"
    v11 = @view vis[:, blt[a1,a1]]
    v22 = @view vis[:, blt[a2,a2]]
    v12 = haskey(blt, a1, a2) ? @view(vis[:, blt[a1,a2]]) :
          haskey(blt, a2, a1) ? @view(vis[:, blt[a2,a1]]) :
          error("no (a1,a2) or (a2,a1) baselines")
    v12 ./ sqrt.(abs.(v11) .* abs.(v22))
end

function corrcoeff(vis::HDF5.Dataset, blt::BLTIndex, a1, a2)
    @assert a1 != a2 "a1 must not equal a2"
    v11 = loadbaseline(vis, blt, a1, a1)
    v22 = loadbaseline(vis, blt, a2, a2)
    v12 = haskey(blt, a1, a2) ? loadbaseline(vis, blt, a1, a2) :
          haskey(blt, a2, a1) ? loadbaseline(vis, blt, a2, a1) :
          error("no (a1,a2) or (a2,a1) baselines")
    @assert size(v11) == size(v12) "a1 auto size != cross size"
    @assert size(v22) == size(v12) "a2 auto size != cross size"
    v12 ./ sqrt.(abs.(v11) .* abs.(v22))
end

function corrcoeff(h5::HDF5.File, blt::BLTIndex, a1, a2)
    corrcoeff(h5["Data/visdata"], blt, a1, a2)
end

"""
Load corrcoefs for all baselines in `vis` involving `refant`.
Returns a Dict{<:Integer,<:AbstractArray}.
"""
function corrcoeff(vis::Union{<:AbstractArray{T,2},<:AbstractArray{T,3},HDF5.Dataset}, blt::BLTIndex, refant::Integer) where T
    allants = ants(blt)
    refant in allants || error("refant is not present")
    xants = filter(!=(refant), allants)
    Dict(xants .=> corrcoeff.(Ref(vis), Ref(blt), refant, xants))
end

function corrcoeff(h5::HDF5.File, blt::BLTIndex, refant::Integer)
    corrcoeff(h5["Data/visdata"], blt, refant)
end

end # module UVH5Files
