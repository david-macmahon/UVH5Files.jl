struct BLTIndices
    all::DataFrame
    autos::GroupedDataFrame
    crosses::GroupedDataFrame

    function BLTIndices(a1s::AbstractVector{<:Integer},
                        a2s::AbstractVector{<:Integer},
                        jds::AbstractVector{Float64})
        # Get all baseline-time-idx named tuples
        all = [(; a1, a2, jd, idx) for (idx, (a1, a2, jd)) in enumerate(zip(a1s, a2s, jds))] |> DataFrame;
        # Group auto and cross subsets by JD
        autos   = groupby(subset(all, [:a1,:a2]=>ByRow((a1,a2)->a1==a2), view=true), :jd)
        crosses = groupby(subset(all, [:a1,:a2]=>ByRow((a1,a2)->a1!=a2), view=true), :jd)

        new(all, autos, crosses)
    end
end

function BLTIndices(h5obj::Union{HDF5.Attribute, HDF5.Dataset, HDF5.Datatype, HDF5.Group})
    BLTIndices(HDF5.file(h5obj))
end

function ants(blts::BLTIndices)
    sort(blts.all.a1 ∪ blts.all.a2)
end

import Base.==
function ==(a::BLTIndices, b::BLTIndices)
    a.all == b.all
end

function Base.show(io::IO, blts::BLTIndices)
    nants = length(ants(blts))
    nautos = nrow(parent(blts.autos))
    ncross = nrow(parent(blts.crosses))
    ntautos = length(blts.autos)
    ntcross = length(blts.crosses)
    print(io, "BLTIndices($nants ants, $nautos autos/$ntautos times, $ncross crosses/$ntcross times)")
end

# blts[[idxs...]]

function getindex(blts::BLTIndices, bltidxs::AbstractArray)
    blts.all[bltidxs, :]
end

# blts[Integer, Integer, Integer]

function haskey(blts::BLTIndices, a1::Integer, a2::Integer, tidx::Integer)
    gdf  = (a1 == a2) ? blts.autos : blts.crosses
    tidx in axes(gdf, 1) || return false
    any(gdf[tidx].a1 .== a1 .&& gdf[tidx].a2 .== a2)
end

function getindex(blts::BLTIndices, a1::Integer, a2::Integer, tidx::Integer)
    sdf = (a1 == a2) ? blts.autos[tidx] : blts.crosses[tidx]
    sdf.idx[sdf.a1.==a1 .&& sdf.a2.==a2]
end

# blts[Integer, Integer, AbstractFloat]

function haskey(blts::BLTIndices, a1::Integer, a2::Integer, jd::AbstractFloat)
    gdf  = (a1 == a2) ? blts.autos : blts.crosses
    haskey(gdf, (jd,)) || return false
    any(gdf[(jd,)].a1 .== a1 .&& gdf[(jd,)].a2 .== a2)
end

function getindex(blts::BLTIndices, a1::Integer, a2::Integer, jd::AbstractFloat)
    sdf = (a1 == a2) ? blts.autos[(jd,)] : blts.crosses[(jd,)]
    sdf.idx[sdf.a1.==a1 .&& sdf.a2.==a2]
end

# blts[Integer, Integer]
# blts[Integer, Integer, :]

function haskey(blts::BLTIndices, a1::Integer, a2::Integer, _::Colon=:)
    df = blts.all
    any(df.a1 .== a1 .&& df.a2 .== a2)
end

function getindex(blts::BLTIndices, a1::Integer, a2::Integer, _::Colon=:)
    df = blts.all
    df.idx[df.a1.==a1 .&& df.a2.==a2]
end

# blts[:, :, Integer]

function haskey(blts::BLTIndices, _::Colon, _::Colon, tidx::Integer)
    tidx in axes(unique(blts.all.jds), 1)
end

function getindex(blts::BLTIndices, _::Colon, _::Colon, tidx::Integer)
    # It is possible for autos and crosses to have a different number of time
    # samples, so we have to bounds check tidx against both.
    gdfa = blts.autos
    idxa = tidx in axes(gdfa, 1) ? gdfa.idx : []
    gdfc = blts.crosses
    idxc = tidx in axes(gdfc, 1) ? gdfc.idx : []
    sort(vcat(idxa, idxc))
end

# blts[:, :, AbstractFloat]

function haskey(blts::BLTIndices, _::Colon, _::Colon, jd::AbstractFloat)
    jd in unique(blts.all.jds)
end

function getindex(blts::BLTIndices, _::Colon, _::Colon, jd::AbstractFloat)
    df = blts.all
    df.idx[df.jd .== jd]
end

# blts[Integer] - All baseline-time indices involving given ant as a1 or a2
# blts[Integer, :]

function haskey(blts::BLTIndices, a::Integer, _::Colon=:)
    df = blts.all
    any(df.a1 .== a .|| df.a2 .== a)
end

function getindex(blts::BLTIndices, a::Integer, _::Colon=:)
    df = blts.all
    df.idx[df.a1.==a .|| df.a2.==a]
end

# blts[:, Integer] - All auto-correlation baseline-time indices for given time index

function haskey(blts::BLTIndices, _::Colon, tidx::Integer)
    haskey(blts.autos, tidx)
end

function getindex(blts::BLTIndices, _::Colon, tidx::Integer)
    blts.autos[tidx].idx
end
