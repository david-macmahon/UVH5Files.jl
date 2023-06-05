"""
Abstract type for HDF5 datasets that have a basline-time index.  Mainly used for
the large datasets in "/Data" that we may not want to load into memory in their
entirety.
"""
abstract type BLTData end

"""
All sub-types of BLTData must implement a type-specific method for
`bltindices()` that returns the BLTIndices corresponding to the BLTData object.
"""
function bltindices(btl::BLTData)
    error("type $(typeof(btl)) has no `bltindices()` method")
end

function Base.haskey(bltdata::BLTData, a1::Integer, a2::Integer, tidx::Integer)
    haskey(bltindices(bltdata), a1, a2, tidx)
end

function Base.getindex(bltdata::BLTData, a1::Integer, a2::Integer, tidx::Integer)
    uv = bltindices(bltdata)
    idxs = uv[a1, a2, tidx]
    getindex(bltdata, idxs)
end

function Base.haskey(bltdata::BLTData, a1::Integer, a2::Integer, jd::AbstractFloat)
    haskey(bltindices(bltdata), a1, a2, jd)
end

function Base.getindex(bltdata::BLTData, a1::Integer, a2::Integer, jd::AbstractFloat)
    uv = bltindices(bltdata)
    idxs = uv[a1, a2, jd]
    getindex(bltdata, idxs)
end

function Base.haskey(bltdata::BLTData, a1::Integer, a2::Integer, _::Colon=:)
    haskey(bltindices(bltdata), a1, a2, :)
end

function Base.getindex(bltdata::BLTData, a1::Integer, a2::Integer, _::Colon=:)
    uv = bltindices(bltdata)
    idxs = uv[a1, a2, :]
    getindex(bltdata, idxs)
end

function Base.haskey(bltdata::BLTData, _::Colon, _::Colon, tidx::Integer)
    haskey(bltindices(bltdata), :, :, tidx)
end

function Base.getindex(bltdata::BLTData, _::Colon, _::Colon, tidx::Integer)
    uv = bltindices(bltdata)
    idxs = uv[:, :, tidx]
    getindex(bltdata, idxs)
end

function Base.haskey(bltdata::BLTData, _::Colon, _::Colon, jd::AbstractFloat)
    haskey(bltindices(bltdata), :, :, jd)
end

function Base.getindex(bltdata::BLTData, _::Colon, _::Colon, jd::AbstractFloat)
    uv = bltindices(bltdata)
    idxs = uv[:, :, jd]
    getindex(bltdata, idxs)
end

function Base.haskey(bltdata::BLTData, a::Integer, _::Colon=:)
    haskey(bltindices(bltdata), a, :)
end

function Base.getindex(bltdata::BLTData, a::Integer, _::Colon=:)
    uv = bltindices(bltdata)
    idxs = uv[a, :]
    getindex(bltdata, idxs)
end

function Base.haskey(bltdata::BLTData, _::Colon, tidx::Integer)
    haskey(bltindices(bltdata), :, tidx)
end

function Base.getindex(bltdata::BLTData, _::Colon, tidx::Integer)
    uv = bltindices(bltdata)
    idxs = uv[:, tidx]
    getindex(bltdata, idxs)
end