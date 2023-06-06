module UVH5Files

using HDF5
using DataFrames

export UVH5File, corrcoeff, loadbaseline

include("BLTData.jl")
include("BLTIndices.jl")
include("Vis.jl")

struct UVH5File
    h5::HDF5.File
    blts::BLTIndices
    vis::Vis

    function UVH5File(h5::HDF5.File)
        blts = BLTIndices(h5)
        vis = Vis(h5, blts)
        new(h5, blts, vis)
    end
end

function UVH5File(h5obj::Union{HDF5.Attribute, HDF5.Dataset,
                               HDF5.Datatype, HDF5.Group})
    UVH5File(HDF5.file(h5obj))
end

function UVH5File(h5name::AbstractString)
    UVH5File(h5open(h5name))
end

function Base.show(io::IO, uv::UVH5File)
    blts = uv.blts
    nants = length(ants(blts))
    nautos = nrow(parent(blts.autos))
    ncross = nrow(parent(blts.crosses))
    ntautos = length(blts.autos)
    ntcross = length(blts.crosses)
    print(io, "UVH5File($(basename(HDF5.filename(uv.h5))): " *
              "$nants ants, $nautos autos/$ntautos times, " *
              "$ncross crosses/$ntcross times)")
end

function Base.close(uv::UVH5File)
    close(uv.h5)
end

# Stuff below here probably belongs elsewhere...

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
