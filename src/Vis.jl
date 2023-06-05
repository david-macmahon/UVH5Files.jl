struct Vis <: BLTData
    vis::HDF5.Dataset
    blts::BLTIndices
end

function Vis(h5::HDF5.File, blts::BLTIndices)
    Vis(h5["Data/visdata"], blts)
end

function bltindices(v::Vis)
    v.blts
end

function Base.getindex(v::Vis, bltidxs::AbstractArray)
    mapreduce(i->v.vis[:,:,i], (a,b)->cat(a,b,dims=3), bltidxs)
end