struct Vis <: BLTData
    vis::HDF5.Dataset
    blts::BLTIndices
end

# We are a BLTData subtype so we have to implement `bltindices`
function bltindices(v::Vis)
    v.blts
end

# The BLTData getindex methods will eventually call this, if needed
function getindex(v::Vis, bltidxs::AbstractArray)
    mapreduce(i->v.vis[:,:,i], (a,b)->cat(a,b,dims=3), bltidxs)
end
