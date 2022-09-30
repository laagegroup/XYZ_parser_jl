module XYZ

using Parameters
import DCD: positions

export load_xyz, natoms, nframes, positions, atomnames


"""
    load_xyz(f::String)

Load a xyz trajectory from file `f`.
"""
function load_xyz(f::String)
    io = open(f)
    load_xyz(io)
end

"""
    load_xyz(io::IO)

Load a xyz trajectory from stream `io`.
"""
function load_xyz(io::IO)
    FileXYZ(io)
end

"""
    @with_kw struct FileXYZ

Stores information about a xyz trajectory file.
"""
@with_kw struct FileXYZ
    io::IO
    natoms::Int64
    nframes::Int64
    atomnames::Array{String,1}
    posbytes::Array{Int64,1}
end

"""
    FileXYZ(io::IO)

Construct a [`XYZ.FileXYZ`](@ref) object from an IO stream `io`.
"""
function FileXYZ(io::IO)
    loader = Loader(io)
    FileXYZ(
        io=io,
        natoms=natoms(loader),
        nframes=nframes(loader),
        atomnames=atomnames(loader),
        posbytes=posbytes(loader),
        )
end

"""
    iostream(f::FileXYZ)

Get the [`IOStream`](@ref) associated with the [`XYZ.FileXYZ`](@ref) object.
"""
iostream(f::FileXYZ) = f.io

"""
    natoms(f::FileXYZ)

Get the number of atoms.
"""
natoms(f::FileXYZ) = f.natoms

"""
    nframes(f::FileXYZ)

Get the number of frames.
"""
nframes(f::FileXYZ) = f.nframes

"""
    atomnames(f::FileXYZ)

Get all atom names, the topology, of a frame.
"""
atomnames(f::FileXYZ) = f.atomnames

"""
    posbytes(f::FileXYZ)

Get the position of each frame in a Array{Int64,1}.
"""
posbytes(f::FileXYZ) = f.posbytes

"""
    Base.eltype(::Type{FileXYZ})

[`XYZ.FileXYZ`](@ref) objects are iterators over [`XYZ.FrameXYZ`](@ref) objects.
"""
Base.eltype(::Type{FileXYZ}) = FrameXYZ

function Base.getindex(f::FileXYZ, i::Int64)
    1 <= i <= f.nframes || throw(BoundsError(f, i))
    return FrameXYZ(f, i)
end
Base.getindex(f::FileXYZ, i::Number) = f[convert(Int64, i)]
Base.getindex(f::FileXYZ, I) = [f[i] for i in I]
Base.firstindex(f::FileXYZ) = 1
Base.lastindex(f::FileXYZ) = f.nframes
Base.length(f::FileXYZ) = f.nframes

function Base.iterate(f::FileXYZ, frame::Int64=1)
    if frame > nframes(f)
        return nothing
    else
        return FrameXYZ(f, frame), frame + 1
    end
end


"""
    @with_kw struct Loader

Stores information from the xyz trajectory file pre-loading.
"""
@with_kw struct Loader
    nframes::Int64
    firststep::Int64
    laststep::Int64
    natoms::Int64
    atomnames::Array{String,1}
    posbytes::Array{Int64,1}
end

"""
    loader(io::IO)

Construct a [`XYZ.Loader`](@ref) object from an [`IOStream`](@ref).
"""
function Loader(io::IO)
    firststep = 1
    natoms = parse(Int64,readline(io))
    readline(io)
    posbytes = Int64[position(io)]
    atomnames = String[]
    for iat=1:natoms
        l = readline(io)
        push!(atomnames,string(split(l)[1]))
    end
    while !eof(io)
        readline(io)
        readline(io)
        append!(posbytes,position(io))
        for iat=1:natoms
            readline(io)
        end
    end
    nframes = length(posbytes)
    laststep = nframes

    Loader(
        nframes=nframes, 
        firststep=firststep, 
        laststep=laststep, 
        natoms=natoms,
        atomnames=atomnames,
        posbytes=posbytes
        )
end

"""
    natoms(h::Loader)

Get the number of atoms.
"""
natoms(h::Loader)::Int64 = h.natoms

"""
    nframes(h::Loader)

Get the number of frames.
"""
nframes(h::Loader)::Int64 = h.nframes

"""
    atomnames(h::Loader)

Get all atom names, the topology, of a frame.
"""
atomnames(h::Loader)::Array{String,1} = h.atomnames

"""
    posbytes(h::Loader)

Get the position of each frame in a Array{Int64,1}.
"""
posbytes(h::Loader)::Array{Int64,1} = h.posbytes

"""
    seekframe(f::FileXYZ, index::Int64)

Move the file's [`IOStream`](@ref) to the position of the indicated frame.
"""
function seekframe(f::FileXYZ, index::Int64)
    1 <= index <= f.nframes || throw(BoundsError(f, index))
    pos = posbytes(f)[index]
    seek(f.io, pos)
end


"""
    @with_kw struct FrameXYZ

Stores positions.
"""
@with_kw struct FrameXYZ
    positions::Array{Float64,2}
end

"""
    FrameXYZ(f::FileXYZ, index::Int64)

Construct a [`XYZ.FrameXYZ`](@ref) object from a [`XYZ.FileXYZ`](@ref) `f` with the specified `index`.
"""
function FrameXYZ(f::FileXYZ, index::Int64)
    seekframe(f, index)
    io = iostream(f)
    positions = Array{Float64,2}(undef, 3, natoms(f))
    for iatom in 1:natoms(f)
        line = split(readline(io))[2:end]
        for dim in 1:3
            positions[dim, iatom] = parse(Float64,line[dim])
        end
    end

    FrameXYZ(positions=positions)
end

"""
    positions(f::FrameXYZ)

Get the current positions.
"""
positions(f::FrameXYZ) = f.positions

end
