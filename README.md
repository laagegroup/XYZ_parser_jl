# XYZ.jl
A Julia package for reading XYZ trajectory files.

This package is in an early stage. It currently supports XYZ files with a constant number of atoms. This package is a copy of the DCD module in my gitlab.

## Installation (with set-up ssh keys)

````julia
pkg> add git@gitlab.com:axel.gomez/xyz-lammps-julia.git
````

## Usage

````julia
using XYZ

xyz = load_xyz("foo.xyz")

na = natoms(xyz)  # Get the number of atoms.
nf = nframes(xyz) # Get the number of frames.

for frame in xyz  # Iterate over each frame.
  r = positions(frame)  # Get an array of all current positions with dimensions (3, na).
  
  # Do stuff.
end
````

Happy coding!
