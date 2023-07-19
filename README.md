<div id="top"></div>

<!-- PROJECT SHIELDS -->

[![GNU AGPL v3.0 License][license-shield]][license-url]
![VERSION](https://img.shields.io/badge/version-0.1.2-blue)

<!-- TABLE OF CONTENTS -->

<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about">XYZ.jl</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation (with set-up ssh keys)</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#acknowledgments">License and acknowledgments</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->
<div id="about"></div>

## XYZ.jl

A Julia package for reading XYZ trajectory files.

It currently supports XYZ files with a constant number of atoms.

<p align="right">(<a href="#top">back to top</a>)</p>

<div id="getting-started"></div>

## Getting Started

<div id="prerequisites"></div>

### Prerequisites

This package is a Julia parser for XYZ files. It works with Julia 1.5 or above.

It requires the DCD package. You should install the DCD parser first!

<p align="right">(<a href="#top">back to top</a>)</p>

<div id="installation"></div>

## Installation (with set-up ssh keys)

Adding private packages to Julia requires a password-less access (see the [manual in this organization for setup](https://github.com/laagegroup/0_HowTo/blob/main/Github_beginner_guide.md#setup-a-password-less-access-over-ssh))

````julia
pkg> add git@github.com:laagegroup/XYZ_parser_jl.git
````
<p align="right">(<a href="#top">back to top</a>)</p>

<div id="usage"></div>

## Installation (with https, RECOMMENDED)

This package can be directly added via https by running the following command from the julia interface:

````julia
pkg> add https://github.com/laagegroup/XYZ_parser_jl.git
````
<p align="right">(<a href="#top">back to top</a>)</p>

<div id="usage"></div>


## Usage

This package can be used in the Julia REPL or within a Julia script. Here is an example of the main features:

````julia
using XYZ

xyz = load_xyz("foo.xyz")

na = natoms(xyz)  # Get the number of atoms.
nf = nframes(xyz) # Get the number of frames.
atoms = atomnames(xyz) # Get the list of atom names.

for frame in xyz  # Iterate over each frame.
  r = positions(frame)  # Get an array of all current positions with dimensions (3, na).

  # Work on the current frame.
end
````

<p align="right">(<a href="#top">back to top</a>)</p>


<div id="acknowledgments"></div>

## License and acknowledgments

This package is a modification of the DCD parser in this organization: https://github.com/laagegroup/DCD_parser_jl.git (The original development is due to Michael von Domaros, see https://github.com/mvondomaros/DCD.jl.git and the MIT license `original_MIT_license`).

Modifications are distributed under the GNU Affero General Public License v3.0. See `LICENSE` for more information.

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[license-shield]: https://img.shields.io/github/license/laagegroup/XYZ_parser_jl.svg?style=for-the-badge
[license-url]: https://github.com/laagegroup/XYZ_parser_jl/blob/main/LICENSE
