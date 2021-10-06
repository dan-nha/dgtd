## Description
This code is a one-dimensional discontinuous Galerkin time-domain code -- a discontinuous finite element method in time domain -- modelled after Hesthaven and Warburton \cite hesthaven2008nodal. The code differs from the given Matlab code in the way that it can be easily extended to a region based assignment of solution scheme parameters.

The code is applied to a mesh of a given 1D structure. The mesh needs to be in the current mesh format given by [Gmsh](https://gmsh.info/).


## Coding Guidelines
Though I do my best to follow the [C++ Core Guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines) by Bjarne Stroustrup and Herb Sutter, there are a few points I would like to an extra emphasis on.
* In some cases there will be a competition between conciseness, readability and performance. In this case, I will put readability (where maintainability is a part of) first, and sustainability second. I will check for performance last as this will not be an issue for a 1D code (given that I have not written mindlessly wasteful pieces of code which need to be addressed). I will not attempt to write extra clever short code, as it might be hard to read. The simpler, the better.
* I preferably write short method units to improve readability (maintainability) and aid testing. Multiple function calls are (at least to me) not an issue.
* Whenever possible and practical, I do not write methods with side effects.
* The code within the files is ordered in a specific way: The most crucial information to the reader is put first. The importance of the implemented methods to the reader descends as the file progresses.
* I distinguish between interface and implementation, where the interface will contain the public members.
* I do my best to choose descriptive names.


## Workflow
You will not see me opening any feature branches in this git repository. Obviously, because this code is only written by me, but also because I do believe that feature branching is counterproductive to CI. I will rather do a dark launch.


## C++ Standard
The standard of this code is C++20, where I mostly make use of the \<set\> library. 

## Requirements
* g++
* Boost
* BLAS
* LAPACK
* Armadillo
* ZLIB

## Installation
In the following you find a command line description of the installation. The installation is pretty standard in this case.

Clone the git repository to obtain the code

    git clone git@github.com:dan-nha/dgtd.git

Create a `build` directory in the code's root directory

    mkdir build

and go into the newly created directory

    cd build

Afterwards type the commands

    cmake ..
and

    make install

into your command line.


## Usage
In the code's root directory you will find an binary file called `dgtd` lying in the folder `bin`. Execute the binary by giving it a mesh file produced by *Gmsh* as argument, e.g.

  ./dgtd my.msh

You can find example meshes under `examples`.


## Documentation
I rather give some motivation for certain design choices, and the context
of certain methods than describe what a part of code is doing.  In this
regard, I hope the given names and the implementation to be verbose enough
to be able to abstain from lengthy documentation. This might not always be
possible though.

To build the documentation go to the `doc` directory and type the following command.

    doxygen Doxyfile.in

Open `doc/doxygen/html/index.html` with a browser to view the documentation.
