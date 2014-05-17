# Julia finite elements

This is a simple two- and three-dimensional P1 finite element method implementation in Julia. Current features are very limited:

* Simple triangle refinement
* Assembly of bilinear forms
* Assembly of linear forms
* Assembly of boundary (bi)linear forms

* Visualization using MATLAB (not included, software is commerical).
* 3D mesh generation using tetgen through iso2mesh (not included, tetgen is free for research purposes).

## Installation

1. Get the newest version of Julia. I develop this using Win7 OS and haven't yet tested it with other OS'es.
2. Install packages PyCall, PyPlot, MATLAB.jl using the command Pkg.add().
3. Make sure that MATLAB is installed (tested with R2013).
4. Add *MATLAB_HOME* system variable for MATLAB.jl which points to MATLAB R2013 base directory.
5. Install iso2mesh 2013 and add all the files under "iso2mesh" subdirectory.
6. Open Julia and go to the base directory, run include("envglue.jl") and then Fenv.init(). The MATLAB prompt should open.
7. Start playing around with example programs in examples.jl (for 2D) and examples3.jl (for 3D).

## Info

The project is something I plan to advance during weekends.
My goal is to learn more about finite element solvers
and the underlying mathematics. Check also my blog
at http://derpsince.blogspot.fi for updates.
