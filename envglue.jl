module Fenv
##
## Finite elements Julia+MATLAB environment
##
## Author: Tom Gustafsson, 17.5.2014
## Licensed under GPLv3
##

using MATLAB
using PyPlot

function init()
    print("Opening MATLAB session... ")
    restart_default_msession()
    print("Changing MATLAB PWD... ")
    matlabchangedir = pwd()
    @mput matlabchangedir
    @matlab cd(matlabchangedir)
    @matlab addpath("iso2mesh")
    print("Ok!\n")
    # Clear MATLAB workspace
    @matlab clear
end

# P1 tools for 2D meshes using Julia
include("meshtools.jl")
include("bilin_asm.jl")
include("lin_asm.jl")
include("bndbilin_asm.jl")
include("bndlin_asm.jl")
include("examples.jl") # Depends on previous

# P1 tools for 3D meshes using Julia+MATLAB
include("meshtools3.jl")
include("bilin_asm3.jl") 
include("lin_asm3.jl") 
include("examples3.jl") # Depends on previous

end
