#
# Simple finite element assembly 3D
# examples using piecewise-linear tetras
#
# Author: Tom Gustafsson, 17.5.2014
# Licensed under GPLv3
#

function test_poisson3()
    ##
    ## Test case solves the Poisson equation
    ##   u_xx + u_yy + u_zz = 1,
    ## with the boundary condition
    ##   u = 0,
    ## on unit cube.
    ##
    mesh = Mesh3(0.0001)
    S = @bilin_asm3(ux*vx+uy*vy+uz*vz,mesh)
    f = @lin_asm3(v,mesh)
    # Find Dirichlet node indices
    D = Array(Int64,0)
    for i=1:size(mesh.p,2)
        if mesh.p[1,i]==0 || mesh.p[1,i]==1 || mesh.p[2,i]==0 || mesh.p[2,i]==1 || mesh.p[3,i]==0 || mesh.p[3,i]==1
            push!(D,i)
        end
    end
    I = setdiff(1:size(mesh.p,2),D)
    # Solve the system
    u = Array(Float64,size(mesh.p,2))
    u = zeros(size(mesh.p,2))
    #u[D] = sin(mesh.p[1,D])
    #f = ones(size(mesh.p,2),1)
    #u[I]=cholfact(S[I,I])\(-S[I,D]*u[D])
    u[I]=cholfact(S[I,I])\(f[I,1])

    tetmesh(mesh)
    tetsurf(mesh,u,"z<0.5")


    # Check that the maximum value is correct
    1
end
