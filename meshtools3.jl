##
## Some 3D mesh tools for finite elements.
## Utilizes iso2mesh.
##
## Author: Tom Gustafsson, 17.5.2014
## Licensed under GPLv3
##

type Mesh3
    ## Mesh3 type for storing
    ## tetrahedral three-dimesional
    ## finite element meshes

    # "TRIMESH"-structure:
    # * p     = node coordinates in 3xN -matrix
    # * f     = faces of tets in 3xN -matrix
    # * t     = tets in 4xN -matrix, contains indices to p
    p::Array{Float64,2}
    f::Array{Int64,2}
    t::Array{Int64,2}

    function Mesh3(a) # Utilizes iso2mesh
        @mput a
        eval_string("[no,fa,el]=meshabox([0 0 0],[1 1 1],a);")
        @mget no fa el 
        new(no',fa',el')
    end
    function Mesh3(p,f,t)
        new(p,f,t)
    end
end

Mesh() = Mesh(1)

function tetmesh(mesh)
    # Plot a finite element tetra mesh (MATLAB)
    mesh
    @mput mesh
    @matlab figure()
    @matlab plotmesh(mesh.p',mesh.f',mesh.t')
end

function tetsurf(mesh,u,setup)
    @mput setup
    @mput mesh
    @mput u
    @matlab figure()
    @matlab plotmesh([mesh.p' u],mesh.f',mesh.t',setup)
    @matlab begin
        shading("interp")
    end
end

function tetsurf(mesh,u)
    @mput mesh
    @mput u
    @matlab figure()
    @matlab plotmesh([mesh.p' u],mesh.f',mesh.t')
    @matlab begin
        shading("interp")
    end
end
