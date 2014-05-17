##
## Some mesh tools for finite elements.
##
## Author: Tom Gustafsson, 27.4.2014
## Licensed under GPLv3
##

# Allows to use the plotting commands
# of matplotlib
using PyPlot

type Mesh
    ## Mesh type for storing
    ## triangular two-dimesional
    ## finite element meshes

    # "TRIMESH"-structure (by M. Juntunen & A. Hannukainen):
    # * p     = node coordinates in 2xN -matrix
    # * t     = triangles in 3xN -matrix, contains indices to p
    # * edges = matrix of all edges in the mesh.
    #           each column is an edge [n1;n2] with n1<n2.
    # * t2e   = connects triangles and edges.
    #           each column corresponds to a triangle and has
    #           the triangle's edges in the order n1->n2, n2->n3, n1->n3
    # * e2t   = vice versa
    p::Array{Float64,2}
    t::Array{Int64,2}
    edges::Array{Int64,2}
    t2e::Array{Int64,2}
    e2t::Array{Int64,2}

    function Mesh(p,t)
        ## Constructor to build the edge
        ## numbering etc

        # Sort the triangle indices column-wise
        t = sort(t,1)
        # Allocate space for edges and e2t-mapping.
        # Must be truncated in the end to correct size.
        edges = Array(Int64,2,3*size(t,2))
        e2t = zeros(2,3*size(t,2))
        # Allocate space for t2e-mapping and for
        # a map-structure for checking whether
        # each edge has already been added.
        t2e = zeros(3,size(t,2))
        check = Dict{(Int64,Int64),Int64}()
        # Compute total amount of edges to N
        N = 1;
        for k=1:size(t,2)
            n1=t[1,k]
            n2=t[2,k]
            n3=t[3,k]
            t2e1 = 0
            t2e2 = 0
            t2e3 = 0
            
            if !haskey(check,(n1,n2))
                edges[1,N]=n1;
                edges[2,N]=n2;
                t2e1 = N;
                check[(n1,n2)]=N
                N += 1;
            else
                t2e1 = check[(n1,n2)]
            end
            if !haskey(check,(n2,n3))
                edges[1,N]=n2;
                edges[2,N]=n3;
                t2e2 = N;
                check[(n2,n3)]=N
                N += 1;
            else
                t2e2 = check[(n2,n3)]
            end
            if !haskey(check,(n1,n3))
                edges[1,N]=n1;
                edges[2,N]=n3;
                t2e3 = N;
                check[(n1,n3)]=N
                N += 1;
            else
                t2e3 = check[(n1,n3)]
            end
            t2e[1,k] = t2e1
            t2e[2,k] = t2e2
            t2e[3,k] = t2e3
            # Build the inverse mappings
            if e2t[1,t2e1]==0
                e2t[1,t2e1]=k
            else
                e2t[2,t2e1]=k
            end
            if e2t[1,t2e2]==0
                e2t[1,t2e2]=k
            else
                e2t[2,t2e2]=k
            end
            if e2t[1,t2e3]==0
                e2t[1,t2e3]=k
            else
                e2t[2,t2e3]=k
            end
        end
        new(p,t,edges[:,1:(N-1)],t2e,e2t[:,1:(N-1)])
    end
end

# Default constructor
Mesh() = Mesh([0 0; 1 0; 0 1; 1 1]',[1 2 3; 2 3 4]')

function trimesh(mesh)
    # Plot a finite element mesh
    figure()
    tx = Array(Float64,size(mesh.t,1)+1,size(mesh.t,2))
    ty = Array(Float64,size(mesh.t,1)+1,size(mesh.t,2))
    for k=1:size(mesh.t,2)
        tx[:,k] = mesh.p'[mesh.t[[1,2,3,1],k],1]
        ty[:,k] = mesh.p'[mesh.t[[1,2,3,1],k],2]
    end
    plot(tx,ty,"k-")
end

function triplot(mesh,u,N=10)
    # Draw a contour plot of "u" on "mesh"
    figure()
    tricontour(mesh.p'[:,1],mesh.p'[:,2],mesh.t'-1,u,N)
end

function trisurf(mesh,u)
    figure()
    # TODO figure a way to add the mesh and some
    # colormap here.
    plot_trisurf(mesh.p'[:,1],mesh.p'[:,2],u,cmap=ColorMap("jet"))
end

function refinetri(mesh::Mesh)
    Np = size(mesh.p,2)
    Ne = size(mesh.edges,2)
    Nt = size(mesh.t,2)
    # New points
    rp = Array(Float64,2,Np+Ne)
    for i=1:Np
        rp[1,i] = mesh.p[1,i]
        rp[2,i] = mesh.p[2,i]
    end
    for j=1:Ne
        rp[1,Np+j] = 0.5*(mesh.p[1,mesh.edges[1,j]]+mesh.p[1,mesh.edges[2,j]])
        rp[2,Np+j] = 0.5*(mesh.p[2,mesh.edges[1,j]]+mesh.p[2,mesh.edges[2,j]])
    end
    # New triangles
    rt = Array(Int64,3,4*Nt)
    for k=1:Nt
        rt[1,k] = mesh.t[1,k]
        rt[1,k+Nt] = mesh.t[2,k]
        rt[1,k+2*Nt] = mesh.t[3,k]
        rt[1,k+3*Nt] = mesh.t2e[1,k]+Np
        rt[2,k] = mesh.t2e[1,k]+Np;
        rt[2,k+Nt] = mesh.t2e[1,k]+Np;
        rt[2,k+2*Nt] = mesh.t2e[3,k]+Np;
        rt[2,k+3*Nt] = mesh.t2e[2,k]+Np;
        rt[3,k] = mesh.t2e[3,k]+Np;
        rt[3,k+Nt] = mesh.t2e[2,k]+Np;
        rt[3,k+2*Nt] = mesh.t2e[2,k]+Np;
        rt[3,k+3*Nt] = mesh.t2e[3,k]+Np;
    end

    return Mesh(rp,rt)
end

function meshtest()
    mesh = Mesh()
    for i=1:4
        mesh = refinetri(mesh)
    end
    f(x,y) = x.*y
    triplot(mesh,f(mesh.p'[:,1],mesh.p'[:,2]),20)
end
