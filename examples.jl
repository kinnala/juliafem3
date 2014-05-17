#
# Simple finite element assembly
# using piecewise linear basis functions
# on triangular mesh.
#
# Author: Tom Gustafsson, 12.4.2014
# Licensed under GPLv3
#

function test_poisson()
    ##
    ## Test case solves the Poisson equation
    ##   u_xx + u_yy = 1,
    ## with the boundary condition
    ##   u = 0,
    ## on unit square.
    ##
    mesh = Mesh()
    for i=1:4
        mesh = refinetri(mesh)
    end
    S = @bilin_asm(ux*vx+uy*vy,mesh)
    f = @lin_asm(v,mesh)
    # Find Dirichlet node indices
    D = Array(Int64,0)
    for i=1:size(mesh.p,2)
        if mesh.p[1,i]==0 || mesh.p[1,i]==1 || mesh.p[2,i]==0 || mesh.p[2,i]==1
            push!(D,i)
        end
    end
    I = setdiff(1:size(mesh.p,2),D)
    # Solve the system
    u = Array(Float64,size(mesh.p,2))
    u = zeros(size(mesh.p,2))
    u[I]=cholfact(S[I,I])\f[I,1]

    # Check that the maximum value is correct
    exactmaxu = 0.07344576676891949
    eps = 0.00001
    if maximum(u)>exactmaxu+eps || maximum(u)<exactmaxu-eps
        return false
    end
    true
end

function test_nonlinear_poisson()
    ##
    ## This test case solves the non-linear Poisson problem
    ## related to the minimal surface area problem.
    ##
    mesh = Mesh()
    for i=1:4
        mesh = refinetri(mesh)
    end

    # Find Dirichlet node indices
    D = Array(Int64,0)
    for i=1:size(mesh.p,2)
        if mesh.p[1,i]==0 || mesh.p[1,i]==1 || mesh.p[2,i]==0 || mesh.p[2,i]==1
            push!(D,i)
        end
    end
    I = setdiff(1:size(mesh.p,2),D)

    # Solution vector
    U = ones(size(mesh.p,2),1)

    ffun = (x,y)->sin(2*pi*x-0.7).*sin(2*pi*y-0.3)

    U=ffun(mesh.p[1,:]',mesh.p[2,:]')

    # Gradient calc
    DX = @bilin_asm(ux*v,mesh)
    DY = @bilin_asm(uy*v,mesh)
    M = @bilin_asm(u*v,mesh)
    CM = cholfact(M)
    dx = w->CM\(DX*w)
    dy = w->CM\(DY*w)

    for i=1:3
        X = 1.0/sqrt(1.0+dx(U).^2+dy(U).^2)
        S = @bilin_asm(a*(ux*vx+uy*vy),mesh,X)
        #Solve the linear system
        U[I,1]=cholfact(S[I,I])\(-S[I,D]*U[D,1])
    end

    exactsumU = -13.562611418335702
    eps = 0.00001
    if sum(U)>exactsumU+eps || sum(U)<exactsumU-eps
        return false
    end
    true
end

function test_poisson_asm_time()
    ##
    ## Test case solves the Poisson equation
    ##   u_xx + u_yy = 1,
    ## with the boundary condition
    ##   u = 0,
    ## on unit square.
    ##
    mesh = Mesh()
    for i=1:8
        mesh = refinetri(mesh)
    end
    N = size(mesh.p,2)
    @time S = @bilin_asm(ux*vx+uy*vy,mesh)
    @time f = @lin_asm(v,mesh)
    @time Sa = @bilin_asm(a*(ux*vx+uy*vy),mesh,ones(N,1))
    @time fa = @lin_asm(v,mesh,ones(N,1))
    print("The benchmark times as of 27.4.2014 are: 0.79, 0.58, 2.76, 0.98.")
    1
end


function runtests()
    print("test_poisson...")
    print(test_poisson())
    print("\n")
    print("test_nonlinear_poisson...")
    print(test_nonlinear_poisson())
    print("\n")
end

function demo_nonlinear_poisson()
    mesh = Mesh()
    for i=1:4
        mesh = refinetri(mesh)
    end

    # Find Dirichlet node indices
    D = Array(Int64,0)
    for i=1:size(mesh.p,2)
        if mesh.p[1,i]==0 || mesh.p[1,i]==1 || mesh.p[2,i]==0 || mesh.p[2,i]==1
            push!(D,i)
        end
    end
    I = setdiff(1:size(mesh.p,2),D)

    # Solution vector
    U = ones(size(mesh.p,2),1)

    ffun = (x,y)->sin(2*pi*x-0.7).*sin(2*pi*y-0.3)

    U=ffun(mesh.p[1,:]',mesh.p[2,:]')

    # Gradient calc
    DX = @bilin_asm(ux*v,mesh)
    DY = @bilin_asm(uy*v,mesh)
    M = @bilin_asm(u*v,mesh)
    CM = cholfact(M)
    dx = w->CM\(DX*w)
    dy = w->CM\(DY*w)

    for i=1:3
        #trisurf(mesh,U[:,1])
        X = 1.0/sqrt(1.0+dx(U).^2+dy(U).^2)
        S = @bilin_asm(a*(ux*vx+uy*vy),mesh,X)
        #Solve the linear system
        U[I,1]=cholfact(S[I,I])\(-S[I,D]*U[D,1])
    end

    #trisurf(mesh,U[:,1])
    @mput mesh
    @mput U
    @matlab trisurf(mesh.t',mesh.p(1,:),mesh.p(2,:),U)
end

function play()
    mesh = Mesh()
    #for i=1:8
    #    mesh = refinetri(mesh)
    #end
    @time M = @bndbilin_asm(u*v,mesh)
    @time T = @bndlin_asm(x*v,mesh)
    1
end
