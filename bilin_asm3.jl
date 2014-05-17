macro bilin_asm3(bilin_form,mesh,nodal...)
    quote
        ##
        ## Assembly of a bilinear form using P1 elements
        ## and resulting in a sparse matrix.
        ##
        ## Usage example:
        ##    To assemble the bilinear form of the Poisson equation, write
        ##
        ##       bilin_asm(ux*uy+vx*vy,mesh)
        ##
        ##    where "mesh" is the mesh construct.
        ##    
        ##    Two additional parameters can be given: Nx1 vectors of nodal
        ##    values. During the assembly procedure, these nodal values
        ##    are then interpolated to the quadrature points and can be
        ##    used in the bilinear form using the variables "a" and "b":
        ##      
        ##      bilin_asm(a*(ux*uy+vx*vy),mesh,ones(N,1))
        ##

        # Store triplets (i,j,K_ij)
        local I = Array(Int64,0)
        local J = Array(Int64,0)
        local V = Array(Float64,0)

        # Local basis functions and their gradients
        local phi_hat = [(x,y,z)->1-x-y-z, (x,y,z)->x, (x,y,z)->y, (x,y,z)->z]
        local dphi_hat = [-1.0 -1.0 -1.0; 1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]'

        # Keast quadrature for tetrahedron (N=4?)
        local qp = [0.5854101966249685 0.1381966011250105 0.1381966011250105 0.1381966011250105;
                    0.1381966011250105 0.1381966011250105 0.1381966011250105 0.5854101966249685;
                    0.1381966011250105 0.1381966011250105 0.5854101966249685 0.1381966011250105]
        local qw = [0.25,0.25,0.25,0.25]/6
        local Nqp = size(qp,2)

        # Array for storing the pre-evaluated basis functions
        # at the quadrature points
        local phi_hat_qp = Array(Float64,4,Nqp)

        # Pre-evaluate basis functions at (local) QP's
        for q=1:Nqp
            for i=1:4
              phi_hat_qp[i,q] = phi_hat[i](qp[1,q],qp[2,q],qp[3,q])
            end
        end

        # Allocate temporary storage for
        #   A           - affine mapping
        #   invAt       - inverse-transpose of affine mapping
        #   n1,n2,n3    - node locations
        #   n1i,n2i,n3i - node indices
        #   gqp         - global quadrature points
        #   dphi        - global basis function values
        #   absdetA     - absolute value of determinant of affine map
        #   detA        - determinant of affine map
        local A = Array(Float64,3,3)
        local invAt = Array(Float64,3,3)
        local n1 = Array(Float64,3)
        local n2 = Array(Float64,3)
        local n3 = Array(Float64,3)
        local n4 = Array(Float64,3)
        local n1ix::Int64 = -1
        local n2ix::Int64 = -1
        local n3ix::Int64 = -1
        local n4ix::Int64 = -1
        local gqp = Array(Float64,3,Nqp)
        local dphi = Array(Float64,size(dphi_hat,1),size(dphi_hat,2))
        local absdetA::Float64 = -1.0
        local detA::Float64 = -1.0
        if length($nodal)==1
            local a1::Float64 = -1.0
            local a2::Float64 = -1.0
            local a3::Float64 = -1.0
            local aset = $(length(nodal)==1 ? nodal[1] : nothing)
        end
        if length($nodal)==2
            local b1::Float64 = -1.0
            local b2::Float64 = -1.0
            local b3::Float64 = -1.0
            local bset = $(length(nodal)==2 ? nodal[2] : nothing)
        end

        for k=1:size($mesh.t,2)
            # Find global node indices
            n1ix = $mesh.t[1,k]
            n2ix = $mesh.t[2,k]
            n3ix = $mesh.t[3,k]
            n4ix = $mesh.t[4,k]
            # Find the parameter values at nodes
            if length($nodal)==1
                a1 = aset[n1ix,1]
                a2 = aset[n2ix,1]
                a3 = aset[n3ix,1]
            end
            if length($nodal)==2
                b1 = bset[n1ix,1]
                b2 = bset[n2ix,1]
                b3 = bset[n3ix,1]
            end
            # Find global node locations
            n1[1] = $mesh.p[1,n1ix]
            n1[2] = $mesh.p[2,n1ix]
            n1[3] = $mesh.p[3,n1ix]
            n2[1] = $mesh.p[1,n2ix]
            n2[2] = $mesh.p[2,n2ix]
            n2[3] = $mesh.p[3,n2ix]
            n3[1] = $mesh.p[1,n3ix]
            n3[2] = $mesh.p[2,n3ix]
            n3[3] = $mesh.p[3,n3ix]
            n4[1] = $mesh.p[1,n4ix]
            n4[2] = $mesh.p[2,n4ix]
            n4[3] = $mesh.p[3,n4ix]
            # Build the affine mapping
            A[1,1] = n2[1]-n1[1]
            A[1,2] = n3[1]-n1[1]
            A[1,3] = n4[1]-n1[1]
            A[2,1] = n2[2]-n1[2]
            A[2,2] = n3[2]-n1[2]
            A[2,3] = n4[2]-n1[2]
            A[3,1] = n2[3]-n1[3]
            A[3,2] = n3[3]-n1[3]
            A[3,3] = n4[3]-n1[3]
            # Calculate global quadrature points
            for l=1:Nqp
                gqp[1,l] = A[1,1]*qp[1,l]+A[1,2]*qp[2,l]+A[1,3]*qp[3,l]+n1[1]
                gqp[2,l] = A[2,1]*qp[1,l]+A[2,2]*qp[2,l]+A[2,3]*qp[3,l]+n1[2]
                gqp[3,l] = A[3,1]*qp[1,l]+A[3,2]*qp[2,l]+A[3,3]*qp[3,l]+n1[3]
            end
            # The determinant of A and its absolute value
            detA=-A[1,3]*A[2,2]*A[3,1]+A[1,2]*A[2,3]*A[3,1]+A[1,3]*A[2,1]*A[3,2]-A[1,1]*A[2,3]*A[3,2]-A[1,2]*A[2,1]*A[3,3]+A[1,1]*A[2,2]*A[3,3]
            if detA<0
                absdetA = -detA
            else
                absdetA = detA
            end
            # Inverse of A for transforming the gradients
            invAt[1,1] = (-A[2,3]*A[3,2]+A[2,2]*A[3,3])/detA
            invAt[2,1] = (A[1,3]*A[3,2]-A[1,2]*A[3,3])/detA
            invAt[3,1] = (-A[1,3]*A[2,2]+A[1,2]*A[2,3])/detA
            invAt[1,2] = (A[2,3]*A[3,1]-A[2,1]*A[3,3])/detA
            invAt[2,2] = (-A[1,3]*A[3,1]+A[1,1]*A[3,3])/detA
            invAt[3,2] = (A[1,3]*A[2,1]-A[1,1]*A[2,3])/detA
            invAt[1,3] = (-A[2,2]*A[3,1]+A[2,1]*A[3,2])/detA
            invAt[2,3] = (A[1,2]*A[3,1]-A[1,1]*A[3,2])/detA
            invAt[3,3] = (-A[1,2]*A[2,1]+A[1,1]*A[2,2])/detA
            # Calculate global gradients
            for l=1:size(dphi_hat,2)
                dphi[1,l] = invAt[1,1]*dphi_hat[1,l]+invAt[1,2]*dphi_hat[2,l]+invAt[1,3]*dphi_hat[3,l]
                dphi[2,l] = invAt[2,1]*dphi_hat[1,l]+invAt[2,2]*dphi_hat[2,l]+invAt[2,3]*dphi_hat[3,l]
                dphi[3,l] = invAt[3,1]*dphi_hat[1,l]+invAt[3,2]*dphi_hat[2,l]+invAt[3,3]*dphi_hat[3,l]
            end

            # Loop over the local stiffness matrix
            for j=1:4
                local ux::Float64 = dphi[1,j]
                local uy::Float64 = dphi[2,j]
                local uz::Float64 = dphi[3,j]
                local jx::Int64 = $mesh.t[j,k]
                for i=1:4
                    local vx::Float64 = dphi[1,i]
                    local vy::Float64 = dphi[2,i]
                    local vz::Float64 = dphi[3,i]
                    local ix::Int64 = $mesh.t[i,k]
                    # Loop over quadrature points
                    for q=1:Nqp
                        local x = gqp[1,q]
                        local y = gqp[2,q]
                        local z = gqp[3,q]
                        local u = phi_hat_qp[j,q]
                        local v = phi_hat_qp[i,q]
                        # Calculate the values of the parameters at
                        # the quadrature points
                        if length($nodal)==1
                            a = a1+(a2-a1)*qp[1,q]+(a3-a1)*qp[2,q]
                        end
                        if length($nodal)==2
                            b = b1+(b2-b1)*qp[1,q]+(b3-b1)*qp[2,q]
                        end
                        # Local stiffness matrix values are
                        # stored in V and the corresponding
                        # global indices in I and J.
                        # Such lists can be directly feeded
                        # to the sparse matrix constructor.
                        push!(I,ix)
                        push!(J,jx)
                        push!(V,($bilin_form)*qw[q]*absdetA)
                    end
                end
            end
        end
        # Finally build the stiffness matrix
        sparse(I,J,V)
    end
end
