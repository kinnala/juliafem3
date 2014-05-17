macro bndlin_asm(bndlin_form,mesh,nodal...)
    quote
        ## TODO Add info
        ## TODO add support for nodal values

        # Store triplets (i,j,K_ij)
        local I = Array(Int64,0)
        local V = Array(Float64,0)

        # Local basis functions and their gradients
        local phi_hat = [(x)->1.0-x, (x)->x]

        local qp = [0.2113;0.7887]
        local qw = [0.5;0.5]
        local Nqp = size(qp,1)

        # Array for storing the pre-evaluated basis functions
        # at the quadrature points
        local phi_hat_qp = Array(Float64,2,Nqp)

        # Pre-evaluate basis functions at (local) QP's
        for q=1:Nqp
            for i=1:2
              phi_hat_qp[i,q] = phi_hat[i](qp[q])
            end
        end

        for k=1:size($mesh.e2t,2)
            if $mesh.e2t[2,k]==0
                # Boundary edge node indices
                local beni1=$mesh.edges[1,k]
                local beni2=$mesh.edges[2,k]
                # Affine mapping
                local A1=$mesh.p[1,beni2]-$mesh.p[1,beni1]
                local A2=$mesh.p[2,beni2]-$mesh.p[2,beni1]
                local b1=$mesh.p[1,beni1]
                local b2=$mesh.p[2,beni1]
                local edgelength=sqrt(A1*A1+A2*A2)
                # Loop over local stiffness matrix
                for i=1:2
                    for q=1:Nqp
                        local x = A1*qp[q]+b1
                        local y = A2*qp[q]+b2
                        local v = phi_hat_qp[i,q]
                        push!(I,$mesh.edges[i,k])
                        push!(V,($bndlin_form)*qw[q]*edgelength)
                    end
                end
            end
        end
        # Finally build the stiffness matrix
        sparsevec(I,V)
    end
end
