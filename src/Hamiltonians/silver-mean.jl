ITensors.space(::SiteType"Trit") = 3

# basis states |0>, |1>, |2>
ITensors.state(::StateName"0", ::SiteType"Trit") = [1, 0, 0]
ITensors.state(::StateName"1", ::SiteType"Trit") = [0, 1, 0]
ITensors.state(::StateName"2", ::SiteType"Trit") = [0, 0, 1]

function sm_MPS(L, A, B)
    sites =  siteinds("Trit", L)
    psi = MPS(sites)
    links = [Index(2, "Link,l=$i") for i in 1:L-1]

    for k in 1:L
        if k == 1
            T = ITensor(sites[k], links[k])
            # b=0: 
            T[sites[k]=>1, links[k]=>1] = 1.0
            # b=1: 
            T[sites[k]=>2, links[k]=>1] = 1.0
            # b=2:
            T[sites[k]=>3, links[k]=>2] = 1.0
            psi[k] = T
        elseif k < L
            T = ITensor(links[k-1], sites[k], links[k])       
            # b=0 (physical index = 1)
            T[links[k-1]=>1, sites[k]=>1, links[k]=>1] = 1.0
            T[links[k-1]=>2, sites[k]=>1, links[k]=>1] = 1.0        
            # b=1 (physical index = 2)
            T[links[k-1]=>1, sites[k]=>2, links[k]=>1] = 1.0
            # b=2 (physical index = 3)
            T[links[k-1]=>1, sites[k]=>3, links[k]=>2] = 1.0
            psi[k] = T
        else
            T = ITensor(links[k-1], sites[k])
            # b=0
            T[links[k-1]=>1, sites[k]=>1] = A
            T[links[k-1]=>2, sites[k]=>1] = A
            # b=1 
            T[links[k-1]=>1, sites[k]=>2] = A
            T[links[k-1]=>2, sites[k]=>2] = 0
            # b=2 
            T[links[k-1]=>1, sites[k]=>3] = B
            T[links[k-1]=>2, sites[k]=>3] = 0
            psi[k] = T
        end
    end
    return psi, sites
end

ITensors.op(::OpName"tau_plus",::SiteType"Trit") =
 [0 0 0
  1 0 0
  0 1 0]

ITensors.op(::OpName"tau_minus",::SiteType"Trit") =
 [0 1 0
  0 0 1
  0 0 0]

#Projection matrix: P01|0> = |0>, P01|1>=|1>, P01|2> = 0
ITensors.op(::OpName"P01",::SiteType"Trit") =
 [1 0 0
  0 1 0
  0 0 0]

#Projection matrix: P0|0> = |0>, P0|1>= P0|2> = 0
ITensors.op(::OpName"P0",::SiteType"Trit") =
 [1 0 0
  0 0 0
  0 0 0]

ITensors.op(::OpName"Id",::SiteType"Trit") =
 [1 0 0
  0 1 0
  0 0 1]

ITensors.op(::OpName"r20",::SiteType"Trit") =
 [0 0 0
  0 0 0
  1 0 0]

ITensors.op(::OpName"r02",::SiteType"Trit") =
 [0 0 1
  0 0 0
  0 0 0]


#Get K matrix in zeckendorf basis.
function kinetic_sm(hopping, L, sites) 
    k1 = OpSum()
    k2 = OpSum()

    #Construct operator k1 (upper diagonal9)
    for i in 1:L
        os = OpSum()
        os += 1,"tau_minus",i

        for j in 1:(i-1)
            os *=  ("Id",j) 
        end

        n = 0
        for k in (i+1):L
            if n == 0 
                os *=  ("r20",k) 
                n = 1
            else
                os *= ("P0", k)
                n = 0
            end
        end        
        k1 += os
    end
    
    #Construct operator k2 (lower diagonal)
    for i in 1:L
        os = OpSum()
        os += 1,"tau_plus",i

        if (i-1) > 0
            os *= ("P01", i-1)
        end

        n = 0
        for k in (i+1):L
            if n == 0 
                os *=  ("r02",k) 
                n = 1
            else
                os *= ("P0", k)
                n = 0
            end
        end     
        k2 += os
    end

    #The opearator K
    k_mpo1 = MPO(k1, sites) 
    k_mpo2 = MPO(k2, sites)
    true_hop_1 = apply(hopping, k_mpo1)
    true_hop_2 = apply(k_mpo2, hopping)
    k_mpo = +(true_hop_1, true_hop_2;  cutoff = 1e-8)

    return k_mpo
end


function sm_Hamiltonian(L, A, B)
    psi, sites = sm_MPS(L, A, B)
    hop = mps_to_diagonal_mpo(psi, sites)
    H = kinetic_sm(hop, L, sites)

    return H, sites
end


