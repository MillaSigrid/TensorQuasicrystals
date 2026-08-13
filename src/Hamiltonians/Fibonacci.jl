#Local operators

function fibonacci_MPS(L, A, B)
    sites = siteinds("Qubit",L,conserve_qns=false)
    psi  =  MPS(sites)
    links = [Index(2, "Link,l=$i") for i in 1:L-1]

    for k in 1:L
        if k == 1
            T = ITensor(sites[k], links[k])
            # b=0: 
            T[sites[k]=>1, links[k]=>1] = 1.0
            # b=1: 
            T[sites[k]=>2, links[k]=>2] = 1.0
            psi[k] = T
        elseif k < L
            T = ITensor(links[k-1], sites[k], links[k])
            # b=0 (physical index = 1)
            T[links[k-1]=>1, sites[k]=>1, links[k]=>1] = 1.0
            T[links[k-1]=>2, sites[k]=>1, links[k]=>1] = 1.0        
            # b=1 (physical index = 2)
            T[links[k-1]=>1, sites[k]=>2, links[k]=>2] = 1.0 
            psi[k] = T
        else
            T = ITensor(links[k-1], sites[k])
            # b=0 
            T[links[k-1]=>1, sites[k]=>1] = A
            T[links[k-1]=>2, sites[k]=>1] = A
            # b=1
            T[links[k-1]=>1, sites[k]=>2] = B
            psi[k] = T
        end
    end
    
    return psi, sites
end


function kinetic_fib(hopping, L, sites; boundary=:OBC) 

    if boundary != :OBC && boundary != :PBC
        error("boundary must be :OBC or :PBC")
    end

    k1 = OpSum()
    k2 = OpSum()

    #Construct operator k1 (upper diagonal)
    for i in 1:L
        os = OpSum()
        os += 1,"sigma_plus",L-(i-1)

        for j in 1:L-i 
            os *=  ("Id",j) 
        end

        n = 0
        for k in L+2-i :L
            if n == 0 
                os *=  ("sigma_minus",k) 
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
        os += 1,"sigma_minus",L-(i-1)

        if i < L
            os *= ("P0", L-i)
        end

        for j in 2:L-i           
            os *=  ("Id",j) 
        end

        n = 0
        for k in L+2-i :L
            if n == 0 
                os *=  ("sigma_plus",k) 
                n = 1
            else
                os *= ("P0", k)
                n = 0
            end
        end        
        k2 += os
    end

    if boundary == :PBC
        #Add PBC part to k1
        os_pbc1 = OpSum()
        os_pbc1 += 1, "sigma_minus", 1

        n = 0
        for k in 2:L
            if n == 0
                os_pbc1 *= ("P0", k)
                n = 1
            else
                os_pbc1 *= ("sigma_minus", k)
                n = 0
            end
        end
        k1 += os_pbc1

        #Add PBC part to k2
        os_pbc2 = OpSum()
        os_pbc2 += 1, "sigma_plus", 1

        n = 0
        for k in 2:L
            if n == 0
                os_pbc2 *= ("P0", k)
                n = 1
            else
                os_pbc2 *= ("sigma_plus", k)
                n = 0
            end
        end
        k2 += os_pbc2
    end
    
    #Construct opearator K
    k_mpo1 = MPO(k1,sites) 
    true_hop_1 = apply(hopping, k_mpo1)
       
    k_mpo2 = MPO(k2, sites)
    true_hop_2 = apply(k_mpo2, hopping)

    k_mpo =  +(true_hop_1, true_hop_2;  cutoff = 1e-8)
    
    return k_mpo
end

function fibonacci_diag_Hamiltonian(L, A, B; boundary=:OBC, const_term = 0)
    psi, sites = fibonacci_MPS(L, A, B)
    H = mps_to_diagonal_mpo(psi, sites)

    if const_term != 0
        c_mps = const_MPS(L, const_term, sites)
        H += c_mps
    end

    return H, sites
end

function fibonacci_off_diag_Hamiltonian(L, A, B; boundary=:OBC, const_term = 0)
    psi, sites = fibonacci_MPS(L, A, B)
    hop = mps_to_diagonal_mpo(psi, sites)
    H = kinetic_fib(hop, L, sites; boundary=boundary)

    if const_term != 0
        pot = constant_MPS(L, const_term, sites)
        H += pot
    end

    return H, sites
end