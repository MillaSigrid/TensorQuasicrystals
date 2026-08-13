function Tribonacci_MPS(L, A, B, C)
    sites =  siteinds("Qubit",L,conserve_qns=false);
    psi = MPS(sites)
    links = [Index(3, "Link,l=$i") for i in 1:L-1]

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
            T[links[k-1]=>3, sites[k]=>1, links[k]=>1] = 1.0      
            # b=1 (physical index = 2)
            T[links[k-1]=>1, sites[k]=>2, links[k]=>2] = 1.0
            T[links[k-1]=>2, sites[k]=>2, links[k]=>3] = 1.0   
            psi[k] = T
        else
            T = ITensor(links[k-1], sites[k])
            # b=0
            T[links[k-1]=>1, sites[k]=>1] = A
            T[links[k-1]=>2, sites[k]=>1] = A
            T[links[k-1]=>3, sites[k]=>1] = A
            # b=1 => zero vector
            T[links[k-1]=>1, sites[k]=>2] = B
            T[links[k-1]=>2, sites[k]=>2] = C
            psi[k] = T
        end
    end
    return psi, sites
end


function kinetic_trib(hopping, L, sites) 
    
    k1 = OpSum()
    k2 = OpSum()

    #Construct operator k1 (upper diagonal)
    for i in 1:L
        os = OpSum()
        os += 1,"sigma_plus", i

        for j in 1:i-1
            os *=  ("Id",j) 
        end

        n = 1
        for k in i+1:L
            if mod(n, 3) == 0
                os *= ("P0", k)
            else
                os *=  ("sigma_minus",k) 
            end
            n += 1
        end        
        k1 += os
    end

    # Construct operator k2 (lower diagonal)
    # Identify "-00" and do addition
    for i in 1:L
        os = OpSum()
        os += 1,"sigma_minus",i

        if (i-1) > 0
            os *= ("P0", i-1)
        end

        for j in 1:i-2           
            os *=  ("Id",j) 
        end

        n = 1
        for k in i+1:L
            if mod(n, 3) == 0
                os *= ("P0", k)
            else
                os *=  ("sigma_plus",k) 
            end
            n += 1
        end        
        k2 += os
    end

    # Identify "010" and do addition (start form 2 to avoid overlap; ie to have two same operators added to k2)
    for i in 2:L
        os = OpSum()
        os += 1,"sigma_minus",i

        if (i-1) > 0
            os *= ("P1", i-1)
        end

        if (i-2) > 0
            os *= ("P0", i-2)
        end

        for j in 1:i-2           
            os *=  ("Id",j) 
        end

        n = 1
        for k in i+1:L
            if mod(n, 3) == 0
                os *= ("P0", k)
            else
                os *=  ("sigma_plus",k) 
            end
            n += 1
        end        
        k2 += os
    end

    #The opearator K
    k_mpo1 = MPO(k1,sites) 
    true_hop_1 = apply(hopping, k_mpo1)
    k_mpo2 = MPO(k2, sites)
    true_hop_2 = apply(k_mpo2, hopping)
    k_mpo =  +(true_hop_1, true_hop_2;  cutoff = 1e-8)

    return k_mpo
end


function fibonacci_Hamiltonian(L, A, B, C)
    psi, sites = tribonacci_MPS(L, A, B, C)
    hop = mps_to_diagonal_mpo(psi, sites)
    H = kinetic_trib(hop, L, sites)

    return H, sites
end