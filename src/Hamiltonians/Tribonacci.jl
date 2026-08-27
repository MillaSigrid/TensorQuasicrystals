# This file contains two types of methods:
# 1) Construction methods for Tribonacci MPO Hamiltonians
# 2) Methods for obtaining the matrix representations of the Hamiltonians


# Common variables
####################################################


# L       = length of the tensor network (i.e, the number of local tensors)
# N       = physical system size
# A, B, C = Tribonacci word parameters


# 1) Construction methods
#####################################################


# Constructs an MPS representation of the Tribonacci word. Returns the MPS and the sites.
function tribonacci_MPS(L::Int, A::Int, B::Int, C::Int)
    sites =  siteinds("Qubit",L,conserve_qns=false);
    psi = MPS(sites)
    links = [Index(3, "Link,l=$i") for i in 1:L-1]

    for k in 1:L
        if k == 1
            T = ITensor(sites[k], links[k])
            # sigma = 0: 
            T[sites[k]=>1, links[k]=>1] = 1.0
            # sigma = 1: 
            T[sites[k]=>2, links[k]=>2] = 1.0
            psi[k] = T
        elseif k < L
            T = ITensor(links[k-1], sites[k], links[k])      
            # sigma = 0:
            T[links[k-1]=>1, sites[k]=>1, links[k]=>1] = 1.0
            T[links[k-1]=>2, sites[k]=>1, links[k]=>1] = 1.0
            T[links[k-1]=>3, sites[k]=>1, links[k]=>1] = 1.0      
            # sigma = 1:
            T[links[k-1]=>1, sites[k]=>2, links[k]=>2] = 1.0
            T[links[k-1]=>2, sites[k]=>2, links[k]=>3] = 1.0   
            psi[k] = T
        else
            T = ITensor(links[k-1], sites[k])
            # sigma = 0:
            T[links[k-1]=>1, sites[k]=>1] = A
            T[links[k-1]=>2, sites[k]=>1] = A
            T[links[k-1]=>3, sites[k]=>1] = A
            # sigma = 1:
            T[links[k-1]=>1, sites[k]=>2] = B
            T[links[k-1]=>2, sites[k]=>2] = C
            psi[k] = T
        end
    end
    return psi, sites
end

# Encodes the hopping-dependent part of the MPO Hamiltonian. Returns TK+(KT)^*, where T is the input MPO and K is a shift tensor.
function kinetic_trib(mpo::MPO, L::Int, sites) 
    
    k1 = OpSum()
    k2 = OpSum()

    #Construct operator k1 (upper-diagonal part of the matrix representation)
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

    # Construct operator k2 (lower-diagonal part of the matrix representation)
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

    # Identify "010" and do addition (start form 2 to avoid overlap, i.e to have two same operators added to k2)
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

    # Construct final result
    k_mpo1 = MPO(k1,sites) 
    true_hop_1 = apply(mpo, k_mpo1)
    k_mpo2 = MPO(k2, sites)
    true_hop_2 = apply(k_mpo2, mpo)
    k_mpo =  +(true_hop_1, true_hop_2;  cutoff = 1e-8)

    return k_mpo
end

# Constructs a diagonal Tribonacci Hamiltonian. The hopping amplitues are assumed to be equal to a constant, given by the optional argument const_term.
function tribonacci_diag_Hamiltonian(L::Int, A::Int, B::Int, C::Int; const_term = 0)
    psi, sites = tribonacci_MPS(L, A, B, C)
    H = mps_to_diagonal_mpo(psi, sites)

    if const_term != 0
        c_mps = constant_MPS(L, const_term, sites)
        c_mpo = mps_to_diagonal_mpo(c_mps, sites)
        H += kinetic_trib(c_mpo, L, sites)
    end

    return H, sites
end

# Constructs a diagonal Fibonacci Hamiltonian. The hopping amplitues are assumed to be equal to a constant, given by the optional argument const_term.
function tribonacci_off_diag_Hamiltonian(L::Int, A::Int, B::Int, C::Int; const_term = 0)
    psi, sites = tribonacci_MPS(L, A, B, C)
    hop = mps_to_diagonal_mpo(psi, sites)
    H = kinetic_trib(hop, L, sites)

    if const_term != 0
        c_mps = constant_MPS(L, const_term, sites)
        c_mpo = mps_to_diagonal_mpo(c_mps, sites)
        H += c_mpo
    end

    return H, sites
end

# Projector MPO that sets the elements corresponding to invalid Tribonacci strings to zero.
function tribonacci_projection_MPO(L::Int, sites)
    proj_MPS = MPS(sites)
    links = [Index(3, "Link,l=$i") for i in 1:L-1]

    for k in 1:L
        if k == 1
            T = ITensor(sites[k], links[k])
            # sigma = 0: 
            T[sites[k]=>1, links[k]=>1] = 1.0
            # sigma = 1: 
            T[sites[k]=>2, links[k]=>2] = 1.0
            proj_MPS[k] = T
        elseif k < L
            T = ITensor(links[k-1], sites[k], links[k])      
            # sigma = 0:
            T[links[k-1]=>1, sites[k]=>1, links[k]=>1] = 1.0
            T[links[k-1]=>2, sites[k]=>1, links[k]=>1] = 1.0
            T[links[k-1]=>3, sites[k]=>1, links[k]=>1] = 1.0      
            # sigma = 1:
            T[links[k-1]=>1, sites[k]=>2, links[k]=>2] = 1.0
            T[links[k-1]=>2, sites[k]=>2, links[k]=>3] = 1.0   
            proj_MPS[k] = T
        else
            T = ITensor(links[k-1], sites[k])
            # sigma = 0:
            T[links[k-1]=>1, sites[k]=>1] = 1.0
            T[links[k-1]=>2, sites[k]=>1] = 1.0
            T[links[k-1]=>3, sites[k]=>1] = 1.0
            # sigma = 1:
            T[links[k-1]=>1, sites[k]=>2] = 1.0
            T[links[k-1]=>2, sites[k]=>2] = 1.0
            proj_MPS[k] = T
        end
    end

    return mps_to_diagonal_mpo(proj_MPS, sites)
end


# Matrix representation methods
######################################################


# Computes the nth Triboancci number, T_n.
function trib_seq(N::Int)
    T_minus_3 = 0
    T_minus_2 = 0
    T_minus_1 = 1

    if (N == 0) || (N == 1)
        return 0
    elseif N==2
        return 1
    end
    
    res = 0

    for _  in 3:N
        res = T_minus_1+T_minus_2+T_minus_3
        T_minus_3 = T_minus_2
        T_minus_2 = T_minus_1
        T_minus_1 = res
    end

    return res
end

# Computes the Tribonacci representation of n.
function zeck_t(n::Int)
    n <= 0 && return 0
    fib = [4,2,1]; while fib[1] < n insert!(fib,1,sum(fib[1:3])) end
    dig = Int[]; for f in fib f <= n ? (push!(dig,1); n = n-f;) : push!(dig,0) end
    return dig[1] == 0 ? dig[2:end] : dig
end

# Returns the Zeckendorf representation of n in a form ["0", "1", ...], where L is the length of the vector. 
function to_zeck_t_vector(n, size)
    # Convert the integer n to a binary string (without leading zeros)
    zeck_str = join(zeck_t(n))
    
    # Pad the binary string with leading zeros on the left
    padded = lpad(zeck_str, size, '0')
    
    # Convert the padded string into a vector of characters,
    # then map each character to a proper String ("0" or "1")
    return collect(padded) |> x -> map(s -> string(s), x)
end

# Returns the Tribonacci MPS representation of integer n.
function zeck_t_to_MPS(n::Integer, L::Integer, sites)
    return MPS(sites, to_zeck_t_vector(n, L))
end

# Returns the matrix representation of a MPO Hamiltonian. Use only for small system sizes.
function get_matrix_trib(mpo::MPO, L::Int, sites)
    size = trib_seq(L+3)
    mat = zeros(size, size) #.+ 0im
 
    # Loop over all possible bra (row) basis states
    for i in 0:size-1
        # Loop over all possible ket (column) basis states
        for j in 0:size-1
            element = inner(
                MPS(sites, to_zeck_t_vector(Int(i), Int(L)))', 
                mpo, 
                MPS(sites, to_zeck_t_vector(Int(j), Int(L)))
            )
            if abs(element) <= 1e-10
                element = 0 #+0im
            end

            mat[i+1, j+1] = real(element)
        end
    end
    
    return mat
end