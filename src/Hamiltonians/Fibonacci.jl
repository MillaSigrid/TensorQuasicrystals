# This file contains two types of methods:
# 1) Construction methods for Fibonacci MPO Hamiltonians
# 2) Methods for obtaining the matrix representations of the Hamiltonians


# Common variables
####################################################


# L    = length of the tensor network (i.e, the number of local tensors)
# N    = physical system size
# A, B = Fibonacci word parameters


# 1) Construction methods
#####################################################


# Constructs an MPS representation of the Fibonacci word. Returns the MPS and the sites.
function fibonacci_MPS(L, A, B)
    sites = siteinds("Qubit", L,conserve_qns=false)
    psi  =  MPS(sites)
    links = [Index(2, "Link,l=$i") for i in 1:L-1]

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
            # sigma = 1:
            T[links[k-1]=>1, sites[k]=>2, links[k]=>2] = 1.0 
            psi[k] = T
        else
            T = ITensor(links[k-1], sites[k])
            # sigma = 0 
            T[links[k-1]=>1, sites[k]=>1] = A
            T[links[k-1]=>2, sites[k]=>1] = A
            # sigma = 1
            T[links[k-1]=>1, sites[k]=>2] = B
            psi[k] = T
        end
    end
    
    return psi, sites
end

# Encodes the hopping-dependent part of the MPO Hamiltonian. Returns TK+(KT)^*, where T is the input MPO and K is a shift tensor.
function kinetic_fib(mpo, L, sites; boundary=:OBC) 

    if boundary != :OBC && boundary != :PBC
        error("boundary must be :OBC or :PBC")
    end

    k1 = OpSum()
    k2 = OpSum()

    # Construct operator k1 (upper-diagonal part of the matrix representation)
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

    # Construct operator k2 (lower diagonal part of the matrix representation)
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
        # Add PBC part to k1 
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

        # Add PBC part to k2
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
    
    # Construct final result
    k_mpo1 = MPO(k1,sites) 
    true_hop_1 = apply(mpo, k_mpo1)
       
    k_mpo2 = MPO(k2, sites)
    true_hop_2 = apply(k_mpo2, mpo)

    k_mpo =  +(true_hop_1, true_hop_2;  cutoff = 1e-8)
    
    return k_mpo
end

# Constructs a diagonal Fibonacci Hamiltonian. Supports open (:OBC) and periodic (:PBC) boundary conditions.
# The hopping amplitues are assumed to be equal to a constant, given by the optional argument const_term.
function fibonacci_diag_Hamiltonian(L, A, B; boundary=:OBC, const_term = 0)
    psi, sites = fibonacci_MPS(L, A, B)
    H = mps_to_diagonal_mpo(psi, sites)

    if const_term != 0
        c_mps = constant_MPS(L, const_term, sites)
        c_mpo = mps_to_diagonal_mpo(c_mps, sites)
        H += kinetic_fib(c_mpo, L, sites; boundary=boundary)
    end

    return H, sites
end

# Constructs an off-diagonal Fibonacci Hamiltonian. Supports open (:OBC) and periodic (:PBC) boundary conditions.
# The on-site potentials are assumed to be equal to a constant, given by the optional argument const_term.
function fibonacci_off_diag_Hamiltonian(L, A, B; boundary=:OBC, const_term = 0)
    psi, sites = fibonacci_MPS(L, A, B)
    hop = mps_to_diagonal_mpo(psi, sites)
    H = kinetic_fib(hop, L, sites; boundary=boundary)

    if const_term != 0
        c_mps = constant_MPS(L, const_term, sites)
        c_mpo = mps_to_diagonal_mpo(c_mps, sites)
        H += c_mpo
    end

    return H, sites
end


# Matrix representation methods
######################################################


# Computes the nth Fiboancci number, F_n.
function fib_seq(n::Int)
    φ = (1 + sqrt(5))/2
    ψ = (1 - sqrt(5))/2
    fib_n = 1/sqrt(5) * (φ^n - ψ^n)
    return Int(round(fib_n))
end

# Computes the Fibonacci (Zeckendorf) representation of n.
function zeck(n)
    n <= 0 && return 0
    fib = [2,1]; while fib[1] < n insert!(fib,1,sum(fib[1:2])) end
    dig = Int[]; for f in fib f <= n ? (push!(dig,1); n = n-f;) : push!(dig,0) end
    return dig[1] == 0 ? dig[2:end] : dig
end

# Returns the Zeckendorf representation of n in a form ["0", "1", ...], where L is the length of the vector. 
function to_zeck_vector(n::Int, L::Int)
    # Convert the integer n to a zeckendorf string
    zeck_str = join(zeck(n))
    
    # Pad the Fibonacci string with leading zeros on the left
    padded = lpad(zeck_str, L, '0')
    
    # Convert the padded string into a vector of characters,
    # then map each character to a proper String ("0" or "1")
    return collect(padded) |> x -> map(s -> string(s), x)
end

# Returns the MPS Fibonacci representation of integer n.
function zeck_to_MPS(n::Int, L::Int, sites)
    return MPS(sites, to_zeck_vector(n, L))
end

# Returns the matrix representation of a MPO Hamiltonian. Use only for small system sizes.
function get_matrix_fib(mpo, L, sites)
    size = fib_seq(L+2)
    mat = zeros(size, size) #+0im
   
    # Loop over all possible bra (row) basis states
    for i in 0:size-1
        # Loop over all possible ket (column) basis states
        for j in 0:size-1
            element = inner(
                zeck_to_MPS(i, L, sites)', 
                mpo, 
                zeck_to_MPS(j, L, sites)
            )
            if abs(element) <= 1e-10
                element = 0 #+0im
            end
            mat[i+1, j+1] = real(element)
        end
    end
    return mat
end