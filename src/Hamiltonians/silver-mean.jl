# This file contains two types of methods:
# 1) Construction methods for silver-mean MPO Hamiltonians
# 2) Methods for obtaining the matrix representations of the Hamiltonians


# Common variables
####################################################


# L    = length of the tensor network (i.e, the number of local tensors)
# N    = physical system size
# A, B = silver-mean word parameters


# Custom SiteType
####################################################


# define custom SiteType
ITensors.space(::SiteType"Trit") = 3

# basis states |0>, |1>, |2>
ITensors.state(::StateName"0", ::SiteType"Trit") = [1, 0, 0]
ITensors.state(::StateName"1", ::SiteType"Trit") = [0, 1, 0]
ITensors.state(::StateName"2", ::SiteType"Trit") = [0, 0, 1]

# local operators
ITensors.op(::OpName"tau_plus",::SiteType"Trit") =
 [0 0 0
  1 0 0
  0 1 0]

ITensors.op(::OpName"tau_minus",::SiteType"Trit") =
 [0 1 0
  0 0 1
  0 0 0]

#Projection opearator: P01|0> = |0>, P01|1>=|1>, P01|2> = 0
ITensors.op(::OpName"P01",::SiteType"Trit") =
 [1 0 0
  0 1 0
  0 0 0]

#Projection opearator: P0|0> = |0>, P0|1>= P0|2> = 0
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


# 1) Construction methods
#####################################################


# Constructs an MPS representation of the silver-mean word. Returns the MPS and the sites.
function sm_MPS(L::Int, A::Int, B::Int)
    sites =  siteinds("Trit", L)
    psi = MPS(sites)
    links = [Index(2, "Link,l=$i") for i in 1:L-1]

    for k in 1:L
        if k == 1
            T = ITensor(sites[k], links[k])
            # sigma = 0: 
            T[sites[k]=>1, links[k]=>1] = 1.0
            # sigma = 1: 
            T[sites[k]=>2, links[k]=>1] = 1.0
            # sigma = 2:
            T[sites[k]=>3, links[k]=>2] = 1.0
            psi[k] = T
        elseif k < L
            T = ITensor(links[k-1], sites[k], links[k])       
            # sigma = 0:
            T[links[k-1]=>1, sites[k]=>1, links[k]=>1] = 1.0
            T[links[k-1]=>2, sites[k]=>1, links[k]=>1] = 1.0        
            # sigma = 1:
            T[links[k-1]=>1, sites[k]=>2, links[k]=>1] = 1.0
            # sigma = 2:
            T[links[k-1]=>1, sites[k]=>3, links[k]=>2] = 1.0
            psi[k] = T
        else
            T = ITensor(links[k-1], sites[k])
            # sigma = 0:
            T[links[k-1]=>1, sites[k]=>1] = A
            T[links[k-1]=>2, sites[k]=>1] = A
            # sigma = 1:
            T[links[k-1]=>1, sites[k]=>2] = A
            T[links[k-1]=>2, sites[k]=>2] = 0
            # sigma = 2:
            T[links[k-1]=>1, sites[k]=>3] = B
            T[links[k-1]=>2, sites[k]=>3] = 0
            psi[k] = T
        end
    end
    return psi, sites
end

# Encodes the hopping-dependent part of the MPO Hamiltonian. Returns TK+(KT)^*, where T is the input MPO and K is a shift tensor.
function kinetic_sm(mpo::MPO, L::Int, sites) 
    k1 = OpSum()
    k2 = OpSum()

    # Construct operator k1 (upper-diagonal part of the matrix representation)
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
    
    # Construct operator k2 (lower-diagonal part of the matrix representation)
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

    # Construct final result
    k_mpo1 = MPO(k1, sites) 
    k_mpo2 = MPO(k2, sites)
    true_hop_1 = apply(mpo, k_mpo1)
    true_hop_2 = apply(k_mpo2, mpo)
    k_mpo = +(true_hop_1, true_hop_2;  cutoff = 1e-8)

    return k_mpo
end

function constant_MPS_trit(L::Int, const_term::InterruptException, sites)
    c_mps = MPS(sites)

    for k in 1:L
        T = ITensor(sites[k])

        T[sites[k]=>1] = const_term^(1/L)
        T[sites[k]=>2] = const_term^(1/L)
        T[sites[k]=>3] = const_term^(1/L)

        c_mps[k] = T
    end

    return c_mps
end

# Constructs a diagonal silver-mean Hamiltonian. The hopping amplitues are assumed to be equal to a constant, given by the optional argument const_term.
function sm_diag_Hamiltonian(L::Int, A::Int, B::Int; const_term = 0)
    psi, sites = sm_MPS(L, A, B)
    H = mps_to_diagonal_mpo(psi, sites)

    if const_term != 0
        c_mps = constant_MPS_trit(L, const_term, sites)
        c_mpo = mps_to_diagonal_mpo(c_mps, sites)
        H += kinetic_sm(c_mpo, L, sites)
    end

    return H, sites
end

# Constructs an off-diagonal Fibonacci Hamiltonian. The on-site potentials are assumed to be equal to a constant, given by the optional argument const_term.
function sm_off_diag_Hamiltonian(L::Int, A::Int, B::Int; const_term = 0)
    psi, sites = sm_MPS(L, A, B)
    hop = mps_to_diagonal_mpo(psi, sites)
    H = kinetic_sm(hop, L, sites)

    if const_term != 0
        c_mps = constant_MPS_trit(L, const_term, sites)
        c_mpo = mps_to_diagonal_mpo(c_mps, sites)
        H += c_mpo
    end

    return H, sites
end

# Projector MPO that sets the elements corresponding to invalid silver-mean strings to zero.
function sm_projection_MPO(L::Int, sites)
    proj_MPS = MPS(sites)
    links = [Index(2, "Link,l=$i") for i in 1:L-1]

    for k in 1:L
        if k == 1
            T = ITensor(sites[k], links[k])
            # sigma = 0: 
            T[sites[k]=>1, links[k]=>1] = 1.0
            # sigma = 1: 
            T[sites[k]=>2, links[k]=>1] = 1.0
            # sigma = 2:
            T[sites[k]=>3, links[k]=>2] = 1.0
            proj_MPS[k] = T
        elseif k < L
            T = ITensor(links[k-1], sites[k], links[k])       
            # sigma = 0:
            T[links[k-1]=>1, sites[k]=>1, links[k]=>1] = 1.0
            T[links[k-1]=>2, sites[k]=>1, links[k]=>1] = 1.0        
            # sigma = 1:
            T[links[k-1]=>1, sites[k]=>2, links[k]=>1] = 1.0
            # sigma = 2:
            T[links[k-1]=>1, sites[k]=>3, links[k]=>2] = 1.0
            proj_MPS[k] = T
        else
            T = ITensor(links[k-1], sites[k])
            # sigma = 0:
            T[links[k-1]=>1, sites[k]=>1] = 1.0
            T[links[k-1]=>2, sites[k]=>1] = 1.0
            # sigma = 1:
            T[links[k-1]=>1, sites[k]=>2] = 1.0
            T[links[k-1]=>2, sites[k]=>2] = 0
            # sigma = 2:
            T[links[k-1]=>1, sites[k]=>3] = 1.0
            T[links[k-1]=>2, sites[k]=>3] = 0
            proj_MPS[k] = T
        end
    end

    return mps_to_diagonal_mpo(proj_MPS, sites)
end


# Matrix representation methods
######################################################

# Computes the nth denominator of the continued-fraction convergents of the silver mean, q_n.
function sm_seq(N::Int)

    if N == 0 return 1 end
    if N == 1 return 3 end

    q_minus_2 = 1
    q_minus_1 = 3
    res = 0

    for i = 2:N
        res = 2*q_minus_1+q_minus_2
        q_minus_2 = q_minus_1
        q_minus_1 = res
    end

    return res
end

# Computes the silver-mean representation of n.
function ostrowski_sm(n::Int)
    n <= 0 && return 0
    qs = [3,1]; while qs[1] < n insert!(qs,1,2*qs[1]+qs[2]) end
    dig = Int[]; for f in qs d = div(n,f); push!(dig,d); n -= d*f end
    return dig[1] == 0 ? dig[2:end] : dig
end

# Returns the silver-mean representation of n in a form [0, 1, ...], where L is the length of the vector. 
function to_sm_vector(n::Int, size::Int)
    # Convert the integer n to a ternary string (without leading zeros)
    sm_str = join(ostrowski_sm(n))
    
    # Pad the ternary string with leading zeros on the left
    padded = lpad(sm_str, size, '0')
    
    # Convert characters '0', '1', '2' to integers 0, 1, 2
    return [parse(Int, c) for c in padded]
end

# Returns the MPS silver-mean representation of integer n.
function sm_to_MPS(n::Int, L::Int, sites)
    sm = to_sm_vector(n, L)
    return MPS([state(sites[j], sm[j] + 1) for j in eachindex(sm)])
end

# Returns the matrix representation of a silver-mean MPO Hamiltonian. Use only for small system sizes.
function get_matrix_sm(mpo::MPO, L::Int, sites)
    size = sm_seq(L)
    mat = zeros(size, size) #.+ 0im
 
    # Loop over all possible bra (row) basis states 
    for i in 0:size-1
        # Loop over all possible ket (column) basis states
        for j in 0:size-1
            element = inner(
                sm_to_MPS(i, L, sites)', 
                mpo, 
                sm_to_MPS(j, L, sites)
            )
            if abs(element) <= 1e-10
                element = 0 #+0im
            end
            mat[i+1, j+1] = real(element)
        end
    end
    
    return mat
end