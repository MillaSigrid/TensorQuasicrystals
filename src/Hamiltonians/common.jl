#This file contains common methods that are used in files Fibonacci.jl, Tribonacci.jl and silver-mean.jl.

# Local operators used in Fibonacci.jl and Tribonacci.jl
ITensors.op(::OpName"P1",::SiteType"Qubit") =
[0 0
0 1]

ITensors.op(::OpName"sigma_plus",::SiteType"Qubit") =
 [0 1
  0 0]

ITensors.op(::OpName"sigma_minus",::SiteType"Qubit") =
 [0 0
  1 0]

#Projection matrix: P0|up> = |up>, P0|down> = 0
ITensors.op(::OpName"P0",::SiteType"Qubit") =
 [1 0
  0 0]

# Promotes an input MPS to a diagonal MPO
function mps_to_diagonal_mpo(mps,sites)
    N = length(mps) 
    mpo_tensors = Vector{ITensor}(undef, N)
    for i in 1:N
        mps_t = mps[i]
        local old_s
        if i == 1
            old_s = uniqueind(mps_t, mps[i+1])
        elseif i == N
            old_s = uniqueind(mps_t, mps[i-1])
        else
            old_s = uniqueind(mps_t, mps[i-1], mps[i+1])
        end
        s = sites[i]      
        s_p = s'          
        s_temp = Index(dim(s), "temp")
        mpo_tensors[i] = replaceind(mps_t, old_s => s_temp) * delta(s_temp, s, s_p)
    end    
    return MPO(mpo_tensors)
end

# Constructs an MPS where every element is equal to const_term.
function constant_MPS(L, const_term, sites)
    c_mps = MPS(sites)

    for k in 1:L
        T = ITensor(sites[k])

        T[sites[k]=>1] = const_term^(1/L)
        T[sites[k]=>2] = const_term^(1/L)

        c_mps[k] = T
    end

    return c_mps
end