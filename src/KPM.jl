
# Computes Fibonacci trace.
function trace_fib(mpo::MPO, sites)
    N = length(mpo)
    T1 = mpo[1]
    s = sites[1]

    res0 = T1*state(s,1)*state(s',1) # strings ending in 0
    res1 = T1*state(s,2)*state(s',2) # strings ending in 1

    for i = 2:N
        Ti = mpo[i]
        s = sites[i]

        new_res0 = (res0+res1)*Ti*state(s,1)*state(s',1)
        new_res1 = res0*Ti*state(s,2)*state(s',2)

        res0 = new_res0
        res1 = new_res1

    end

    return scalar(res0)+scalar(res1)
end

# Computes Tribonacci trace.
function trace_trib(mpo::MPO, sites)
    N = length(mpo)
    T1 = mpo[1]
    s = sites[1]

    resx0 = T1*state(s,1)*state(s',1)  # strings ending in 0
    res01 = T1*state(s,2)*state(s',2)  # strings ending in 01
    res11 = 0*res01                    # strings ending in 11

    for i = 2:N
        Ti = mpo[i]
        s = sites[i]

        new_resx0 = (resx0+res01+res11)*Ti*state(s,1)*state(s',1)
        new_res01 = resx0*Ti*state(s,2)*state(s',2)
        new_res11 = res01*Ti*state(s,2)*state(s',2)

        res11 = new_res11
        resx0 = new_resx0
        res01 = new_res01
    end

    return scalar(resx0)+scalar(res01)+scalar(res11)
end

# Computes silver-mean trace.
function trace_sm(mpo::MPO, sites)
    N = length(mpo)
    T1 = mpo[1]
    s = sites[1]

    res0 = T1*state(s,1)*state(s',1)  # strings ending in 0
    res1 = T1*state(s,2)*state(s',2)  # strings ending in 1
    res2 = T1*state(s,3)*state(s',3)  # strings ending in 11

    for i = 2:N
        Ti = mpo[i]
        s = sites[i]

        new_res0 = (res0+res1+res2)*Ti*state(s,1)*state(s',1)
        new_res1 = (res0+res1)*Ti*state(s,2)*state(s',2)
        new_res2 = (res0+res1)*Ti*state(s,3)*state(s',3)

        res0 = new_res0
        res1 = new_res1
        res2 = new_res2
    end

    return scalar(res0)+scalar(res1)+scalar(res2)
end

# Returns a trace function corresponding to the used model (:fib, :trib or :sm).
function get_trace_function(model::Symbol)
    if model == :fib
        return trace_fib
    elseif model == :trib
        return trace_trib
    elseif model == :sm
        return trace_sm
    else
        error("Unknown model: $model")
    end
end

# Returns the system size given the used model and the number of local tensors L.
function get_physical_system_size(L::Int, model::Symbol)
    if model == :fib
        return fib_seq(L+2)
    elseif model == :trib
        return trib_seq(L+3)
    elseif model == :sm
        return sm_seq(L)
    else
        error("Unknown model: $model")
    end
end

# Kernel Polynomial Method for computing the Tn polynomials of the Hamiltonian H. Nmu is the number of Tn polynomials to compute
# E is the energy grid (physical units). Model can be :fib, :sm or :trib.
function KPM_mu(H::MPO, sites, model::Symbol, Nmu::Int, E; cutoff=10^-6, maxdim=200)
    apply_kwargs = (cutoff=cutoff, maxdim=maxdim)
    #sites = getindex.(siteinds(H), 2)
    #sites = [siteind(H, i) for i in 1:length(H)]   # ITensorMPS convenience: grabs the plev=0 site index

    # Select trace function
    trace_function = get_trace_function(model)

    Id_op = MPO(sites, "Id")

    # Normalize H using the provided energy window E
    e  = 1e-2                              # small epsilon: keeps spec away from ±1
    W2 = (maximum(E) - minimum(E)) / 2     # half-width of the energy range
    a  = (maximum(E) + minimum(E)) / 2     # center of the energy range
    # Shift & rescale so Ĥ ≈ (H - a I)/(W2 + e) has spectrum in (-1, 1)
    Ham_n = (H - a*Id_op) / (W2 + e)

    Tn_minus_2 = Id_op
    Tn_minus_1 = Ham_n

    mu = Vector{Float64}(undef, Nmu)

    mu[1] = trace_function(Tn_minus_2, sites)
    mu[2] = trace_function(Tn_minus_1, sites)


    # Recurrence: T_k = 2 Ĥ T_{k-1} - T_{k-2}
    for k in 3:Nmu
        Tn_k = +(2 * apply(Ham_n, Tn_minus_1; apply_kwargs...),
               -Tn_minus_2; apply_kwargs...)
        Tn_k = ITensorMPS.truncate!(Tn_k; apply_kwargs...)

	    if ((k+200)%200)==0
		    println("k=$(k), maxlinkdim=$(maxlinkdim(Tn_k))")
	    end

        Tn_minus_2 = Tn_minus_1
        Tn_minus_1 = Tn_k

        mu[k] = trace_function(Tn_k, sites)
    end

    return mu
end

# Computes the DOS given the moments and rescaled energies.
function get_DOS_from_moments(mu, Nmu::Int, E)
    jackson_kernel = [(Nmu - n) * cos(π * n / Nmu) + sin(π * n / Nmu) / tan(π / Nmu) for n in 0:Nmu-1] / Nmu

    A = mu[1] * cos(0 * acos(E)) * jackson_kernel[1]
    for n in 2:Nmu
        A += 2 * mu[n] * cos((n-1) * acos(E)) * jackson_kernel[n]
    end

    return A / (π * sqrt(1 - E^2)) 
end

# Computes the DOS of the input Hamiltonian by using KMP with Nmu moments. The used model can be :fib, :sm or :trib.
function get_DOS(H::MPO, sites, model::Symbol, Nmu::Int, E)
    L = length(H)
    N = get_physical_system_size(L, model)
    mu_list = KPM_mu(H, sites, model, Nmu, E)/N
    E_scaled = range(-1,1, length(E))
    DOS = [get_DOS_from_moments(mu_list, Nmu, E) for E in E_scaled]

    return DOS
end