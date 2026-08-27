module TensorQuasicrystals

using ITensors
using ITensorMPS

include("Hamiltonians/common.jl")
include("Hamiltonians/Fibonacci.jl")
include("Hamiltonians/Tribonacci.jl")
include("Hamiltonians/silver-mean.jl")
include("KPM.jl")

end