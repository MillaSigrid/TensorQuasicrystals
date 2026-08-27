# Quasicrystals-as-tensor-networks

Methods for representing quasicrystals as tensor networks

## Repository structure

### `TensorQuasicrystals.jl/src`

**`Hamiltonians`**  
* tools for constructing Fibonacci, Tribonacci and silver-mean tight-binding Hamiltonians as matrix product operators
* choose between off-diagonal and diagonal model and specify used sequence parameters
* all models are implemented with open boundary conditions by default but Fibonacci Hamiltonian has an option for periodic boundary conditions as well

**`KPM.jl`**  
* tools for computing the DOS using tenderized kernel polynomial method

### `TensorQuasicrystals.jl/examples/DOS_example.jl`
* an example of how to use the methods in this repository to compute the DOS

## Dependencies

The code is written in Julia and depends on the following packages:
* ITensorMPS.jl
* ITensors.jl
* Plots.jl (only used in examples)

## Acknowledgement

