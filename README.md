# Quasicrystals-as-tensor-networks

This repository contains methods for representing quasicrystals as tensor networks and accompanies the manuscript "One-dimensional quasicrystals with tensor-network finite-state automata" by Milla Kolehmainen, Jose L. Lado and Anouar Moustaj (https://arxiv.org/pdf/2609.06040). 

The repository includes:
* all construction methods for quasicrystal matrix product operator (MPO) Hamiltonians
* an example workflow for using the constructions to compute the density of states (DOS)

## Repository structure

### `TensorQuasicrystals.jl/src`

**`Hamiltonians`**  
* tools for constructing Fibonacci, Tribonacci and silver-mean Hamiltonians as MPOs
* choose between off-diagonal and diagonal model and specify the sequence parameters to use
* all models are implemented with open boundary conditions by default but the Fibonacci Hamiltonian has an option for periodic boundary conditions as well

**`KPM.jl`**  
* tools for computing the DOS using tensorized kernel polynomial method

### `TensorQuasicrystals.jl/examples/DOS_example.ipynb`
* an example of how to use the methods in this repository to compute the DOS

## Dependencies

The code is written in Julia and depends on the following packages:
* ITensorMPS.jl
* ITensors.jl
* Plots.jl (only used in examples)

## Acknowledgement

The code is built upon the package [TensorBinding.jl](https://tensorbinding.github.io/TensorBinding.jl/) by Anouar Moustaj, Tiago Antão and Yitao Sun. 

