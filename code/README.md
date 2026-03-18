# Numerical experiments
The test cases have been run with Julia v1.10.6. It is strongly recommended to start julia with multiple threads, as some simulations take hours to be completed. To do that, start julia with the following command:
```bash
julia --project=. --threads=num
```
where `num` is the number of threads that you want to use for that session. Alternatively, you can also set `--threads=auto` to use a reasonable default number of threads for your system.

The equations are defined inside the folder `equations/`,  as well as other utility functions.

Before running the simulations, please activate the project and instantiate, i.e., start Julia as described above and run the following code in the Julia REPL:
```julia
julia> using Pkg
julia> Pkg.activate(".")
julia> Pkg.instantiate()
```

## Test case 1: Variable-coefficient advection equation

Run the following command to reproduce the results
```julia
julia> include("numerical_1.jl")
```

## Test case 2: General polynomial equation
Run the following commands to reproduce the results
```julia
julia> include("numerical_2.jl")
```

## Test case 3: Hyperbolized Sainte-Marie equations
Run the following command to reproduce the results about the entropy conservation errors in 1D
```julia
julia> include("escalante_sainte_marie.jl")
```
and
```julia
julia> include("escalante_sainte_marie_well_balanced_2d.jl")
```
to reproduce the results about the well-balancedness on curved mesh.

## Test case 4: Compressible Euler equations in nonconservative form
Run the following command to reproduce the results about convergence analysis with SBP operators
```julia
julia> include("sbp_eoc.jl")
```
and
```julia
julia> include("baroclinic.jl")
```
to reproduce the baroclinic instability test case.
Run the following command to create the plot shown in the paper
```julia
julia> include("plot_baroclinic.jl")
```