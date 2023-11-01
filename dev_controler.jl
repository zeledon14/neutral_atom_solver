include("Grid.jl")
using .Grid


r_min=-2.0
r_max=2.0
N=10
unif_grid= Grid.uniform_grid(r_min, r_max, N)