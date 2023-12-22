
using CSV
using Plots
using DataFrames

include("Grids.jl");
using .Grids;
include("Potentials.jl")
using .Potentials;
include("Hydrogen.jl")
using .Hydrogen;


include("Integral_Numerical_Methods.jl");
using .Integral_Numerical_Methods;


r_max::Float32=10.0;
l::Int32=1;
Z::Int32=1;
w::Float32=1;
m::Float32=1;

expo_grid= Grids.exponential_grid(r_max, Z);
v_angu= Potentials.angular_potential(l, expo_grid);
v_colu= Potentials.coulomb_potential(Z, expo_grid);
#v_harm= Potentials.harmoic_oscilator_potential(w, m, unif_grid);

u_s1_hydr= Hydrogen.u_s1_hydrogen(expo_grid);

E::Float32= -0.5000;
l::Int32=0;
v_effe= v_colu;
init_valu1::Float32=u_s1_hydr[1];#expo_grid[1]^(l+1.0)
init_valu2::Float32=u_s1_hydr[2];#expo_grid[2]^(l+1.0)
f::Vector{Float32}= 2.0.*(v_effe .- E);
g=zeros(Float32, size(f)[1]);

u_pred= Integral_Numerical_Methods.integrate_second_order_DE(expo_grid,g,f,init_valu1,init_valu2)
