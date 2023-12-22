using CSV
using Plots
using DataFrames

include("Grid.jl")
using .Grid

include("Density.jl")
using .Density

include("Process_Quantum_Numbers.jl")
using .Process_Quantum_Numbers: Process_Quantum_Numbers as PQN

include("Atomic_Schrodinger_Equation.jl")
using .Atomic_Schrodinger_Equation: Atomic_Schrodinger_Equation as  ASE

include("Density_Dependent_Atomic_potentials.jl")
using .Density_Dependent_Atomic_potentials: Density_Dependent_Atomic_potentials as DDAP
#include("Atomic_Poisson_Equation.jl")
#using .Atomic_Poisson_Equation: Atomic_Poisson_Equation as APE

include("Hydrogen.jl")
using .Hydrogen

include("Utils.jl")
using .Utils


config=Dict([
    #space grid parameters
    ("r_max",50.0),
    #elements parameters
    ("Z",18),
    ("Quantum_Numbers_file_path","/home/arturo_hernandez/Desktop/physics_coding_projects/Data/Quantum_Numbers.csv")
]);
quan_numb_df=DataFrame(CSV.File(config["Quantum_Numbers_file_path"]));

grid= Grid.NIST_exponential_grid(config["r_max"], config["Z"])
#delta=0.001
#grid=[delta*i for i in 1:10000]
grid_sqrt= grid.^2.0
grid_bwrd= reverse(grid)
occupation_by_shell, n_by_shell, l_by_shell= Process_Quantum_Numbers.get_shells_occupations_and_quantum_numbers(config["Z"], quan_numb_df);

Z=config["Z"]
eigen_guesses= [-0.5*(Float32(config["Z"])^2.0)/((n)^2.0) for n in n_by_shell];

N=size(grid)[1]
N_search=50
step=0
density_in= zeros(Float32, N)
v_h= zeros(Float32, N)
v_x= zeros(Float32, N)
v_c= zeros(Float32, N)
E_total=1.0
E_total_before=2.0;

while abs(E_total - E_total_before) > 10.0e-8

    #if step == 0
    #    density_in= zeros(Float32, N)
    #    v_h= zeros(Float32, N)
    #    v_x= zeros(Float32, N)
    #    v_c= zeros(Float32, N)
    #end
    eigen_interval_by_shell=[]
    Ea=0.0
    Ep=0.0
    for (i_shell,l) in enumerate(l_by_shell)
        i_guess= eigen_guesses[i_shell]
        delta= 0.05*abs(i_guess)
        start= i_guess - delta
        Ea=float(start)
        _, _, _, ua, _= ASE.integrate_SE(grid,grid_bwrd,v_h,v_x,v_c,l,Z,Ea);
        for n in 1:N_search
            Ep=start + n*delta
            _, _, _, up, _= ASE.integrate_SE(grid,grid_bwrd,v_h,v_x,v_c,l,Z,Ep);
            if Integer(sign(ua)) != Integer(sign(up))
                append!(eigen_interval_by_shell, (Ea,Ep))
                break
            end
            Ea= float(Ep)
            ua= float(up)
        end
    end

    temp_u_basis_set=Matrix{Dict{String, Any}}(undef, length(l_by_shell),5)
    for (i_shell,l) in enumerate(l_by_shell)
        nodes=eigen_interval_by_shell[(2*i_shell -1):(2*i_shell)]
        u, ub, uf, node, E= ASE.illinois_eigen_finder(nodes,grid,grid_bwrd,v_h,v_x,v_c,l,Z);
        #println("l ", l)
        #println("E ", E)
        #println("node ", node)
        #display(plot(grid[node:node+1000],ub[node:node+1000]))
        #display(plot(grid[1:node+100],uf[1:node+100]))
        #display(plot(grid[10:7500],u[10:7500]))
        eigen_guesses[i_shell]=float(E)
        #display(plot(grid,u))
        temp_u_basis_set[i_shell]= Dict([("E",E),("Z",Z),("u",u),("l",l),("occu",occupation_by_shell[i_shell])]);
    end

    density_out= Density.calculate_density(grid, temp_u_basis_set);
    global density_in= Density.linear_mixing(density_in, density_out, alpha=0.25);
    v_hp,  v_xp, E_xp, v_cp,  E_cp= DDAP.density_potentials(grid, density_in, Z);
    global v_h = v_hp
    global v_x = v_xp
    global E_x = E_xp
    global v_c = v_cp
    global E_c = E_cp
     

    #v_xp= 4.0*pi*(Utils.integral(grid, (v_x.*density_in.*grid_sqrt)))

    #v_cp= 4.0*pi*(Utils.integral(grid, (v_c.*density_in.*grid_sqrt)))

    v_xc= v_x .+ v_c;
    v_xc= 4.0*pi*(Utils.integral(grid, (v_xc.*density_in.*grid_sqrt)))
    E_hartree= 0.5*4.0*pi*(Utils.integral(grid, (v_h.*density_in.*grid_sqrt)))
    E_x= 4.0*pi*(Utils.integral(grid, (E_x.*density_in.*grid_sqrt)))
    E_c= 4.0*pi*(Utils.integral(grid, (E_c.*density_in.*grid_sqrt)))


    global E_total_before= float(E_total)
    global E_total=0.0
    E_eigen=0.0
    for i in axes(temp_u_basis_set,1)
        E_eigen += temp_u_basis_set[i]["E"]*temp_u_basis_set[i]["occu"]

        #uu= temp_u_basis_set[i]["u"].^2.0
        #display(plot(grid[1:7000],[uu[1:7000],temp_u_basis_set[i]["u"][1:7000]]))
    end    
    global step+=1
    global E_total= E_eigen - E_hartree + E_x + E_c - v_xc
    C_in= Utils.integral(grid,(4.0.*pi.*density_in.*grid.^2))
    C_out= Utils.integral(grid,(4.0.*pi.*density_out.*grid.^2))
    E_total= E_eigen - E_hartree + E_x + E_c - v_xc
    println("C_in ", C_in)
    println("C_out ", C_out)
    println("eigen ", eigen_guesses)
    println("E_hartree ", E_hartree)
    println("E_xc ", (E_c + E_x))
    println("E_total ", E_total)
    println("step ", step)
    println("************************************")
end

