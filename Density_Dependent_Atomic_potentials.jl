module Density_Dependent_Atomic_potentials
    
    include("Atomic_Poisson_Equation.jl")
    using .Atomic_Poisson_Equation: Atomic_Poisson_Equation as APE

    include("Exchange_Correlation.jl")
    using .Exchange_Correlation: Exchange_Correlation as EC

    function density_potentials(grid, density, Z)
        V_hartree= APE.integrate_PE(grid,density,Z)
        V_x, E_x, V_c, E_c= EC.Exchange_Correlation_Potentials(density)
        return V_hartree, V_x, E_x, V_c, E_c
    end

end