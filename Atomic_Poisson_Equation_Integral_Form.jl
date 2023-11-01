module Atomic_Poisson_Equation_Integral_Form

    include("Integrator.jl")
    using .Integrator

    include("Utils.jl")
    using .Utils

    function do_integration(grid, g1, y0_0)
        #integrates U_hartree
        #V_hartree= U_hartree/r
        N=size(grid)[1]

        y0=zeros(Float64, N)
        y1=zeros(Float64, N)
        g0=zeros(Float64, 6)
        f1=zeros(Float64, 6)
        f0=ones(Float64, 6)

        y0[1]= y0_0
        y1[1]= y0_0/grid[1]#0.0
        #integration loop to generate first 4 elements using rk4
        #for i in 1:3
        #integration loop to generate first 5 elements using rk4
        for i in 1:4
            h= grid[i+1] - grid[i]
            y0[i+1], y1[i+1]= Integrator.RK4(g0,f0,g1[i:i+1],f1,y0[i], y1[i],h)
        end
        #integration loop using prediction correction adams moulton degree 4
        #for i in 5:N
        #    h= grid[i] - grid[i-1]
        #    y0[i], y1[i]= Integrator.PCAM4(g0,f0,g1[i-4:i],f1,
        #                                            y0[i-4:i-1], y1[i-4:i-1],h)
        #end
        #integration loop using prediction correction adams moulton degree 5
        for i in 6:N
            h= grid[i] - grid[i-1]
            y0[i], y1[i]= Integrator.PCABM5(g0,f0,g1[i-5:i],f1,
                                                    y0[i-5:i-1], y1[i-5:i-1],h)
        end
        return y0, y1

    end

    function integrate_PE(grid,density,Z)
    
        #g1=[-4.0*pi*density[i]*xi for (i,xi) in enumerate(grid)]
        g1=-4.0.*pi.*density.*grid
        y0_0= -1.0*grid[1]*(Utils.integral(grid,g1))
        #U_hartree integration
        U_hartree,_ =do_integration(grid, g1, y0_0)
        #set boudary condtion over U_hartree
        #a= (Z - U_hartree[end])/grid[end]
        #U_hartree= U_hartree .+ a.*grid
        #for (i,xi) in enumerate(grid)
        #    U_hartree[i]+= a*xi
        #end
        #transform into V_hartree
        V_hartree=U_hartree./grid#[U_hartree[i]/xi for (i,xi) in enumerate(grid)]
        #return U_hartree
        return V_hartree
    end

end