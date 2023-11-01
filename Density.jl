module Density

    function  calculate_density(grid, u_basis_set)
        density= zero(grid)
        for i in axes(u_basis_set,1)
            i_bs=u_basis_set[i]
            u=i_bs["u"]
            c=i_bs["occu"]/(4.0*pi)
            density .= density .+ c.*(u.^2.0)./(grid.^2.0)
            #for (j,xj) in enumerate(grid)
            #    density[j]= density[j] + c*(u[j]^2.0)/(xj^2.0)
            #end
        end
        return density
    end    

    function linear_mixing(density_in, density_out; alpha=0.3)
        temp1= alpha.*density_out
        temp2= (1.0 - alpha).*density_in
        return temp1 .+ temp2
    end


end