
module Atomic_Schrodinger_Equation

    include("Utils.jl")
    using .Utils

    include("Integrator.jl")
    using .Integrator


    #function v_l(l,x)
    #    return l*(l+1.0)/(2.0*x^2.0)
    #end

    #function v_z(Z,x)
    #    return -1.0*Z/x
    #end

    function v_atom_array(grid,v_hart,v_xchg,v_corr,v_angu, v_ext)
        #out=[ v_xchg[i] + v_corr[i] + v_hart[i] + v_l(l,xi) + v_z(Z,xi) for (i,xi) in enumerate(grid)]
        out=[ v_xchg[i] + v_corr[i] + v_hart[i] + v_angu(l,xi) + v_ext(Z,xi) for (i,xi) in enumerate(grid)]
    end

    function do_integration(grid, f1, l, E, forward)
        
        N=size(grid)[1]

        y0=zeros(Float32, N)
        y1=zeros(Float32, N)
        g0=zeros(Float32, 6)
        g1=zeros(Float32, 6)
        f0=ones(Float32, 6)

        if forward
            #set initial conditions
            y0[1]= grid[1]^(l+1.0)
            y1[1]= (grid[2]^(l+1.0) - y0[1])/(grid[2]-grid[1])
            #integration loop to generate first 5 elements using rk4
            for i in 1:4
                h= grid[i+1] - grid[i]
                y0[i+1], y1[i+1]= Integrator.RK4(g0,f0,g1,f1[i:i+1],y0[i], y1[i],h)
            end
            #integration loop using prediction correction adams moulton degree 5
            for i in 6:N
                h= grid[i] - grid[i-1]
                y0[i], y1[i]= Integrator.PCABM5(g0,f0,g1,f1[i-5:i],
                                                        y0[i-5:i-1], y1[i-5:i-1],h)
            end
            return y0, y1          
        else
            if sign(E) < 0.0
                lambda= (-2.0*E)^0.5
            end
            for i in 1:N
                temp = grid[i]*exp(-1.0*lambda*grid[i])
                if temp < 10.0e-100 || i < 6
                    y0[i]= temp
                    y1[i]= (grid[i+1]*exp(-1.0*lambda*grid[i+1]) - y0[i])/(grid[i+1]-grid[i])
                else
                    h= grid[i] - grid[i-1]
                    y0[i], y1[i]= Integrator.PCABM5(g0,f0,g1,f1[i-5:i],
                                                            y0[i-5:i-1], y1[i-5:i-1],h)

                end
            end
            return y0, y1
        end

        #integration loop to generate first 4 elements using rk4
        #for i in 1:3
        #integration loop to generate first 5 elements using rk4
        #for i in 1:4
        #    h= grid[i+1] - grid[i]
        #    y0[i+1], y1[i+1]= Integrator.RK4(g0,f0,g1,f1[i:i+1],y0[i], y1[i],h)
        #end
        #integration loop using prediction correction adams moulton degree 4
        #for i in 5:N
        #    h= grid[i] - grid[i-1]
        #    y0[i], y1[i]= Integrator.PCAM4(g0,f0,g1,f1[i-4:i],
        #                                            y0[i-4:i-1], y1[i-4:i-1],h)
        #end
        #integration loop using prediction correction adams moulton degree 5
        #for i in 6:N
        #    h= grid[i] - grid[i-1]
        #    y0[i], y1[i]= Integrator.PCABM5(g0,f0,g1,f1[i-5:i],
        #                                            y0[i-5:i-1], y1[i-5:i-1],h)
        #end
        

        #return y0, y1

    end

    function get_f1(grid,v_hart,v_xchg,v_corr,v_angu, v_ext,E)
        f1=[2.0*(v_xchg[i] + v_corr[i] + v_hart[i] + v_l(l,xi) + v_z(Z,xi) - E) for (i,xi) in enumerate(grid)]
        return f1
    end

    function integrate_SE(grid,grid_bwrd,v_hart,v_xchg,v_corr,v_angu, v_ext,E)
        #println("E ", E)
        #println("l ", l)
        f1=get_f1(grid,v_hart,v_xchg,v_corr,v_angu, v_ext,E)#[2.0*(v_xchg[i] + v_corr[i] + v_hart[i] + v_l(l,xi) + v_z(Z,xi) - E) for (i,xi) in enumerate(grid)]
        node= Utils.find_index_function_nodes(f1)
        #println("node ", node)
        if l == 0
            node= node[1]
        else
            #node= node[2]
            node= node[1]
        end
        #forward Integration
        u_fwrd, _=do_integration(grid, f1,l, E, true)
        #print(u_fwrd[1:100])
        #backward Integration
        u_bwrd, _=do_integration(grid_bwrd, reverse(f1), l, E, false)
        u_bwrd= reverse(u_bwrd)#now in forward order
        #print(u_bwrd[1:100])
        #scale u_fwrd or u_bwrd for merging
        Af=u_fwrd[node]
        Ab=u_bwrd[node]
        if abs(Ab) > abs(Af)
            u_bwrd=[(Af/Ab)*ui for ui in u_bwrd]
        else
            u_fwrd=[(Ab/Af)*ui for ui in u_fwrd]
        end
        #merge value for eigenvalue search
        deri_b= Utils.three_point_derivative(u_bwrd[node-1:node+1],grid[node-1:node+1])
        deri_f= Utils.three_point_derivative(u_fwrd[node-1:node+1],grid[node-1:node+1])
        merge_value= deri_b - deri_f
        #merge forward and backward integration
        N=size(grid)[1]
        u=zeros(Float32, N)
        for i in 1:node-3
            u[i]=u_fwrd[i]
        end
        for i in node-2:node+2
            u[i]= 0.5*(u_fwrd[i] + u_bwrd[i])
        end
        for i in node+3:N
            u[i]=u_bwrd[i]
        end
        return u, u_bwrd, u_fwrd, merge_value, node
    end

    function illinois_eigen_finder(nodes,grid,grid_bwrd,v_hart,v_xchg,v_corr,v_angu, v_ext,N_max=300, tolerance=10.0e-10)
        i=0
        Ec_befo=10.0e2
        Ea=nodes[1]
        Eb=nodes[2]
        Ec=0.0
        _, _, _, u0a, _= integrate_SE(grid,grid_bwrd,v_hart,v_xchg,v_corr,v_angu, v_ext,Ea)
        _, _, _, u0b, _= integrate_SE(grid,grid_bwrd,v_hart,v_xchg,v_corr,v_angu, v_ext,Eb)
        while i < N_max
            Ec=(Ea*u0b -Eb*u0a)/(u0b - u0a)
            if abs(Ec-Ec_befo) < tolerance
                break
            end
            _, _, _, u0c, _= integrate_SE(grid,grid_bwrd,v_hart,v_xchg,v_corr,v_angu, v_ext,Ec)
            if Integer(sign(u0c)) == Integer(sign(u0a))
                Ea=float(Ec)
                u0a=float(u0c)
                u0b=0.5*u0b
            else
                Eb=float(Ec)
                u0b=float(u0c)
                u0a=0.5*u0a
            end
            Ec_befo=Ec
            i+=1
        end
        u, ub, uf, _, node= integrate_SE(grid,grid_bwrd,v_hart,v_xchg,v_corr,v_angu, v_ext,Ec)
        u= Utils.normalize(grid, u)
        return u, ub, uf, node, Ec
    end
end