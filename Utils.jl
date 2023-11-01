module Utils
    function find_index_function_nodes(func)
        #func: vector containing the values of the function
        temp= Int[]
        sign_befo= Integer(sign(func[1]))
        for (i,yi) in enumerate(func[2:end])
            if sign_befo != Integer(sign(func[i]))
                append!(temp, i)
                sign_befo = Integer(sign(func[i]))
            end
        end
        return temp
    end

    function three_point_derivative(func,grid)
        nume= func[3] - func[1]
        deno= (grid[3] - grid[2]) + (grid[2] - grid[1])
        temp= nume/deno
        return temp
    end

    function integral(grid,func)
        I= 0.5.*(grid[2:end] .- grid[1:end-1]).*(func[2:end] .+ func[1:end-1])
        I= sum(I)
        #I=0.0
        #for (i,xi) in enumerate(grid[1:end-1])
        #    I+= 0.5*(grid[i+1]-xi)*(func[i] + func[i+1])
        #end
        return I
    end

    function normalize(grid,func)
        func_sqrt= func.^2.0
        I= integral(grid, func_sqrt)
        #I=0.0
        #for (i,xi) in enumerate(grid[1:end-1])
        #    I+= 0.5*(grid[i+1]-xi)*(func[i]^2.0 + func[i+1]^2)
        #end 
        I=I^(0.5)
        #out= [temp/I for temp in func]
        out= func./I
        return out
    end
end