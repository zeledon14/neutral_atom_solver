module Math_Utils

    """
    turning_points_indices(func::Vector{Float32})::Vector{Int32}

    Return the indices where the function func changes sign, the index 
    returned is the one before the change in sign. As func is of the form
    V_effe - E, the returned indices mark the classical turning points.

    **Input:**
        -func::Vector{Float32} the function to find the turning points
    
    **Output:**
        -indi::Vector{Int32} a vector with the indices of the turning points
            such that func[indi[i]] has a different sign than func[indi[i+1]]
    """
    function turning_points_indices(func::Vector{Float32})::Vector{Int32}
            indi= Int32[]
            sign_befo= Integer(sign(func[1]))
            for (i,f_i) in enumerate(func[2:end])
                if sign_befo != Integer(sign(f_i))
                    append!(temp, i)
                    sign_befo = Integer(sign(f_i))
                end
            end
            return indi
            
    end

end