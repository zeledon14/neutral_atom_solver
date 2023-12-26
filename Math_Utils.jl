using AutomaticDocstrings
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
            append!(indi, i)
            sign_befo = Integer(sign(f_i))
        end
    end
    return indi
        
end

"""
    rescale!(solution1::Vector{Float32}, 
                   solution2::Vector{Float32}, 
                   turning_point::Int32)::Tuple{Vector{Float32},Vector{Float32}}
rescales the absolute biggest function to the smallest function
suth that solution1(turning_point) = solution2(turning_point). 

**Inputs:**
- `solution1::Vector{Float32}`: Vector with the values of the function1 to be scaled
- `solution2::Vector{Float32}`: Vector with the values of the function2 to be scaled
- `turning_point::Int32`: Point at wich v_effe - E = 0, there may be many of this points
                          the one use is lower.
**Output:**
- `solution1, solution2::Tuple{Vector{Float32},Vector{Float32}}`: Rescaled solutions.
"""
function  rescale!(solution1::Vector{Float32}, 
                   solution2::Vector{Float32}, 
                   turning_point::Int32)::Tuple{Vector{Float32},Vector{Float32}}
    A1=solution1[turning_point];
    A2=solution2[turning_point];
    if  abs(A1) > abs(A2)
        solution1= (A2/A1).*solution1
    else
        solution2= (A1/A2).*solution2
    end
    return solution1, solution2
end

"""
    three_point_derivative(func::Vector{Float32}, grid::Vector{Float32}, turning_point::Float32)

calculate the derivative of func over grid at the point turning_point,
the derivative uses a three point derivative.

**Inputs:**
- `func::Vector{Float32}`: the function for which the derivative is calculated
- `grid::Vector{Float32}`: the grid over which the function is defined
- `turning_point::Float32`: the point where the derivative is calculated
"""
function three_point_derivative(func::Vector{Float32},
                                grid::Vector{Float32},
                                turning_point::Float32)::Float32
    nume= func[turning_point+1] - func[turning_point-1]
    deno= (grid[turning_point+1] - grid[turning_point]) + (grid[turning_point] - grid[turning_point-1])
    temp= nume/deno
    return temp
end

"""
    merge_solutions(forward::Vector{Float32}, backward::Vector{Float32}, turning_point::Float32)

merges the backward and forward integrated solutions into a merged one,
the merging happens at the turning_point.

**Inputs:**
- `forward::Vector{Float32}`: vextor with the forward integrated solution.
- `backward::Vector{Float32}`: vector with the backward integreated solution.
- `turning_point::Float32`: turning point where the solutions are merged.

**Inputs:**
- `u_merged::Vector{Float32}`: vector with the merged function.
"""
function merge_solutions(forward::Vector{Float32}, backward::Vector{Float32}, 
                         turning_point::Int32)::Vector{Float32}
    u_merged=zeros(Float32, size(forward)[1]);
    u_merged[1:turning_point-2]= forward[1:turning_point-2];
    u_merged[turning_point-1:turning_point+1]= (0.5).*(forward[turning_point-1:turning_point+1] 
                                                       .+ backward[turning_point-1:turning_point+1])
    u_merged[turning_point+2:end]= backward[turning_point+2:end];
    return u_merged
end
end
