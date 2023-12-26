module Tests

    include("../Grids.jl");
    using .Grids;
    include("../Potentials.jl")
    using .Potentials;
    include("../Hydrogen.jl")
    using .Hydrogen;


    function test_hydrogen_u_s1_integration()
        return true
        
    end

#testing Math_Utils.merge_solutions
frwd::Vector{Float32}=[1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0];
brwd::Vector{Float32}=[9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0];
#expected result
resu::Vector{Float32}=[1.0,2.0,7.0,8.0,9.0,14.0,15.0,16.0];
point::Int32=4
u_merged= Math_Utils.merge_solutions(frwd, brwd, point)
end