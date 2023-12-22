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
end