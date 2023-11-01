module Hydrogen

    function u(grid)
        u=[xi*exp(-1.0*xi) for xi in grid]
        return u 
    end

    function U_hartree(grid)
        out=[(1.0 -1.0*(xi+1.0)*exp(-2.0*xi))  for (i,xi) in enumerate(grid)]
    end
end