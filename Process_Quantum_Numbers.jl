module Process_Quantum_Numbers

    function get_shells_occupations_and_quantum_numbers(basis_set_size, quan_numb_data_frame)
        """returns the quantun numbers n, l and occupation by shell"""
        comu_occu=0 #comulative occupancy
        occupation_by_shell=[]
        l_by_shell=[]
        n_by_shell=[]
        remain= Integer(basis_set_size)
        for row in eachrow(quan_numb_data_frame)
            if remain - row.number_of_electrons >= 0
                append!(occupation_by_shell,row.number_of_electrons)
                append!(l_by_shell,row.l)
                append!(n_by_shell,row.n)
            else
                if remain > 0
                    append!(occupation_by_shell,remain)
                    append!(l_by_shell,row.l)
                    append!(n_by_shell, row.n)
                end
            end
            remain-= row.number_of_electrons
            comu_occu+=row.number_of_electrons
            if comu_occu >= basis_set_size
                break
            end
        end
        return occupation_by_shell, n_by_shell, l_by_shell
    end
    
end