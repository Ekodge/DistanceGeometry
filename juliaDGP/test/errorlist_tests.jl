@testset "ErrorList constructors" begin

    # constructor from a DGP instance and one possible Realization instance (a possible solution)
    @testset "constructor from a DGP instance and a Realization instance" begin
        atoms = create_sample_PDBTool_atoms()
        realization = Realization(atoms)
        dgp = DGP(realization)
        errorlist = ErrorList(dgp,realization)
        #Errorlist constructor does not care about arguments order
        errorlist2 = ErrorList(realization,dgp)
        #Atoms dimension
        @test errorlist.K == 3
        @test errorlist2.K == 3
        #Test without Errors
        @test errorlist.distlist == Dict{Tuple{myAtom, myAtom}, Distance}()
        @test errorlist2.distlist == Dict{Tuple{myAtom, myAtom}, Distance}()
        #Test with Errors
        #TODO
    end

    @testset "ErrorList type various overrides and functions" begin

        # overriding Base show function
        @testset "overriding Base show function" begin
            atoms = create_sample_PDBTool_atoms()
            realization = Realization(atoms)
            dgp = DGP(realization)
            errorlist = ErrorList(dgp,realization)
            io = IOBuffer()
            show(io,errorlist)
            output = String(take!(io))
            expected_output = "ErrorList{myAtom} (K = 3 with 0 elements and 0 errors)"
            @test output == expected_output
        end

        # looking for the smallest error
        @testset "minerror" begin
            atoms = create_sample_PDBTool_atoms()
            realization = Realization(atoms)
            dgp = DGP(realization)
            errorlist = ErrorList(dgp,realization)
            output = minerror(errorlist)
            #case : no error
            @test output == nothing
            #case : error
            #TODO
        end

        # looking for the largest error
        @testset "maxerror" begin
            atoms = create_sample_PDBTool_atoms()
            realization = Realization(atoms)
            dgp = DGP(realization)
            errorlist = ErrorList(dgp,realization)
            output = maxerror(errorlist)
            #case : no error
            @test output == nothing
            #case : error
            #TODO
        end

        # computing the mde
        @testset "mde" begin
            atoms = create_sample_PDBTool_atoms()
            realization = Realization(atoms)
            dgp = DGP(realization)
            errorlist = ErrorList(dgp,realization)
            output = mde(errorlist)
            #TODO
            #case : no error
            #@test output == NaN
        end
    end
end