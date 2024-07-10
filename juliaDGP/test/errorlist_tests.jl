@testset "ErrorList constructors" begin

    #defining the necessary atoms
    atoms = create_sample_PDBTool_atoms()
    atom1 = myAtom("N", 1, "ALA")
    atom2 = myAtom("CA", 1, "ALA")
    atom3 = myAtom("C", 1, "ALA")

    # constructor from a DGP instance and one possible Realization instance (a possible solution)
    @testset "constructor from a DGP instance and a Realization instance" begin
        realization = Realization(atoms)
        dgp = DGP(realization)
        errorlist = ErrorList(dgp,realization)
        #Errorlist constructor does not care about arguments order
        errorlist2 = ErrorList(realization,dgp)
        #Atoms dimension
        @test errorlist.K == 3
        @test errorlist2.K == 3

        #Case : no error
        @test errorlist.distlist == Dict{Tuple{myAtom, myAtom}, Distance}()
        @test errorlist2.distlist == Dict{Tuple{myAtom, myAtom}, Distance}()

        #Case : errors
        errorlist.distlist[(atom1,atom2)] = Distance(1.0)
        errorlist.distlist[(atom1,atom3)] = Distance(8.0)
        errorlist2.distlist[(atom1,atom2)] = Distance(1.0)
        errorlist2.distlist[(atom1,atom3)] = Distance(8.0)
        expectedDictlist = Dict{Tuple{myAtom, myAtom}, Distance}((atom1,atom2) => Distance(1.0), (atom1,atom3) => Distance(8.0))
        @test errorlist.distlist == expectedDictlist
        @test errorlist2.distlist == expectedDictlist
    end

    @testset "ErrorList type various overrides and functions" begin

        # overriding Base show function
        @testset "overriding Base show function" begin
            realization = Realization(atoms)
            dgp = DGP(realization)
            errorlist = ErrorList(dgp,realization)
            io = IOBuffer()
            show(io,errorlist)
            output = String(take!(io))

            #Case : no error
            expected_output = "ErrorList{myAtom} (K = 3 with 0 elements and 0 errors)"
            @test output == expected_output

            #Case : errors
            errorlist.distlist[(atom1,atom2)] = Distance(1.0)
            errorlist.distlist[(atom1,atom3)] = Distance(8.0)
            io = IOBuffer()
            show(io,errorlist)
            output = String(take!(io))
            expected_output = "ErrorList{myAtom} (K = 3 with 0 elements and 2 errors)"
        end

        # looking for the largest error / smallest error
        @testset "maxerror minerror" begin
            realization = Realization(atoms)
            dgp = DGP(realization)
            errorlist = ErrorList(dgp,realization)
            output = maxerror(errorlist)

            #Case : no error
            @test output == nothing

            #Case : error
            errorlist.distlist[(atom1,atom2)] = Distance(1.0)
            maxOutput = maxerror(errorlist)
            minOutput = minerror(errorlist)
            @test maxOutput == ((atom1,atom2) => 1.0)
            @test minOutput == ((atom1,atom2) => 1.0)

            #Case : errors
            errorlist.distlist[(atom1,atom3)] = Distance(8.0)
            maxOutput = maxerror(errorlist)
            minOutput = minerror(errorlist)
            @test maxOutput == ((atom1,atom3) => 8.0)
            @test minOutput == ((atom1,atom2) => 1.0)
        end

        # computing the mde
        @testset "mde" begin
            realization = Realization(atoms)
            dgp = DGP(realization)
            errorlist = ErrorList(dgp,realization)
            output = mde(errorlist)

            #Case : no error
            @test isnan(output)

            #Case : errors
            errorlist.distlist[(atom1, atom2)] = Distance(1.0)
            errorlist.distlist[(atom2, atom3)] = Distance(2.0)
            errorlist.distlist[(atom1, atom3)] = Distance(3.0)
            output = mde(errorlist)
            expected_mde = (1.0 + 2.0 + 3.0) / 3  # Sum of errors divided by number of errors
        end
    end
end