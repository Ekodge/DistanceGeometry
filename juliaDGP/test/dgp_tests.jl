#defining the necessary atoms
atoms = create_sample_PDBTool_atoms()
atom1 = myAtom("N", 1, "ALA")
atom2 = myAtom("CA", 1, "ALA")
atom3 = myAtom("C", 1, "ALA")

@testset "dgp various constructors" begin

    @testset "constructor for random fully-connected instance from given Realization" begin
        atoms = create_sample_PDBTool_atoms()
        realization = Realization(atoms)
        dgp = DGP(realization)
        #Check dimension
        @test dgp.K == 3
        #Check the vertex order
        @test all(x -> x in dgp.order, [myAtom("N", 1, "ALA"), myAtom("CA", 1, "ALA"), myAtom("C", 1, "ALA")])
        #Check the distances
        @test length(dgp.distlist) == 3
        @test dgp.distlist[(atom2, atom1)] == Distance([4.0, 5.0, 6.0], [1.0, 2.0, 3.0])
        @test dgp.distlist[(atom3, atom2)] == Distance([7.0, 8.0, 9.0], [4.0, 5.0, 6.0])
        @test dgp.distlist[(atom3, atom1)] == Distance([7.0, 8.0, 9.0], [1.0, 2.0, 3.0])
    end

    @testset "constructor from a given Realization instance and distance cut-off" begin
        atoms = create_sample_PDBTool_atoms()
        realization = Realization(atoms)
        dgp = DGP(realization,8.0)
        #Error : The specified threshold is 0.0 or even negative
        @test_throws ArgumentError DGP(realization,0.0)
        @test_throws ArgumentError DGP(realization,-8.0)
        #Check dimension
        @test dgp.K == 3
        #Check the vertex order
        @test all(x -> x in dgp.order, [myAtom("N", 1, "ALA"), myAtom("CA", 1, "ALA"), myAtom("C", 1, "ALA")])
        #Check the distances
        @test length(dgp.distlist) == 2
        @test dgp.distlist[(atom2, atom1)] == Distance([4.0, 5.0, 6.0], [1.0, 2.0, 3.0])
        @test dgp.distlist[(atom3, atom2)] == Distance([7.0, 8.0, 9.0], [4.0, 5.0, 6.0])
        @test_throws KeyError dgp.distlist[(atom3, atom1)]
    end

    @testset "constructor for a paradoxical DGP from a given Realization instance" begin
        atoms = create_sample_PDBTool_atoms()
        realization = Realization(atoms)
        #Case : cut = realization size
        dgp = DGP(realization,3)
        #Error : The specified cut value is 0 or even negative
        @test_throws ArgumentError DGP(realization,0)
        #Error : The value of cut cannot be larger than the Realization size
        @test_throws ArgumentError DGP(realization,6)
        #Check dimension
        @test dgp.K == 3
        #Check the vertex order
        @test all(x -> x in dgp.order, [myAtom("N", 1, "ALA"), myAtom("CA", 1, "ALA"), myAtom("C", 1, "ALA")])
        #Check the distances
        @test dgp.distlist[(atom2, atom1)] == Distance([4.0, 5.0, 6.0], [1.0, 2.0, 3.0])
        @test dgp.distlist[(atom3, atom2)] == Distance([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])
        @test dgp.distlist[(atom3, atom1)] == Distance([7.0, 8.0, 9.0], [1.0, 2.0, 3.0])
        #Case : cut < realization size
        dgp = DGP(realization,1)
        #Check dimension
        @test dgp.K == 3
        #Check the vertex order
        @test all(x -> x in dgp.order, [myAtom("N", 1, "ALA"), myAtom("CA", 1, "ALA"), myAtom("C", 1, "ALA")])
        #Check the distances
        @test length(dgp.distlist) == 3
        @test dgp.distlist[(atom1, atom2)] == Distance([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])
        @test dgp.distlist[(atom2, atom3)] == Distance([4.0, 5.0, 6.0], [7.0, 8.0, 9.0])
        @test dgp.distlist[(atom3, atom1)] == Distance([7.0, 8.0, 9.0], [1.0, 2.0, 3.0])
    end

    @testset "constructor from a given Realization instance and its connection graph" begin
        atoms = create_sample_PDBTool_atoms()
        g = UndirectedGraph{myAtom}()
        realization = Realization(atoms)
        add!(g,atom2,atom1)
        add!(g,atom3,atom2)
        #Case : Didn't specify all the edges
        dgp = DGP(realization,g)
        #Check dimension
        @test dgp.K == 3
        #Check the vertex order
        @test all(x -> x in dgp.order, [myAtom("N", 1, "ALA"), myAtom("CA", 1, "ALA"), myAtom("C", 1, "ALA")])
        #Check the distances
        @test length(dgp.distlist) == 2
        @test dgp.distlist[(atom2, atom1)] == Distance([4.0, 5.0, 6.0], [1.0, 2.0, 3.0])
        @test dgp.distlist[(atom3, atom2)] == Distance([7.0, 8.0, 9.0], [4.0, 5.0, 6.0])
        @test_throws KeyError dgp.distlist[(atom3, atom1)]
        #Case : specified all the edges
        add!(g,atom1,atom3)
        dgp = DGP(realization,g)
        #Check dimension
        @test dgp.K == 3
        #Check the vertex order
        @test all(x -> x in dgp.order, [myAtom("N", 1, "ALA"), myAtom("CA", 1, "ALA"), myAtom("C", 1, "ALA")])
        #Check the distances
        @test length(dgp.distlist) == 3
        @test dgp.distlist[(atom2, atom1)] == Distance([4.0, 5.0, 6.0], [1.0, 2.0, 3.0])
        @test dgp.distlist[(atom3, atom2)] == Distance([7.0, 8.0, 9.0], [4.0, 5.0, 6.0])
        @test dgp.distlist[(atom3, atom1)] == Distance([7.0, 8.0, 9.0], [1.0, 2.0, 3.0])
    end

    @testset "constructor from STAR file" begin
        prefix = "test/fileContainer/"
        atoms = create_sample_PDBTool_atoms()
        realization = Realization(atoms)
        dgp = DGP(realization)
        #Check dimension
        @test dgp.K == 3
        #Check the vertex order
        @test all(x -> x in dgp.order, [myAtom("N", 1, "ALA"), myAtom("CA", 1, "ALA"), myAtom("C", 1, "ALA")])
        #Check the distances
        @test length(dgp.distlist) == 3
        @test dgp.distlist[(atom2, atom1)] == Distance([4.0, 5.0, 6.0], [1.0, 2.0, 3.0])
        @test dgp.distlist[(atom3, atom2)] == Distance([7.0, 8.0, 9.0], [4.0, 5.0, 6.0])
        @test dgp.distlist[(atom3, atom1)] == Distance([7.0, 8.0, 9.0], [1.0, 2.0, 3.0])
        #Create the star file corresponding to dgp
        starfile = prefix * "test_output.str"
        create_pseudo_STAR(dgp,starfile)
        #Create a dgp using the created file
        dgpFromStarfile = DGP(starfile)
        #Check dimension
        @test dgpFromStarfile.K == 3
        #Check the vertex order
        @test all(x -> x in dgpFromStarfile.order, [myAtom("N", 1, "ALA"), myAtom("CA", 1, "ALA"), myAtom("C", 1, "ALA")])
        #Check the distances
        @test length(dgpFromStarfile.distlist) == 3
        @test dgpFromStarfile.distlist[(atom2, atom1)] == Distance([4.0, 5.0, 6.0], [1.0, 2.0, 3.0])
        @test dgpFromStarfile.distlist[(atom3, atom2)] == Distance([7.0, 8.0, 9.0], [4.0, 5.0, 6.0])
        @test dgpFromStarfile.distlist[(atom3, atom1)] == Distance([7.0, 8.0, 9.0], [1.0, 2.0, 3.0])
    end
end

@testset "DGP type various overrides and functions" begin

    @testset "overriding Base show function" begin
        atoms = create_sample_PDBTool_atoms()
        realization = Realization(atoms)
        dgp = DGP(realization)
        io = IOBuffer()
        show(io,dgp)
        output = String(take!(io))
        expected_output = "DGP{myAtom} (K = 3 with 3 elements and 3 distances)"
        @test output == expected_output
    end

    @testset "details" begin
        atoms = create_sample_PDBTool_atoms()
        realization = Realization(atoms)
        dgp = DGP(realization)
        @test print_to_string_stdout(details,dgp) == "DGP{myAtom} (K = 3) {\n  (ALA1-CA,ALA1-N) => Distance(5.196152422706632)\n  (ALA1-C,ALA1-N) => Distance(10.392304845413264)\n  (ALA1-C,ALA1-CA) => Distance(5.196152422706632)\n}"
    end

    @testset "references" begin
        atoms = create_sample_PDBTool_atoms()
        realization = Realization(atoms)
        dgp = DGP(realization)
        #Case : the given element does not belong to the DGP instance
        @test_throws ArgumentError references(dgp,myAtom("O", 1, "ALA"))
        #Case : the given element belongs to the DGP instance
        @test references(dgp,myAtom("C", 1, "ALA")) == [myAtom("N", 1, "ALA"), myAtom("CA", 1, "ALA")]
        @test references(dgp,myAtom("CA", 1, "ALA")) == [myAtom("N", 1, "ALA")]
        @test references(dgp,myAtom("N", 1, "ALA")) == []
    end

    @testset "MDjeep_format and DGP from given data : K (dimension), distlist and order" begin
        atoms = create_sample_PDBTool_atoms()
        realization = Realization(atoms)
        dgp = DGP(realization)
        #Check dimension
        @test dgp.K == 3
        #Check the vertex order
        @test all(x -> x in dgp.order, [myAtom("N", 1, "ALA"), myAtom("CA", 1, "ALA"), myAtom("C", 1, "ALA")])
        #Check the distances
        @test length(dgp.distlist) == 3
        @test dgp.distlist[(atom2, atom1)] == Distance([4.0, 5.0, 6.0], [1.0, 2.0, 3.0])
        @test dgp.distlist[(atom3, atom2)] == Distance([7.0, 8.0, 9.0], [4.0, 5.0, 6.0])
        @test dgp.distlist[(atom3, atom1)] == Distance([7.0, 8.0, 9.0], [1.0, 2.0, 3.0])
        suppress_stdout(MDjeep_format,dgp)
        dgpFromMDjeepFormat = create_DGP_from_MDjeep("test_instance.nmr")
        #Check dimension
        @test dgpFromMDjeepFormat.K == 3
        #Check the vertex order
        @test all(x -> x in dgpFromMDjeepFormat.order, [myAtom("N", 1, "ALA"), myAtom("CA", 1, "ALA"), myAtom("C", 1, "ALA")])
        #Check the distances
        @test length(dgpFromMDjeepFormat.distlist) == 3
        @test dgpFromMDjeepFormat.distlist[(atom2, atom1)] == Distance([4.0, 5.0, 6.0], [1.0, 2.0, 3.0])
        @test dgpFromMDjeepFormat.distlist[(atom3, atom2)] == Distance([7.0, 8.0, 9.0], [4.0, 5.0, 6.0])
        @test dgpFromMDjeepFormat.distlist[(atom3, atom1)] == Distance([7.0, 8.0, 9.0], [1.0, 2.0, 3.0])
        #Delete the created file
        rm("test_instance.nmr")
    end
end