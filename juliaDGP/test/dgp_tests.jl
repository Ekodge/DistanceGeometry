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

    @testset "the different equals operators" begin

        @testset "override of the '==' operator" begin
            atoms = create_sample_PDBTool_atoms()
            realization = Realization(atoms)
            dgp = DGP(realization)
            dgp2 = DGP(realization)

            #Case : the two instances are the same
            @test dgp == dgp2

            #Case : the two instances have different dimensions
            atoms = [Atom(index=1, index_pdb=1, name="N", resname="ALA", chain="A", resnum=1, x=1.0, y=2.0, z=3.0),
                 Atom(index=2, index_pdb=2, name="CA", resname="CYS", chain="A", resnum=2, x=4.0, y=5.0, z=6.0),
                 Atom(index=3, index_pdb=3, name="C", resname="GLY", chain="A", resnum=3, x=7.0, y=8.0, z=9.0),
                 Atom(index=4, index_pdb=4, name="O", resname="GLY", chain="A", resnum=4, x=10.0, y=11.0, z=12.0)]
            realization2 = Realization(atoms)
            dgp2 = DGP(realization2)
            @test dgp != dgp2

            #Case : the two instances have different distances
            atoms = [Atom(index=1, index_pdb=1, name="N", resname="ALA", chain="A", resnum=1, x=7.0, y=2.5, z=5.0),
                 Atom(index=2, index_pdb=2, name="CA", resname="CYS", chain="A", resnum=2, x=16.0, y=5.0, z=5.0),
                 Atom(index=3, index_pdb=3, name="C", resname="GLY", chain="A", resnum=3, x=17.0, y=17.0, z=9.0)]
            realization3 = Realization(atoms)
            dgp2 = DGP(realization3)
            @test dgp != dgp2

            #Case : the two instances have different atoms but same dimensions,structures and distances
            atoms = [Atom(index=1, index_pdb=1, name="A", resname="ALA", chain="A", resnum=1, x=1.0, y=2.0, z=3.0),
            Atom(index=2, index_pdb=2, name="B", resname="CYS", chain="A", resnum=2, x=4.0, y=5.0, z=6.0),
            Atom(index=3, index_pdb=3, name="C", resname="GLY", chain="A", resnum=3, x=7.0, y=8.0, z=9.0)]
            realization3 = Realization(atoms)
            dgp2 = DGP(realization3)
            @test dgp == dgp2

            #Case : the two instances does not have the same order
            atoms=[Atom(index=1, index_pdb=1, name="N", resname="ALA", chain="A", resnum=1, x=1.0, y=2.0, z=3.0),
                   Atom(index=2, index_pdb=2, name="C", resname="ALA", chain="A", resnum=1, x=4.0, y=5.0, z=6.0),
                   Atom(index=3, index_pdb=3, name="CA", resname="ALA", chain="A", resnum=1, x=7.0, y=8.0, z=9.0)]
            realization4 = Realization(atoms)
            dgp2 = DGP(realization4)
            @test dgp == dgp2
        end

        @testset " (===) equalOrder function" begin
            atoms = create_sample_PDBTool_atoms()
            realization = Realization(atoms)
            dgp = DGP(realization)
            dgp2 = DGP(realization)

            #Case : the two instances are the same
            @test equalOrder(dgp,dgp2)

            #Case : the two instances have different dimensions
            atoms = [Atom(index=1, index_pdb=1, name="N", resname="ALA", chain="A", resnum=1, x=1.0, y=2.0, z=3.0),
                 Atom(index=2, index_pdb=2, name="CA", resname="CYS", chain="A", resnum=2, x=4.0, y=5.0, z=6.0),
                 Atom(index=3, index_pdb=3, name="C", resname="GLY", chain="A", resnum=3, x=7.0, y=8.0, z=9.0),
                 Atom(index=4, index_pdb=4, name="O", resname="GLY", chain="A", resnum=4, x=10.0, y=11.0, z=12.0)]
            realization2 = Realization(atoms)
            dgp2 = DGP(realization2)
            @test notEqualOrder(dgp,dgp2)

            #Case : the two instances have different distances
            atoms = [Atom(index=1, index_pdb=1, name="N", resname="ALA", chain="A", resnum=1, x=7.0, y=2.5, z=5.0),
                 Atom(index=2, index_pdb=2, name="CA", resname="CYS", chain="A", resnum=2, x=16.0, y=5.0, z=5.0),
                 Atom(index=3, index_pdb=3, name="C", resname="GLY", chain="A", resnum=3, x=17.0, y=17.0, z=9.0)]
            realization3 = Realization(atoms)
            dgp2 = DGP(realization3)
            @test notEqualOrder(dgp,dgp2)

            #Case : the two instances have different atoms but same dimensions,structures and distances
            atoms = [Atom(index=1, index_pdb=1, name="A", resname="ALATTEST", chain="A", resnum=1, x=1.0, y=2.0, z=3.0),
            Atom(index=2, index_pdb=2, name="B", resname="CYSTEST", chain="A", resnum=2, x=4.0, y=5.0, z=6.0),
            Atom(index=3, index_pdb=3, name="C", resname="GLYTEST", chain="A", resnum=3, x=7.0, y=8.0, z=9.0)]
            realization3 = Realization(atoms)
            dgp2 = DGP(realization3)
            @test equalOrder(dgp,dgp2)

            #Case : the two instances does not have the same order
            atoms=[Atom(index=1, index_pdb=1, name="N", resname="ALA", chain="A", resnum=1, x=1.0, y=2.0, z=3.0),
                   Atom(index=3, index_pdb=3, name="CA", resname="ALA", chain="A", resnum=1, x=7.0, y=8.0, z=9.0),
                   Atom(index=2, index_pdb=2, name="C", resname="ALA", chain="A", resnum=1, x=4.0, y=5.0, z=6.0)]
            realization4 = Realization(atoms)
            dgp2 = DGP(realization4)
            @test notEqualOrder(dgp,dgp2)
        end

        @testset " (====) equals function" begin
            atoms = create_sample_PDBTool_atoms()
            realization = Realization(atoms)
            dgp = DGP(realization)
            dgp2 = DGP(realization)

            #Case : the two instances are the same
            @test equals(dgp,dgp2)

            #Case : the two instances have different dimensions
            atoms = [Atom(index=1, index_pdb=1, name="N", resname="ALA", chain="A", resnum=1, x=1.0, y=2.0, z=3.0),
                 Atom(index=2, index_pdb=2, name="CA", resname="CYS", chain="A", resnum=2, x=4.0, y=5.0, z=6.0),
                 Atom(index=3, index_pdb=3, name="C", resname="GLY", chain="A", resnum=3, x=7.0, y=8.0, z=9.0),
                 Atom(index=4, index_pdb=4, name="O", resname="GLY", chain="A", resnum=4, x=10.0, y=11.0, z=12.0)]
            realization2 = Realization(atoms)
            dgp2 = DGP(realization2)
            @test notEquals(dgp,dgp2)

            #Case : the two instances have different distances
            atoms = [Atom(index=1, index_pdb=1, name="N", resname="ALA", chain="A", resnum=1, x=7.0, y=2.5, z=5.0),
                 Atom(index=2, index_pdb=2, name="CA", resname="CYS", chain="A", resnum=2, x=16.0, y=5.0, z=5.0),
                 Atom(index=3, index_pdb=3, name="C", resname="GLY", chain="A", resnum=3, x=17.0, y=17.0, z=9.0)]
            realization3 = Realization(atoms)
            dgp2 = DGP(realization3)
            @test notEquals(dgp,dgp2)

            #Case : the two instances have different atoms but same dimensions,structures and distances
            atoms = [Atom(index=1, index_pdb=1, name="A", resname="ALATTEST", chain="A", resnum=1, x=1.0, y=2.0, z=3.0),
            Atom(index=2, index_pdb=2, name="B", resname="CYSTEST", chain="A", resnum=2, x=4.0, y=5.0, z=6.0),
            Atom(index=3, index_pdb=3, name="C", resname="GLYTEST", chain="A", resnum=3, x=7.0, y=8.0, z=9.0)]
            realization3 = Realization(atoms)
            dgp2 = DGP(realization3)
            @test notEquals(dgp,dgp2)

            #Case : the two instances does not have the same order
            atoms=[Atom(index=1, index_pdb=1, name="N", resname="ALA", chain="A", resnum=1, x=1.0, y=2.0, z=3.0),
                   Atom(index=3, index_pdb=3, name="CA", resname="ALA", chain="A", resnum=1, x=7.0, y=8.0, z=9.0),
                   Atom(index=2, index_pdb=2, name="C", resname="ALA", chain="A", resnum=1, x=4.0, y=5.0, z=6.0)]
            realization4 = Realization(atoms)
            dgp2 = DGP(realization4)
            @test notEquals(dgp,dgp2)
        end

    end
end