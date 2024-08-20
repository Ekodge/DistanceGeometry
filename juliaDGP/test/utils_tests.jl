@testset "Various functions from Utils.jl" begin

#defining the necessary atoms
    atoms = create_sample_PDBTool_atoms()
    atom1 = myAtom("N", 1, "ALA")
    atom2 = myAtom("CA", 1, "ALA")
    atom3 = myAtom("C", 1, "ALA")

    @testset "common_elements" begin
        #Error : the two instances differ in dimensionality
        realization = Realization(3,3)
        dgp = DGP(realization)
        realization2 = Realization(2,2)
        @test_throws ArgumentError common_elements(dgp, realization2)

        #Common_elements does not care about arguments order
        @test_throws ArgumentError common_elements(realization2, dgp)

        #Case : all atoms are common
        realization = Realization(atoms)
        dgp = DGP(realization)
        @test common_elements(dgp, realization) == [atom1, atom2, atom3]
        @test common_elements(realization, dgp) == [atom1, atom2, atom3]

        #Case : no common elements
        atoms2 = [Atom(index=1, index_pdb=1, name="CB", resname="ALA", chain="A", resnum=1, x=1.0, y=2.0, z=3.0)]
        realization2 = Realization(atoms2)
        @test common_elements(dgp, realization2) == []
        @test common_elements(realization2, dgp) == []

        #Case : one common element
        atoms2 = [Atom(index=1, index_pdb=1, name="N", resname="ALA", chain="A", resnum=1, x=1.0, y=2.0, z=3.0)]
        realization2 = Realization(atoms2)
        @test common_elements(dgp, realization2) == [atom1]
        @test common_elements(realization2, dgp) == [atom1]
    end

    @testset "count_common_elements" begin
        #Case : all atoms are common
        realization = Realization(atoms)
        dgp = DGP(realization)
        @test count_common_elements(dgp, realization) == 3

        #Common_elements does not care about arguments order
        @test count_common_elements(realization, dgp) == 3

        #Case : no common elements
        atoms2 = [Atom(index=1, index_pdb=1, name="CB", resname="ALA", chain="A", resnum=1, x=1.0, y=2.0, z=3.0)]
        realization2 = Realization(atoms2)
        @test count_common_elements(dgp, realization2) == 0
        @test count_common_elements(realization2, dgp) == 0

        #Case : one common element
        atoms2 = [Atom(index=1, index_pdb=1, name="N", resname="ALA", chain="A", resnum=1, x=1.0, y=2.0, z=3.0)]
        realization2 = Realization(atoms2)
        @test count_common_elements(dgp, realization2) == 1
        @test count_common_elements(realization2, dgp) == 1
    end

    @testset "get_extension" begin
        #Case : .pdb extension
        @test get_extension("test.pdb") == "pdb"

        #Case : .txt extension
        @test get_extension("test.txt") == "txt"

        #Case : .xyz extension
        @test get_extension("test.xyz") == "xyz"

        #Case : no extension
        ## TODO : should you return test or nothing or error ?
        #@test get_extension("test") == ""
    end

    @testset "sequence" begin
        #Error : Input Vector{Atom} has zero length
        atoms = Vector{Atom}()
        @test_throws ArgumentError sequence(atoms)
        #Case: three residues with different resnums
        atoms = [Atom(index=1, index_pdb=1, name="N", resname="ALA", chain="A", resnum=1, x=1.0, y=2.0, z=3.0),
                 Atom(index=2, index_pdb=2, name="CA", resname="CYS", chain="A", resnum=2, x=4.0, y=5.0, z=6.0),
                 Atom(index=3, index_pdb=3, name="C", resname="GLY", chain="A", resnum=3, x=7.0, y=8.0, z=9.0)]
        @test sequence(atoms) == "ACG"

        #Case: multiple atoms within the same residue
        atoms = [Atom(index=1, index_pdb=1, name="N", resname="ALA", chain="A", resnum=1, x=1.0, y=2.0, z=3.0),
                 Atom(index=2, index_pdb=2, name="CA", resname="ALA", chain="A", resnum=1, x=4.0, y=5.0, z=6.0),
                 Atom(index=3, index_pdb=3, name="C", resname="ALA", chain="A", resnum=1, x=7.0, y=8.0, z=9.0)]
        @test sequence(atoms) == "A"

        #Case: mixed residues with multiple atoms and different order
        atoms = [Atom(index=1, index_pdb=1, name="N", resname="ALA", chain="A", resnum=2, x=1.0, y=2.0, z=3.0),
                 Atom(index=3, index_pdb=3, name="N", resname="CYS", chain="A", resnum=3, x=7.0, y=8.0, z=9.0),
                 Atom(index=2, index_pdb=2, name="CA", resname="ALA", chain="A", resnum=2, x=4.0, y=5.0, z=6.0),
                 Atom(index=4, index_pdb=4, name="CA", resname="CYS", chain="A", resnum=3, x=10.0, y=11.0, z=12.0),
                 Atom(index=5, index_pdb=5, name="C", resname="GLY", chain="A", resnum=1, x=13.0, y=14.0, z=15.0)]
        @test sequence(atoms) == "GAC"

        @testset "sequence with PDB file as argument" begin
            prefix = "test/fileContainer/"

            #Case: sequence from a PDB file
            PDBfile = prefix * "test.pdb"
            @test sequence(PDBfile) == "INFYGELVDLGVKEKLIEKAGAWYSYKGEKIGQGKANATAWLKDNPETAKEIEKKVRELLLSN"

            #Error : Input file is supposed to be a PDB file; check if extension is coherent with content
            @test_throws ArgumentError sequence("test.txt")

            #Error : Something went wrong while reading the PDB file
            PDBfile = prefix * "test_empty.pdb"

            ## TODO : does not work if .pdb is empty or does not have any atoms ?
            #@test_throws ArgumentError sequence(PDBfile, 1)

            # Error : does not work if .pdb does not exist
            PDBfile = "doesnotexist.pdb"
            @test_throws ArgumentError sequence(PDBfile)
            end
    end

    @testset "dot" begin
    g = UndirectedGraph{myAtom}()

    # Define sample atoms
    atom1 = myAtom("N", 1, "ALA")
    atom2 = myAtom("CA", 1, "ALA")
    atom3 = myAtom("C", 1, "ALA")

    # Add atoms and edges to the graph
    add!(g, atom1, atom2)
    add!(g, atom2, atom3)
    add!(g, atom1, atom3)

    # Define the expected dot string
    expected_dot = "graph {\n  \"ALA1-N\" -- \"ALA1-C\";\n  \"ALA1-N\" -- \"ALA1-CA\";\n  \"ALA1-CA\" -- \"ALA1-C\";\n}"

    # Test the dot function
    @test dot(g) == expected_dot
    end
end