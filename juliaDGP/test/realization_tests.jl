@testset "Realization type various overrides and functions" begin
    #length function
    @testset "length" begin
        K = 3
        N = 5
        realization = Realization(K, N)
        @test length(realization) == N
        N = 2
        realization = Realization(K, N)
        @test length(realization) == N
    end

    #elements function
    @testset "elements" begin
        sample_atoms = create_sample_PDBTool_atoms()
        realization = Realization(sample_atoms)
        atom_elements = elements(realization)

        @test length(elements(realization)) == 3
        @test length(atom_elements) == length(sample_atoms)

        for atom in sample_atoms
            myatom = myAtom(atom)
            @test myatom in atom_elements
        end
    end

    #override of the show function
    @testset "override of the Base show function" begin
        atoms = create_sample_PDBTool_atoms()
        realization = Realization(atoms)
        io = IOBuffer()
        show(io,realization)
        output = String(take!(io))
        expected_output = "Realization{myAtom} (3 elements in dimension 3)"
        @test output == expected_output
    end

    @testset "Details Function" begin
        atoms = create_sample_PDBTool_atoms()
        realization = Realization(atoms)
        #TODO
        #@test print_to_string(details,realization) == ""
    end
end

@testset "Realization constructors" begin
    # Test the random Realization constructor
    @testset "Random Realization" begin
        K = 3
        N = 5
        realization = Realization(K, N)
        @test realization.K == K
        @test realization.N == N
        @test length(realization) == N
        @test length(elements(realization)) == N
        # Test error cases :
        # Error : The dimension cannot be nonpositive
        @test_throws ArgumentError Realization(-1, 5)
        @test_throws ArgumentError Realization(0, -8)
        # Error : The size of the DGP instance cannot be nonpositive
        @test_throws ArgumentError Realization(3, 0)
    end

    # Test the constructor from a list of Atom instances
    @testset "Realization from Atoms" begin
        atoms = create_sample_PDBTool_atoms()
        realization = Realization(atoms)
        @test realization.K == 3
        @test realization.N == length(atoms)
        @test length(realization) == length(atoms)
        @test length(elements(realization)) == length(atoms)
        # Error : Input Vector{Atom} has zero length
        empty_atoms = Vector{Atom}()
        @test_throws ArgumentError Realization(empty_atoms)
    end

    # Test the constructor from a PDB file
    @testset "Realization from PDB file" begin
        prefix = "test/fileContainer/"
        PDBfile = prefix * "test.pdb"
        modelId = 1
        realization = Realization(PDBfile, modelId)
        @test realization.K == 3
        @test realization.N == 1018
        @test length(realization) == 1018
        @test length(elements(realization)) == 1018
        # Test error cases :
        # Error : Input file is supposed to be a PDB file; check if extension is coherent with content
        @test_throws ArgumentError Realization("test.txt", 1)
        # Error : Model ID cannot be nonpositive
        @test_throws ArgumentError Realization(PDBfile, 0)
        # Error : Something went wrong while reading the PDB file
        PDBfile = prefix * "doesnotexist.pdb"
        ## TODO : does not work if .pdb does not exist
        #@test_throws ArgumentError Realization(PDBfile, 1)
        # Error : Something went wrong while selecting the requested model
        PDBfile = prefix * "test_empty.pdb"
        ## TODO : does not work if .pdb doesn't not have atom in it
        #@test_throws ArgumentError Realization(PDBfile, 1)
    end
end