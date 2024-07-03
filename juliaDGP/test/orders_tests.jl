@testset "vertex_order" begin

    # extracting the vertex order from a Realization instance
    @testset "vertex order from a Realization instance" begin
        atoms = create_sample_PDBTool_atoms()
        realization = Realization(atoms)
        result = vertex_order(realization)
        N = myAtom("N",1,"ALA")
        CA = myAtom("CA",1,"ALA")
        C = myAtom("C",1,"ALA")
        expected = [N,CA,C]
        @test result == expected
    end

    # extracting the vertex order from a DGP instance
    @testset "vertex order from a DGP instance" begin
        atoms = create_sample_PDBTool_atoms()
        realization = Realization(atoms)
        dgp = DGP(realization)
        result = vertex_order(dgp)
        N = myAtom("N",1,"ALA")
        CA = myAtom("CA",1,"ALA")
        C = myAtom("C",1,"ALA")
        expected = [N,CA,C]
        @test result == expected
    end
end

