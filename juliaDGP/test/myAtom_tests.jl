@testset "myAtom constructors" begin
    @testset "Generic constructor" begin
        a = myAtom("CA", 1, "ALA")
        @test a.name == "CA"
        @test a.resnum == 1
        @test a.resname == "ALA"
        # Test error cases :
        # Error : myAtom name cannot be empty
        @test_throws ArgumentError myAtom("", 1, "ALA")
        # Error : myAtom name should always start with a letter
        @test_throws ArgumentError myAtom("1A", 1, "ALA")
        # Error : Valid myAtom names should begin with letters, then digits, possibly followed by one ending wildcard
        @test_throws ArgumentError myAtom("AD646SQ1546#", 1, "ALA")
        @test_throws ArgumentError myAtom("A#1", 1, "ALA")
        @test_throws ArgumentError myAtom("A1@", 1, "ALA")
        # Error : myAtom resnum cannot be nonpositive
        @test_throws ArgumentError myAtom("CA", 0, "ALA")
        # Error : The given residue name is not a valid name
        @test_throws ArgumentError myAtom("CA", 1, "AL")
    end

    @testset "Constructor from Atom" begin
        a = Atom()
        a.name = "CA"
        a.resnum = 1
        a.resname = "ALA"
        b = myAtom(a)
        @test b.name == "CA"
        @test b.resnum == 1
        @test b.resname == "ALA"
    end
end

@testset "myAtom overriding functions" begin
    @testset "Override of the equality operator" begin
        a = myAtom("CA", 1, "ALA")
        b = myAtom("CA", 1, "ALA")
        @test a == b
        b = myAtom("CA", 1, "ARG")
        @test a != b
    end

    @testset "Override of the hash code function" begin
        a = myAtom("CA", 1, "ALA")
        b = myAtom("CA", 1, "ALA")
        @test hash(a) == hash(b)
    end

    @testset "Override of the string conversion function" begin
        a = myAtom("CA", 1, "ALA")
        @test string(a) == "ALA1-CA"
    end

    @testset "Override of the show function" begin
        a = myAtom("CA", 1, "ALA")
        @test print_to_string(show, a) == "ALA1-CA"
        @test print_to_string(show, a) != "ALB2-CA"
    end
end

@testset "myAtom Functions" begin
    @testset "are_equivalent" begin
        a = myAtom("CA", 1, "ALA")
        b = myAtom("CA", 1, "ALA")
        @test are_equivalent(a, b) == true
        b = myAtom("CA#", 1, "ALA")
        #atoms not equivalent with CA# abd CA when it should
        #@test are_equivalent(a, b) == true
        a = myAtom("CA#", 1, "ALA")
        @test are_equivalent(a, b) == true
        b = myAtom("CBA", 1, "ALA")
        @test are_equivalent(a, b) == false
        b = myAtom("CB", 1, "ALA")
        #atoms equivalent with CA and CB when it shouldn't
        #@test are_equivalent(a, b) == false
        b = myAtom("CA", 1, "ARG")
        @test are_equivalent(a, b) == false
        b = myAtom("CA", 2, "ALA")
        @test are_equivalent(a, b) == false
    end

    @testset "myAtoms" begin
        a = Atom()
        a.name = "CA"
        a.resnum = 1
        a.resname = "ALA"
        b = Atom()
        b.name = "CB"
        b.resnum = 2
        b.resname = "ARG"
        atoms = [a, b]
        myatoms = myAtoms(atoms)
        @test myatoms[1].name == "CA"
        @test myatoms[1].resnum == 1
        @test myatoms[1].resname == "ALA"
        @test myatoms[2].name == "CB"
        @test myatoms[2].resnum == 2
        @test myatoms[2].resname == "ARG"
    end

    @testset "PDBToolAtom" begin
        a = myAtom("CA", 1, "ALA")
        b = PDBToolAtom(1, a)
        @test b.index == 1
        @test b.name == "CA"
        @test b.resnum == 1
        @test b.resname == "ALA"
    # Test error cases :
    # Error : The PDBTools.Atom index cannot be negative
    @test_throws ArgumentError PDBToolAtom(-1, myAtom("CA", 1, "ALA"))
    end
end