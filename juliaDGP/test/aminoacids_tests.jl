@testset "codes_bijection" begin
    map = codes_bijection()
    @test map['A'] == "ALA"
    @test map['R'] == "ARG"
    @test map['N'] == "ASN"
    @test map['D'] == "ASP"
    @test map['B'] == "ASX"
    @test map['C'] == "CYS"
    @test map['E'] == "GLU"
    @test map['Q'] == "GLN"
    @test map['Z'] == "GLX"
    @test map['G'] == "GLY"
    @test map['H'] == "HIS"
    @test map['I'] == "ILE"
    @test map['L'] == "LEU"
    @test map['K'] == "LYS"
    @test map['M'] == "MET"
    @test map['F'] == "PHE"
    @test map['P'] == "PRO"
    @test map['S'] == "SER"
    @test map['T'] == "THR"
    @test map['W'] == "TRP"
    @test map['Y'] == "TYR"
    @test map['V'] == "VAL"
    # Test error cases :
    # Error : key not found
    @test_throws KeyError map['@']
end

@testset "fun_dict" begin
    dict = fun_dict()
    @test dict['A'] == "alanine"
    @test dict['R'] == "arginine"
    @test dict['N'] == "asparagine"
    @test dict['D'] == "aspartic"
    @test dict['C'] == "cysteine"
    @test dict['E'] == "glutamic"
    @test dict['Q'] == "glutamine"
    @test dict['G'] == "glycine"
    @test dict['H'] == "histidine"
    @test dict['I'] == "isoleucine"
    @test dict['L'] == "leucine"
    @test dict['K'] == "lysine"
    @test dict['M'] == "methionine"
    @test dict['F'] == "phenylalanine"
    @test dict['P'] == "proline"
    @test dict['S'] == "serine"
    @test dict['T'] == "threonine"
    @test dict['W'] == "tryptophan"
    @test dict['Y'] == "tyrosine"
    @test dict['V'] == "valine"
    # Test error cases :
    # Error : key not found
    @test_throws KeyError dict['@']
end

@testset "aa_name_convert" begin
    # Tests for converting 1-letter to 3-letter codes
    @test aa_name_convert("A") == "ALA"
    @test aa_name_convert("R") == "ARG"
    @test aa_name_convert("N") == "ASN"
    @test aa_name_convert("D") == "ASP"
    @test aa_name_convert("C") == "CYS"
    @test aa_name_convert("E") == "GLU"
    @test aa_name_convert("Q") == "GLN"
    @test aa_name_convert("G") == "GLY"
    @test aa_name_convert("H") == "HIS"
    @test aa_name_convert("I") == "ILE"
    @test aa_name_convert("L") == "LEU"
    @test aa_name_convert("K") == "LYS"
    @test aa_name_convert("M") == "MET"
    @test aa_name_convert("F") == "PHE"
    @test aa_name_convert("P") == "PRO"
    @test aa_name_convert("S") == "SER"
    @test aa_name_convert("T") == "THR"
    @test aa_name_convert("W") == "TRP"
    @test aa_name_convert("Y") == "TYR"
    @test aa_name_convert("V") == "VAL"
    # Tests for converting 3-letter to 1-letter codes
    @test aa_name_convert("ALA") == 'A'
    @test aa_name_convert("ARG") == 'R'
    @test aa_name_convert("ASN") == 'N'
    @test aa_name_convert("ASP") == 'D'
    @test aa_name_convert("CYS") == 'C'
    @test aa_name_convert("GLU") == 'E'
    @test aa_name_convert("GLN") == 'Q'
    @test aa_name_convert("GLY") == 'G'
    @test aa_name_convert("HIS") == 'H'
    @test aa_name_convert("ILE") == 'I'
    @test aa_name_convert("LEU") == 'L'
    @test aa_name_convert("LYS") == 'K'
    @test aa_name_convert("MET") == 'M'
    @test aa_name_convert("PHE") == 'F'
    @test aa_name_convert("PRO") == 'P'
    @test aa_name_convert("SER") == 'S'
    @test aa_name_convert("THR") == 'T'
    @test aa_name_convert("TRP") == 'W'
    @test aa_name_convert("TYR") == 'Y'
    @test aa_name_convert("VAL") == 'V'
    # Test error cases :
    # Error : Unknown 1-letter code for amino acid
    @test_throws ArgumentError aa_name_convert("@")
    # Error : Unknown 3-letter code for amino acid
    @test_throws ArgumentError aa_name_convert("@.A")
    # Error : Amino acid codes can be composed either by 1 or 3 letters
    @test_throws ArgumentError aa_name_convert("XT")
end

@testset "backbone! and backbone" begin
    g = UndirectedGraph{myAtom}()
    resnum = 1
    resname = "ALA"
    backbone!(g, resnum, resname)

    # Create expected atoms
    N = myAtom("N",resnum, resname)
    H = myAtom("H",resnum, resname)
    CA = myAtom("CA",resnum, resname)
    C = myAtom("C",resnum, resname)
    O = myAtom("O",resnum, resname)

    # Create expected bonds
    expected_edges = Set([
    (N,H),
    (N,CA),
    (CA,C),
    (C,O)
    ])

    #TODO
end