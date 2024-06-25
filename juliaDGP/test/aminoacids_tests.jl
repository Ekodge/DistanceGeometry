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

@testset "backbone and backbone!" begin
    resnum = rand(1:1000)
    resname = "ALA"
    g = backbone(resnum, resname)

    # Create expected atoms corresponding to backgone
    N = myAtom("N",resnum, resname)
    H = myAtom("H",resnum, resname)
    CA = myAtom("CA",resnum, resname)
    C = myAtom("C",resnum, resname)
    O = myAtom("O",resnum, resname)

    expected_atoms_backbone = Set([N,H,CA,C,O])

    # Create expected bonds corresponding to backbone
    expected_edges_backbone = Set([
    (N,H),
    (N,CA),
    (CA,C),
    (C,O)
    ])

    @test counter(vlist(g)) == counter(expected_atoms_backbone)
    @test sets_of_tuples_equal(Set(elist(g)), expected_edges_backbone)
end

@testset "alanine and alanine!" begin
    resnum = rand(1:1000)
    resname = "ALA"
    g = alanine(resnum)

    # Create expected atoms corresponding to backgone
    N = myAtom("N",resnum, resname)
    H = myAtom("H",resnum, resname)
    CA = myAtom("CA",resnum, resname)
    C = myAtom("C",resnum, resname)
    O = myAtom("O",resnum, resname)

    expected_atoms_backbone = Set([N,H,CA,C,O])

    # Create expected bonds corresponding to backbone
    expected_edges_backbone = Set([
    (N,H),
    (N,CA),
    (CA,C),
    (C,O)
    ])

    # Create expected atoms
    CA = myAtom("CA",resnum,resname);
    HA = myAtom("HA",resnum,resname);
    CB = myAtom("CB",resnum,resname);
    HB = myAtom("HB#",resnum,resname);

    expected_atoms = union!(expected_atoms_backbone,Set([CA,HA,CB,HB]))

    # Create expected bonds
    expected_edges = union!(expected_edges_backbone,Set([
    (CA,HA),
    (CA,CB),
    (CB,HB)
    ]))

    @test counter(vlist(g)) == counter(expected_atoms)
    @test sets_of_tuples_equal(Set(elist(g)), expected_edges)
end

@testset "arginine and arginine!" begin
    resnum = rand(1:1000)
    resname = "ARG"
    g = arginine(resnum)
    backbone = backbone_atoms_bonds(resnum,resname)
    expected_atoms_backbone = backbone.atoms
    expected_edges_backbone = backbone.bonds

    # Create expected atoms
    CA = myAtom("CA",resnum,resname);
    HA = myAtom("HA",resnum,resname);
    CB = myAtom("CB",resnum,resname);
    CG = myAtom("CG",resnum,resname);
    CD = myAtom("CD",resnum,resname);
    NE = myAtom("NE",resnum,resname);
    CZ = myAtom("CZ",resnum,resname);
    NH1 = myAtom("NH1",resnum,resname);
    NH2 = myAtom("NH2",resnum,resname);
    HB = myAtom("HB#",resnum,resname);
    HG = myAtom("HG#",resnum,resname);
    HD = myAtom("HD#",resnum,resname);
    HE = myAtom("HE",resnum,resname);
    HH1 = myAtom("HH1#",resnum,resname);
    HH2 = myAtom("HH2#",resnum,resname);

    expected_atoms = union!(expected_atoms_backbone,Set([CA,HA,CB,CG,CD,NE,CZ,NH1,NH2,HB,HG,HD,HE,HH1,HH2]))

    # Create expected bonds
    expected_edges = union!(expected_edges_backbone,Set([
    (CA,HA),
    (CA,CB),
    (CB,CG),
    (CG,CD),
    (CD,NE),
    (NE,CZ),
    (CZ,NH1),
    (CZ,NH2),
    (CB,HB),
    (CG,HG),
    (CD,HD),
    (NE,HE),
    (NH1,HH1),
    (NH2,HH2)
    ]))

    @test counter(vlist(g)) == counter(expected_atoms)
    @test sets_of_tuples_equal(Set(elist(g)), expected_edges)
end

@testset "asparagine and asparagine!" begin
    resnum = rand(1:1000)
    resname = "ASN"
    g = asparagine(resnum)
    backbone = backbone_atoms_bonds(resnum,resname)
    expected_atoms_backbone = backbone.atoms
    expected_edges_backbone = backbone.bonds

    # Create expected atoms
    CA = myAtom("CA",resnum,resname);
    HA = myAtom("HA",resnum,resname);
    CB = myAtom("CB",resnum,resname);
    CG = myAtom("CG",resnum,resname);
    OD1 = myAtom("OD1",resnum,resname);
    ND2 = myAtom("ND2",resnum,resname);
    HB = myAtom("HB#",resnum,resname);
    HD2 = myAtom("HD2#",resnum,resname);

    expected_atoms = union!(expected_atoms_backbone,Set([CA,HA,CB,CG,OD1,ND2,HB,HD2]))

    # Create expected bonds
    expected_edges = union!(expected_edges_backbone,Set([
    (CA,HA),
    (CA,CB),
    (CB,CG),
    (CG,OD1),
    (CG,ND2),
    (CB,HB),
    (ND2,HD2)
    ]))

    @test counter(vlist(g)) == counter(expected_atoms)
    @test sets_of_tuples_equal(Set(elist(g)), expected_edges)
end

@testset "aspartic and aspartic!" begin
    resnum = rand(1:1000)
    resname = "ASP"
    g = aspartic(resnum)
    backbone = backbone_atoms_bonds(resnum,resname)
    expected_atoms_backbone = backbone.atoms
    expected_edges_backbone = backbone.bonds

    # Create expected atoms
    CA = myAtom("CA",resnum,resname);
    HA = myAtom("HA",resnum,resname);
    CB = myAtom("CB",resnum,resname);
    CG = myAtom("CG",resnum,resname);
    OD = myAtom("OD#",resnum,resname);
    HB = myAtom("HB",resnum,resname);

    expected_atoms = union!(expected_atoms_backbone,Set([CA,HA,CB,CG,OD,HB]))

    # Create expected bonds
    expected_edges = union!(expected_edges_backbone,Set([
    (CA,HA),
    (CA,CB),
    (CB,CG),
    (CG,OD),
    (CB,HB)
    ]))

    @test counter(vlist(g)) == counter(expected_atoms)
    @test sets_of_tuples_equal(Set(elist(g)), expected_edges)
end

@testset "cysteine and cysteine!" begin
    resnum = rand(1:1000)
    resname = "CYS"
    g = cysteine(resnum)
    backbone = backbone_atoms_bonds(resnum,resname)
    expected_atoms_backbone = backbone.atoms
    expected_edges_backbone = backbone.bonds

    # Create expected atoms
    CA = myAtom("CA",resnum,resname);
    HA = myAtom("HA",resnum,resname);
    CB = myAtom("CB",resnum,resname);
    SG = myAtom("SG",resnum,resname);
    HB = myAtom("HB#",resnum,resname);
    HG = myAtom("HG",resnum,resname);

    expected_atoms = union!(expected_atoms_backbone,Set([CA,HA,CB,SG,HB,HG]))

    # Create expected bonds
    expected_edges = union!(expected_edges_backbone,Set([
    (CA,HA),
    (CA,CB),
    (CB,SG),
    (CB,HB),
    (SG,HG)
    ]))

    @test counter(vlist(g)) == counter(expected_atoms)
    @test sets_of_tuples_equal(Set(elist(g)), expected_edges)
end

@testset "glutamine and glutamine!" begin
    resnum = rand(1:1000)
    resname = "GLN"
    g = glutamine(resnum)
    backbone = backbone_atoms_bonds(resnum,resname)
    expected_atoms_backbone = backbone.atoms
    expected_edges_backbone = backbone.bonds

    # Create expected atoms
    CA = myAtom("CA",resnum,resname);
    HA = myAtom("HA",resnum,resname);
    CB = myAtom("CB",resnum,resname);
    CG = myAtom("CG",resnum,resname);
    CD = myAtom("CD",resnum,resname);
    OE1 = myAtom("OE1",resnum,resname);
    NE2 = myAtom("NE2",resnum,resname);
    HB = myAtom("HB#",resnum,resname);
    HG = myAtom("HG#",resnum,resname);
    HE2 = myAtom("HE2",resnum,resname);

    expected_atoms = union!(expected_atoms_backbone,Set([CA,HA,CB,CG,CD,OE1,NE2,HB,HG,HE2]))

    # Create expected bonds
    expected_edges = union!(expected_edges_backbone,Set([
    (CA,HA),
    (CA,CB),
    (CB,CG),
    (CG,CD),
    (CD,OE1),
    (CD,NE2),
    (CB,HB),
    (CG,HG),
    (NE2,HE2),
    ]))

    @test counter(vlist(g)) == counter(expected_atoms)
    @test sets_of_tuples_equal(Set(elist(g)), expected_edges)
end

@testset "glutamic and glutamic!" begin
    resnum = rand(1:1000)
    resname = "GLU"
    g = glutamic(resnum)
    backbone = backbone_atoms_bonds(resnum,resname)
    expected_atoms_backbone = backbone.atoms
    expected_edges_backbone = backbone.bonds

    # Create expected atoms
    CA = myAtom("CA",resnum,resname);
    HA = myAtom("HA",resnum,resname);
    CB = myAtom("CB",resnum,resname);
    CG = myAtom("CG",resnum,resname);
    CD = myAtom("CD",resnum,resname);
    OE = myAtom("OE",resnum,resname);
    HB = myAtom("HB#",resnum,resname);
    HG = myAtom("HG#",resnum,resname);

    expected_atoms = union!(expected_atoms_backbone,Set([CA,HA,CB,CG,CD,OE,HB,HG]))

    # Create expected bonds
    expected_edges = union!(expected_edges_backbone,Set([
    (CA,HA),
    (CA,CB),
    (CB,CG),
    (CG,CD),
    (CD,OE),
    (CB,HB),
    (CG,HG)
    ]))

    @test counter(vlist(g)) == counter(expected_atoms)
    @test sets_of_tuples_equal(Set(elist(g)), expected_edges)
end

@testset "glycine and glycine!" begin
    resnum = rand(1:1000)
    resname = "GLY"
    g = glycine(resnum)
    backbone = backbone_atoms_bonds(resnum,resname)
    expected_atoms_backbone = backbone.atoms
    expected_edges_backbone = backbone.bonds

    # Create expected atoms
    CA = myAtom("CA",resnum,resname);
    HA = myAtom("HA#",resnum,resname);

    expected_atoms = union!(expected_atoms_backbone,Set([CA,HA]))

    # Create expected bonds
    expected_edges = union!(expected_edges_backbone,Set([
    (CA,HA)
    ]))

    @test counter(vlist(g)) == counter(expected_atoms)
    @test sets_of_tuples_equal(Set(elist(g)), expected_edges)
end

@testset "histidine and histidine!" begin
    resnum = rand(1:1000)
    resname = "HIS"
    g = histidine(resnum)
    backbone = backbone_atoms_bonds(resnum,resname)
    expected_atoms_backbone = backbone.atoms
    expected_edges_backbone = backbone.bonds

    # Create expected atoms
    CA = myAtom("CA",resnum,resname);
    HA = myAtom("HA",resnum,resname);
    CB = myAtom("CB",resnum,resname);
    CG = myAtom("CG",resnum,resname);
    ND1 = myAtom("ND1",resnum,resname);
    CE1 = myAtom("CE1",resnum,resname);
    NE2 = myAtom("NE2",resnum,resname);
    CD2 = myAtom("CD2",resnum,resname);
    HB = myAtom("HB#",resnum,resname);
    HD1 = myAtom("HD1",resnum,resname);
    HE1 = myAtom("HE1",resnum,resname);
    HE2 = myAtom("HE2",resnum,resname);
    HD2 = myAtom("HD2",resnum,resname);

    expected_atoms = union!(expected_atoms_backbone,Set([CA,HA,CB,CG,ND1,CE1,NE2,CD2,HB,HD1,HE1,HE2,HD2]))

    # Create expected bonds
    expected_edges = union!(expected_edges_backbone,Set([
    (CA,HA),
    (CA,CB),
    (CB,CG),
    (CG,ND1),
    (ND1,CE1),
    (CE1,NE2),
    (NE2,CD2),
    (CD2,CG),
    (CB,HB),
    (ND1,HD1),
    (CE1,HE1),
    (NE2,HE2),
    (CD2,HD2)
    ]))

    @test counter(vlist(g)) == counter(expected_atoms)
    @test sets_of_tuples_equal(Set(elist(g)), expected_edges)
end

@testset "isoleucine and isoleucine!" begin
    resnum = rand(1:1000)
    resname = "ILE"
    g = isoleucine(resnum)
    backbone = backbone_atoms_bonds(resnum,resname)
    expected_atoms_backbone = backbone.atoms
    expected_edges_backbone = backbone.bonds

    # Create expected atoms
    CA = myAtom("CA",resnum,resname);
    HA = myAtom("HA",resnum,resname);
    CB = myAtom("CB",resnum,resname);
    HB = myAtom("HB#",resnum,resname);
    CG1 = myAtom("CG1",resnum,resname);
    CG2 = myAtom("CG2",resnum,resname);
    CD1 = myAtom("CD1",resnum,resname);
    HG1 = myAtom("HG1#",resnum,resname);
    HG2 = myAtom("HG2#",resnum,resname);
    HD1 = myAtom("HD1#",resnum,resname);

    expected_atoms = union!(expected_atoms_backbone, Set([CA, HA, CB, HB, CG1, CG2, CD1, HG1, HG2, HD1]))

    # Create expected bonds
    expected_edges = union!(expected_edges_backbone, Set([
    (CA, HA),
    (CA, CB),
    (CB, CG1),
    (CB, CG2),
    (CG1, CD1),
    (CB, HB),
    (CG1, HG1),
    (CG2, HG2),
    (CD1, HD1)
    ]))

    @test counter(vlist(g)) == counter(expected_atoms)
    @test sets_of_tuples_equal(Set(elist(g)), expected_edges)
end

@testset "leucine and leucine!" begin
    resnum = rand(1:1000)
    resname = "LEU"
    g = leucine(resnum)
    backbone = backbone_atoms_bonds(resnum,resname)
    expected_atoms_backbone = backbone.atoms
    expected_edges_backbone = backbone.bonds

    # Create expected atoms
    CA = myAtom("CA",resnum,resname);
    HA = myAtom("HA",resnum,resname);
    CB = myAtom("CB",resnum,resname);
    CG = myAtom("CG",resnum,resname);
    CD1 = myAtom("CD1",resnum,resname);
    CD2 = myAtom("CD2",resnum,resname);
    HB = myAtom("HB#",resnum,resname);
    HG = myAtom("HG",resnum,resname);
    HD1 = myAtom("HD1#",resnum,resname);
    HD2 = myAtom("HD2#",resnum,resname);

    expected_atoms = union!(expected_atoms_backbone, Set([CA, HA, CB, CG, CD1, CD2, HB, HG, HD1, HD2]))

    # Create expected bonds
    expected_edges = union!(expected_edges_backbone, Set([
    (CA, HA),
    (CA, CB),
    (CB, CG),
    (CG, CD1),
    (CG, CD2),
    (CB, HB),
    (CG, HG),
    (CD1, HD1),
    (CD2, HD2)
    ]))

    @test counter(vlist(g)) == counter(expected_atoms)
    @test sets_of_tuples_equal(Set(elist(g)), expected_edges)
end

@testset "lysine and lysine!" begin
    resnum = rand(1:1000)
    resname = "LYS"
    g = lysine(resnum)
    backbone = backbone_atoms_bonds(resnum,resname)
    expected_atoms_backbone = backbone.atoms
    expected_edges_backbone = backbone.bonds

    # Create expected atoms
    CA = myAtom("CA",resnum,resname);
    HA = myAtom("HA",resnum,resname);
    CB = myAtom("CB",resnum,resname);
    CG = myAtom("CG",resnum,resname);
    CD = myAtom("CD",resnum,resname);
    CE = myAtom("CE",resnum,resname);
    NZ = myAtom("NZ",resnum,resname);
    HB = myAtom("HB#",resnum,resname);
    HG = myAtom("HG#",resnum,resname);
    HD = myAtom("HD#",resnum,resname);
    HE = myAtom("HE#",resnum,resname);
    HZ = myAtom("HZ#",resnum,resname);

    expected_atoms = union!(expected_atoms_backbone, Set([CA, HA, CB, CG, CD, CE, NZ, HB, HG, HD, HE, HZ]))

    # Create expected bonds
    expected_edges = union!(expected_edges_backbone, Set([
    (CA, HA),
    (CA, CB),
    (CB, CG),
    (CG, CD),
    (CD, CE),
    (CE, NZ),
    (CB, HB),
    (CG, HG),
    (CD, HD),
    (CE, HE),
    (NZ, HZ)
    ]))

    @test counter(vlist(g)) == counter(expected_atoms)
    @test sets_of_tuples_equal(Set(elist(g)), expected_edges)
end

@testset "methionine and methionine!" begin
    resnum = rand(1:1000)
    resname = "MET"
    g = methionine(resnum)

    # Create expected backbone atoms
    N = myAtom("N",resnum,resname)
    H = myAtom("H#",resnum,resname)
    CA = myAtom("CA",resnum,resname)
    HA = myAtom("HA",resnum,resname)
    C = myAtom("C",resnum,resname)
    O = myAtom("O",resnum,resname)

    expected_atoms_backbone = Set([N,H,CA,HA,C,O])

    # Create expected atoms
    CB = myAtom("CB",resnum,resname);
    HB = myAtom("HB#",resnum,resname);
    CG = myAtom("CG",resnum,resname);
    SD = myAtom("SD",resnum,resname);
    CE = myAtom("CE",resnum,resname);
    HG = myAtom("HG#",resnum,resname);
    HE = myAtom("HE#",resnum,resname);

    expected_atoms = union!(expected_atoms_backbone, Set([CB, HB, CG, SD, CE, HG, HE]))

    expected_edges_backbone = Set([
    (N,H),
    (N,CA),
    (CA,C),
    (C,O)
    ])

    # Create expected bonds
    expected_edges = union!(expected_edges_backbone, Set([
    (CA, HA),
    (CA, CB),
    (CB, CG),
    (CG, SD),
    (SD, CE),
    (CB, HB),
    (CG, HG),
    (CE, HE)
    ]))

    @test counter(vlist(g)) == counter(expected_atoms)
    @test sets_of_tuples_equal(Set(elist(g)), expected_edges)
end

@testset "phenylalanine and phenylalanine!" begin
    resnum = rand(1:1000)
    resname = "PHE"
    g = phenylalanine(resnum)
    backbone = backbone_atoms_bonds(resnum,resname)
    expected_atoms_backbone = backbone.atoms
    expected_edges_backbone = backbone.bonds

    # Create expected atoms
    CA = myAtom("CA",resnum,resname);
    HA = myAtom("HA",resnum,resname);
    CB = myAtom("CB",resnum,resname);
    CG = myAtom("CG",resnum,resname);
    CD1 = myAtom("CD1",resnum,resname);
    CD2 = myAtom("CD2",resnum,resname);
    CE1 = myAtom("CE1",resnum,resname);
    CE2 = myAtom("CE2",resnum,resname);
    CZ = myAtom("CZ",resnum,resname);
    HB = myAtom("HB#",resnum,resname);
    HD1 = myAtom("HD1",resnum,resname);
    HD2 = myAtom("HD2",resnum,resname);
    HE1 = myAtom("HE1",resnum,resname);
    HE2 = myAtom("HE2",resnum,resname);
    HZ = myAtom("HZ",resnum,resname);

    expected_atoms = union!(expected_atoms_backbone, Set([CA, HA, CB, CG, CD1, CD2, CE1, CE2, CZ, HB, HD1, HD2, HE1, HE2, HZ]))

    # Create expected bonds
    expected_edges = union!(expected_edges_backbone, Set([
    (CA, HA),
    (CA, CB),
    (CB, CG),
    (CG, CD1),
    (CG, CD2),
    (CD1, CE1),
    (CD2, CE2),
    (CE1, CZ),
    (CE2, CZ),
    (CB, HB),
    (CD1, HD1),
    (CD2, HD2),
    (CE1, HE1),
    (CE2, HE2),
    (CZ, HZ)
    ]))

    @test counter(vlist(g)) == counter(expected_atoms)
    @test sets_of_tuples_equal(Set(elist(g)), expected_edges)
end

@testset "proline and proline!" begin
    resnum = rand(1:1000)
    resname = "PRO"
    g = proline(resnum)
    backbone = backbone_atoms_bonds(resnum,resname)
    expected_atoms_backbone = backbone.atoms
    expected_edges_backbone = backbone.bonds

    # Create expected atoms
    N = myAtom("N",resnum,resname);
    CA = myAtom("CA",resnum,resname);
    HA = myAtom("HA",resnum,resname);
    CB = myAtom("CB",resnum,resname);
    HB = myAtom("HB#",resnum,resname);
    CG = myAtom("CG",resnum,resname);
    HG = myAtom("HG#",resnum,resname);
    CD = myAtom("CD",resnum,resname);
    HD = myAtom("HD",resnum,resname);

    expected_atoms = union!(expected_atoms_backbone, Set([N, CA, HA, CB, HB, CG, HG, CD, HD]))

    # Create expected bonds
    expected_edges = union!(expected_edges_backbone, Set([
    (CA, HA),
    (CA, CB),
    (CB, CG),
    (CG, CD),
    (CD, N),
    (CB, HB),
    (CG, HG),
    (CD, HD)
    ]))

    @test counter(vlist(g)) == counter(expected_atoms)
    @test sets_of_tuples_equal(Set(elist(g)), expected_edges)
end

@testset "serine and serine!" begin
    resnum = rand(1:1000)
    resname = "SER"
    g = serine(resnum)
    backbone = backbone_atoms_bonds(resnum,resname)
    expected_atoms_backbone = backbone.atoms
    expected_edges_backbone = backbone.bonds

    # Create expected atoms
    CA = myAtom("CA",resnum,resname);
    HA = myAtom("HA",resnum,resname);
    CB = myAtom("CB",resnum,resname);
    HB = myAtom("HB",resnum,resname);
    OG = myAtom("OG",resnum,resname);
    HG = myAtom("HG",resnum,resname);

    expected_atoms = union!(expected_atoms_backbone, Set([CA, HA, CB, HB, OG, HG]))

    # Create expected bonds
    expected_edges = union!(expected_edges_backbone, Set([
    (CA, HA),
    (CA, CB),
    (CB, OG),
    (CB, HB),
    (OG, HG)
    ]))

    @test counter(vlist(g)) == counter(expected_atoms)
    @test sets_of_tuples_equal(Set(elist(g)), expected_edges)
end

@testset "threonine and threonine!" begin
    resnum = rand(1:1000)
    resname = "THR"
    g = threonine(resnum)
    backbone = backbone_atoms_bonds(resnum,resname)
    expected_atoms_backbone = backbone.atoms
    expected_edges_backbone = backbone.bonds

    # Create expected atoms
    CA = myAtom("CA",resnum,resname);
    HA = myAtom("HA",resnum,resname);
    CB = myAtom("CB",resnum,resname);
    HB = myAtom("HB",resnum,resname);
    OG1 = myAtom("OG1",resnum,resname);
    CG2 = myAtom("CG2",resnum,resname);
    HG1 = myAtom("HG1",resnum,resname);
    HG2 = myAtom("HG2#",resnum,resname);

    expected_atoms = union!(expected_atoms_backbone, Set([CA, HA, CB, HB, OG1, CG2, HG1, HG2]))

    # Create expected bonds
    expected_edges = union!(expected_edges_backbone, Set([
    (CA, HA),
    (CA, CB),
    (CB, OG1),
    (CB, CG2),
    (CB, HB),
    (OG1, HG1),
    (CG2, HG2)
    ]))

    @test counter(vlist(g)) == counter(expected_atoms)
    @test sets_of_tuples_equal(Set(elist(g)), expected_edges)
end

@testset "tryptophan and tryptophan!" begin
    resnum = rand(1:1000)
    resname = "TRP"
    g = tryptophan(resnum)
    backbone = backbone_atoms_bonds(resnum,resname)
    expected_atoms_backbone = backbone.atoms
    expected_edges_backbone = backbone.bonds

    # Create expected atoms
    CA = myAtom("CA",resnum,resname);
    HA = myAtom("HA",resnum,resname);
    CB = myAtom("CB",resnum,resname);
    HB = myAtom("HB",resnum,resname);
    CG = myAtom("CG",resnum,resname);
    CD1 = myAtom("CD1",resnum,resname);
    CD2 = myAtom("CD2",resnum,resname);
    NE1 = myAtom("NE1",resnum,resname);
    CE2 = myAtom("CE2",resnum,resname);
    CE3 = myAtom("CE3",resnum,resname);
    CZ2 = myAtom("CZ2",resnum,resname);
    CZ3 = myAtom("CZ3",resnum,resname);
    CH2 = myAtom("CH2",resnum,resname);
    HD1 = myAtom("HD1",resnum,resname);
    HE1 = myAtom("HE1",resnum,resname);
    HE3 = myAtom("HE3",resnum,resname);
    HZ2 = myAtom("HZ2",resnum,resname);
    HZ3 = myAtom("HZ3",resnum,resname);
    HH2 = myAtom("HH2",resnum,resname);

    expected_atoms = union!(expected_atoms_backbone, Set([CA, HA, CB, HB, CG, CD1, CD2, NE1, CE2, CE3, CZ2, CZ3, CH2, HD1, HE1, HE3, HZ2, HZ3, HH2]))
    # Create expected bonds
    expected_edges = union!(expected_edges_backbone, Set([
    (CA, HA),
    (CA, CB),
    (CB, CG),
    (CG, CD1),
    (CG, CD2),
    (CD1, NE1),
    (CD2, CE2),
    (NE1, CE2),
    (CD2, CE3),
    (CE2, CZ2),
    (CE3, CZ3),
    (CZ3, CH2),
    (CH2, CZ2),
    (CB, HB),
    (CD1, HD1),
    (NE1, HE1),
    (CE3, HE3),
    (CZ3, HZ3),
    (CH2, HH2),
    (CZ2, HZ2)
    ]))

    @test counter(vlist(g)) == counter(expected_atoms)
    @test sets_of_tuples_equal(Set(elist(g)), expected_edges)
end

@testset "tyrosine and tyrosine!" begin
    resnum = rand(1:1000)
    resname = "TYR"
    g = tyrosine(resnum)
    backbone = backbone_atoms_bonds(resnum,resname)
    expected_atoms_backbone = backbone.atoms
    expected_edges_backbone = backbone.bonds

    # Create expected atoms
    CA = myAtom("CA",resnum,resname);
    HA = myAtom("HA",resnum,resname);
    CB = myAtom("CB",resnum,resname);
    HB = myAtom("HB",resnum,resname);
    CG = myAtom("CG",resnum,resname);
    CD1 = myAtom("CD1",resnum,resname);
    CD2 = myAtom("CD2",resnum,resname);
    CE1 = myAtom("CE1",resnum,resname);
    CE2 = myAtom("CE2",resnum,resname);
    CZ = myAtom("CZ",resnum,resname);
    OH = myAtom("OH",resnum,resname);
    HD1 = myAtom("HD1",resnum,resname);
    HD2 = myAtom("HD2",resnum,resname);
    HE1 = myAtom("HE1",resnum,resname);
    HE2 = myAtom("HE2",resnum,resname);
    HH = myAtom("HH",resnum,resname);

    expected_atoms = union!(expected_atoms_backbone, Set([CA, HA, CB, HB, CG, CD1, CD2, CE1, CE2, CZ, OH, HD1, HD2, HE1, HE2, HH]))

    # Create expected bonds
    expected_edges = union!(expected_edges_backbone, Set([
    (CA, HA),
    (CA, CB),
    (CB, CG),
    (CG, CD1),
    (CG, CD2),
    (CD1, CE1),
    (CD2, CE2),
    (CE1, CZ),
    (CE2, CZ),
    (CZ, OH),
    (CB, HB),
    (CD1, HD1),
    (CD2, HD2),
    (CE1, HE1),
    (CE2, HE2),
    (OH, HH)
    ]))

    @test counter(vlist(g)) == counter(expected_atoms)
    @test sets_of_tuples_equal(Set(elist(g)), expected_edges)
end

@testset "valine and valine!" begin
    resnum = rand(1:1000)
    resname = "VAL"
    g = valine(resnum)
    backbone = backbone_atoms_bonds(resnum,resname)
    expected_atoms_backbone = backbone.atoms
    expected_edges_backbone = backbone.bonds

    # Create expected atoms
    CA = myAtom("CA",resnum,resname);
    HA = myAtom("HA",resnum,resname);
    CB = myAtom("CB",resnum,resname);
    HB = myAtom("HB",resnum,resname);
    CG1 = myAtom("CG1",resnum,resname);
    CG2 = myAtom("CG2",resnum,resname);
    HG1 = myAtom("HG1",resnum,resname);
    HG2 = myAtom("HG2",resnum,resname);

    expected_atoms = union!(expected_atoms_backbone, Set([CA, HA, CB, HB, CG1, CG2, HG1, HG2]))

    # Create expected bonds
    expected_edges = union!(expected_edges_backbone, Set([
    (CA, HA),
    (CA, CB),
    (CB, CG1),
    (CB, CG2),
    (CB, HB),
    (CG1, HG1),
    (CG2, HG2)
    ]))

    @test counter(vlist(g)) == counter(expected_atoms)
    @test sets_of_tuples_equal(Set(elist(g)), expected_edges)
end

@testset "Protein" begin
    @testset "Empty sequence" begin
        @test_throws ArgumentError protein("")
    end

    @testset "Single amino acid" begin
        sequence = "G"  # Glycine
        g = protein(sequence)
        expected_atoms, expected_bonds = expected_protein_graph(sequence)
        @test counter(vlist(g)) == counter(expected_atoms)
        @test sets_of_tuples_equal(Set(elist(g)), expected_bonds)
    end

    @testset "Multiple amino acids" begin
        sequences = ["GA", "GV", "GAV", "GAAGVV", "G"^10]  # Different sequences
        for sequence in sequences
            g = protein(sequence)
            expected_atoms, expected_bonds = expected_protein_graph(sequence)
            @test counter(vlist(g)) == counter(expected_atoms)
            @test sets_of_tuples_equal(Set(elist(g)), expected_bonds)
        end
    end

    @testset "Complex sequences" begin
        sequences = ["MTEYKLVVVG", "ADLQGKSSS", "MKTFARKVV"]  # More complex sequences
        for sequence in sequences
            g = protein(sequence)
            expected_atoms, expected_bonds = expected_protein_graph(sequence)
            @test counter(vlist(g)) == counter(expected_atoms)
            @test sets_of_tuples_equal(Set(elist(g)), expected_bonds)
        end
    end
end

@testset "enrich!" begin
    resnum = rand(1:1000)
    resname = "GLU"
    g = glutamic(resnum)
    backbone = backbone_atoms_bonds(resnum,resname)
    expected_atoms_backbone = backbone.atoms
    expected_edges_backbone = backbone.bonds

    # Create expected backbone atoms
    N = myAtom("N",resnum,resname)
    H = myAtom("H",resnum,resname)
    CA = myAtom("CA",resnum,resname)
    HA = myAtom("HA",resnum,resname)
    C = myAtom("C",resnum,resname)
    O = myAtom("O",resnum,resname)

    # Create expected atoms
    HA = myAtom("HA",resnum,resname);
    CB = myAtom("CB",resnum,resname);
    CG = myAtom("CG",resnum,resname);
    CD = myAtom("CD",resnum,resname);
    OE = myAtom("OE",resnum,resname);
    HB = myAtom("HB#",resnum,resname);
    HG = myAtom("HG#",resnum,resname);

    # Create expected bonds
    expected_edges = union!(expected_edges_backbone,Set([
    (CA,HA),
    (CA,CB),
    (CB,CG),
    (CG,CD),
    (CD,OE),
    (CB,HB),
    (CG,HG)
    ]))

    expected_enriched_edges = Set([
    (OE,CD),
    (C,CB),
    (CB,HG),
    (H,C),
    (CD,CB),
    (CA,H),
    (CD,HG),
    (HA,H),
    (N,HB),
    (HB,C),
    (N,HA),
    (CA,HB),
    (CG,CB),
    (CD,CG),
    (HA,C),
    (N,CA),
    (HB,CG),
    (OE,HG),
    (HB,CB),
    (O,C),
    (CG,HG),
    (HA,CG),
    (CA,CD),
    (HB,HA),
    (CD,HB),
    (CA,CG),
    (HB,HG),
    (CA,HA),
    (O,CB),
    (N,C),
    (CG,C),
    (N,CB),
    (N,CG),
    (N,H),
    (CA,CB),
    (HA,CB),
    (CA,HG),
    (N,O),
    (CA,C),
    (OE,CB),
    (H,CB),
    (CA,O),
    (HA,O),
    (OE,CG)
    ])

    # Check that the graph is not enriched
    @test sets_of_tuples_equal(Set(elist(g)), expected_edges)

    # Check that the graph is enriched
    enrich!(g)
    @test sets_of_tuples_equal(Set(elist(g)), expected_enriched_edges)
end