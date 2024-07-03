# This file is used to define some utility functions for testing

# print_to_string function
# This function is used to print the output of a function to a string
print_to_string(f::Function, arg::Any) = (io = IOBuffer(); f(io,arg); String(take!(io)))

# sets_of_vectors_equal function
# This function is used to compare two sets of tuples{T,T} (not caring about order)
function sets_of_tuples_equal(set1::Set{Tuple{T,T}}, set2::Set{Tuple{T,T}}) where T
    set1_converted = Set([Set(t) for t in set1])
    set2_converted = Set([Set(t) for t in set2])
    return set1_converted == set2_converted
end

# set_of_backbone_atoms function
# This function is used to get the set of backbone atoms and bonds from a resnum and a resname
function backbone_atoms_bonds(resnum::Int, resname::String)
    N = myAtom("N",resnum, resname)
    H = myAtom("H",resnum, resname)
    CA = myAtom("CA",resnum, resname)
    C = myAtom("C",resnum, resname)
    O = myAtom("O",resnum, resname)
    return (atoms = Set([N,H,CA,C,O]), bonds = Set([
    (N,H),
    (N,CA),
    (CA,C),
    (C,O)
    ]))
end

# Convert amino acid sequences to expected sets of atoms and bonds
function expected_protein_graph(sequence::String)
    n = length(sequence)
    if n == 0
        throw(ArgumentError("The protein sequence has zero length"))
    end

    map = fun_dict()
    rescode = sequence[1]
    g = getfield(Main, Symbol(map[rescode]))(1)
    expected_atoms = Set(vlist(g))
    expected_bonds = Set(elist(g))

    for i in 2:n
        prev = rescode
        rescode = sequence[i]
        getfield(Main, Symbol(map[rescode] * "!"))(g, i)

        C = myAtom("C", i - 1, aa_name_convert(string(prev)))
        N = myAtom("N", i, aa_name_convert(string(rescode)))
        add!(g, C, N)

        expected_atoms = union!(expected_atoms, Set(vlist(g)))
        expected_bonds = union!(expected_bonds, Set(elist(g)))
    end

    return (expected_atoms, expected_bonds)
end

#Create Atoms following the PDBTool atom struct
function create_sample_PDBTool_atoms()
    return [
        Atom(index=1, index_pdb=1, name="N", resname="ALA", chain="A", resnum=1, x=1.0, y=2.0, z=3.0),
        Atom(index=2, index_pdb=2, name="CA", resname="ALA", chain="A", resnum=1, x=4.0, y=5.0, z=6.0),
        Atom(index=3, index_pdb=3, name="C", resname="ALA", chain="A", resnum=1, x=7.0, y=8.0, z=9.0)
    ]
end


