# This file is used to define some utility functions for testing

# print_to_string function
# This function is used to print the output of a function with io as first parameter to a string
print_to_string(f::Function, arg::Any) = (io = IOBuffer(); f(io,arg); String(take!(io)))

# print_to_string_stdout function
# This function is used to print the output of a function using the standardOutput to a string
function print_to_string_stdout(f::Function, args::Any...)
    stdout = Base.stdout
    redirect_stdout = IOBuffer()
    try
        # Redirect stdout to the buffer
        Base.stdout = redirect_stdout
        f(args)  # Call the function with its argument
    catch e
        # If an error occurs, ensure the original stdout is restored
        Base.stdout = stdout
        rethrow(e)  # Re-throw the caught error
    finally
        # Always restore the original stdout, even if an error occurs
        Base.stdout = stdout
    end
    output = String(take!(redirect_stdout))
    return output
end

function print_to_string_stdout(f::Function, arg::Any)
    stdout = Base.stdout
    redirect_stdout = IOBuffer()
    try
        # Redirect stdout to the buffer
        Base.stdout = redirect_stdout
        f(arg)  # Call the function with its argument
    catch e
        # If an error occurs, ensure the original stdout is restored
        Base.stdout = stdout
        rethrow(e)  # Re-throw the caught error
    finally
        # Always restore the original stdout, even if an error occurs
        Base.stdout = stdout
    end
    output = String(take!(redirect_stdout))
    return output
end

# suppress_stdout function
# This function calls another function but doesn't allow it to print anything to the standardOutput
function suppress_stdout(f::Function, args::Any...)
    stdout = Base.stdout
    redirect_stdout = IOBuffer()
    try
    Base.stdout = redirect_stdout
    f(args...)
    catch e
    Base.stdout = stdout
    rethrow(e)
    finally
    Base.stdout = stdout
    end
end

function suppress_stdout(f::Function, arg::Any)
    stdout = Base.stdout
    redirect_stdout = IOBuffer()
    try
    Base.stdout = redirect_stdout
    f(arg)
    catch e
    Base.stdout = stdout
    rethrow(e)
    finally
    Base.stdout = stdout
    end
end

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

#Create a STARfile with a readable format for DGP(STARfile::String) from a DGP instance with the given filename
function create_pseudo_STAR(dgp::DGP{myAtom}, starfile::String)
    open(starfile, "w") do f
        # Write the header
        println(f, "loop_")
        # Write the identifiers
        println(f, "_Gen_dist_constraint.PDB_residue_no_1")
        println(f, "_Gen_dist_constraint.PDB_residue_no_2")
        println(f, "_Gen_dist_constraint.PDB_residue_name_1")
        println(f, "_Gen_dist_constraint.PDB_residue_name_2")
        println(f, "_Gen_dist_constraint.Auth_atom_ID_1")
        println(f, "_Gen_dist_constraint.Auth_atom_ID_2")
        println(f, "_Gen_dist_constraint.Distance_lower_bound_val")
        println(f, "_Gen_dist_constraint.Distance_upper_bound_val")
        # Write the distances
        for ((a, b), dist) in dgp.distlist
            println(f, join([
                a.resnum,
                b.resnum,
                a.resname,
                b.resname,
                a.name,
                b.name,
                dist.lb,
                dist.ub
            ], " "))
        end
        # End the loop
        println(f, "stop_")
    end
end

#Create a DGP from a MDjeep file
function create_DGP_from_MDjeep(filename::String)
    # Check if the file exists
    if !isfile(filename)
        throw(ArgumentError("The file $filename does not exist"))
    end
    # Initialize data structures
    distlist = Dict{Tuple{myAtom, myAtom}, Distance}()
    order = Set{myAtom}()
    # Open the MDjeep file and read its content
    open(filename, "r") do file
        for line in eachline(file)
            # Parse the line
            parts = split(line)
            if length(parts) != 10
                throw(ArgumentError("FORMAT ERROR: $line"))
            end
            # Extract values
            j = parse(Int, parts[1])
            i = parse(Int, parts[2])
            lb = parse(Float64, parts[3])
            ub = parse(Float64, parts[4])
            name1 = String(parts[5])
            name2 = String(parts[6])
            resnum1 = parse(Int64, parts[7])
            resnum2 = parse(Int64, parts[8])
            resname1 = String(parts[9])
            resname2 = String(parts[10])
            # Create myAtom instances
            atom1 = myAtom(name1, resnum1, resname1)
            atom2 = myAtom(name2, resnum2, resname2)
            # Add atoms to the order set
            push!(order, atom1)
            push!(order, atom2)
            # Create and add the distance
            D = Distance(lb, ub)
            distlist[(atom2, atom1)] = D
        end
    end
    # Convert the set to a vector to maintain order
    ordered_atoms = collect(order)
    # Create and return the DGP instance
    return DGP(3, distlist, ordered_atoms)
end

