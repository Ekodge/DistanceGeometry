
# DGP Julia project (see JuliaDGP.jl)
#
# DGP type
#
# importing Base.show
# using DataStructures
# including DistList, Orders, Realization, myAtom, Utils
#
# AM

struct DGP{T} <: DistList{T}
   K::Int64;  # dimension
   distlist::Dict{Tuple{T,T},Distance};  # list of distances
   order::Vector{T};  # order for the elements of type T

   # constructor for random fully-connected instance from given Realization
   function DGP(R::Realization{T}) where {T}
      order = vertex_order(R);

      # computing all the distances
      n = length(order);
      distlist = Dict{Tuple{T,T},Distance}();
      for i in 1:n
         u = order[i];
         for j in 1:i - 1
            v = order[j];
            push!(distlist,(u,v) => Distance(R.coords[u],R.coords[v]));
         end
      end

      # DGP instance is ready
      new{T}(R.K,distlist,order);
   end

   # constructor from a given Realization instance and distance cut-off
   function DGP(R::Realization{T},threshold::Float64) where {T}
      if (threshold <= 0.0) throw(ArgumentError("The specified threshold is 0.0 or even negative")) end
      order = vertex_order(R);

      # computing and selecting distances
      n = length(order);
      distlist = Dict{Tuple{T,T},Distance}();
      for i = 1:n
         u = order[i];
         for j in 1:i - 1
            v = order[j];
            dist = Distance(R.coords[u],R.coords[v]);
	    if dist.ub <= threshold
               push!(distlist,(u,v) => dist);
            end
         end
      end

      # DGP instance is ready
      new{T}(R.K,distlist,order);
   end

   # constructor for a paradoxical DGP from a given Realization instance
   # -> "cut" is the (constant) number of edges crossing DGP elements
   function DGP(R::Realization{T},cut::Int64) where {T}
      n = R.N;
      if (cut <= 0) throw(ArgumentError("THe specified cut value is 0 or even negative")) end
      if (cut > n) throw(ArgumentError("The value of cut cannot be larger than the Realization size")) end
      if (cut == n) return DGP(R) end  # it's going to be a complete DGP instance

      # preparing the vectors of elements and references
      elements = vertex_order(R);
      refs = Deque{T}();
      for i in n-cut+1:n
         push!(refs,elements[i]);
      end

      # constructing the paradoxical instance
      distlist = Dict{Tuple{T,T},Distance}();
      for i in 1:n
         v = elements[i];
         for u in refs
            push!(distlist,(u,v) => Distance(R.coords[u],R.coords[v]));
         end
         popfirst!(refs);
         push!(refs,v);
      end

      # paradoxical DGP is ready
      new{T}(R.K,distlist,elements);
   end

   # constructor from a given Realization instance and its connection graph
   function DGP(R::Realization{T},g::UndirectedGraph{T}) where {T}
      elements = vertex_order(R);  # initial order with *all* vertices

      # selecting only distances for which an edge is in the graph g
      n = length(elements);
      order = OrderedDict{Int64,T}();
      distlist = Dict{Tuple{T,T},Distance}();
      for i = 1:n
         u = elements[i];
         for j in 1:i - 1
            v = elements[j];
            if (u,v) in g.E || (v,u) in g.E
               push!(order,i => u);
               push!(order,j => v);
               push!(distlist,(u,v) => Distance(R.coords[u],R.coords[v]));
            end
         end
      end

      # DGP instance is ready
      new{T}(R.K,distlist,collect(values(order)));
   end

   # constructor from STAR file (NMR data, we ignore uncertain distances)
   function DGP(STARfile::String)
      ex = get_extension(STARfile);
      if ex != "str"
         throw(ArgumentError("Input file is supposed to be a STAR file; check if extension is coherent with content"));
      end

      # preparing the data structures
      order = Vector{myAtom}();
      distlist = Dict{Tuple{myAtom,myAtom},Distance}();

      # opening the STAR file
      star = open(STARfile);

      # reading the lines of the STAR file
      extracting = false;
      identifiers = Vector{String}();
      for line in readlines(star)
         if !extracting
            # looking for the blocks of data
            if contains(line,"loop_")
               extracting = true;
            end
         else
            if contains(line,"_Gen_dist_constraint.")
               # reading the identifiers
               splitline = split(line,'.');
               id = String(splitline[2]);
               push!(identifiers,id);
            elseif !isempty(identifiers)
               # reading the data if the searched identifiers were found
               n = length(identifiers);
               splitline = split(line);
               if n == length(splitline)
                  # preparing some local variables
                  resnum1 = nothing;  resnum2 = nothing;
                  resname1 = nothing; resname2 = nothing;
                  name1 = nothing;    name2 = nothing;
                  lb = nothing;       ub = nothing;
                  uncertain = false;

                  # reading
                  for i in 1:n
                     if contains(identifiers[i],"Member_logic_code")
                        if (String(splitline[i]) == "OR") uncertain = true; break end  # we skip it!
                     elseif contains(identifiers[i],"PDB_residue_no_1")
                        resnum1 = parse(Int64,String(splitline[i]));
                     elseif contains(identifiers[i],"PDB_residue_no_2")
                        resnum2 = parse(Int64,String(splitline[i]));
                     elseif contains(identifiers[i],"PDB_residue_name_1")
                        resname1 = String(splitline[i]);
                     elseif contains(identifiers[i],"PDB_residue_name_2")
                        resname2 = String(splitline[i]);
                     elseif contains(identifiers[i],"Auth_atom_ID_1")
                        name1 = String(splitline[i]);
                     elseif contains(identifiers[i],"Auth_atom_ID_2")
                        name2 = String(splitline[i]);
                     elseif contains(identifiers[i],"Distance_lower_bound_val")
                        lb = parse(Float64,String(splitline[i]));
                     elseif contains(identifiers[i],"Distance_upper_bound_val")
                        ub = parse(Float64,String(splitline[i]));
                     end
                  end

                  # defining the data structures
                  if !uncertain
                     if (resnum1 == nothing || resnum2 == nothing) throw(ArgumentError("Something went wrong...")) end
                     if (resname1 == nothing || resname2 == nothing) throw(ArgumentError("Something went wrong...")) end
                     if (name1 == nothing || name2 == nothing) throw(ArgumentError("Something went wrong...")) end
                     if (name1 == "HN")  name1 = "H" end
                     a = myAtom(name1,resnum1,resname1);
                     if !(a in order) push!(order,a) end
                     if (name2 == "HN")  name2 = "H" end
                     b = myAtom(name2,resnum2,resname2);
                     if !(b in order) push!(order,b) end
                     if (lb == nothing || ub == nothing) throw(ArgumentError("Something went wrong...")) end
                     D = Distance(lb,ub);
                     push!(distlist,(a,b) => D);
                  end
               end
            elseif contains(line,"stop_")
               # stopping
               identifiers = Vector{String}();
               extracting = false;
            end
         end
      end

      # closing the STAR file
      close(star);

      # DGP instance is ready
      new{myAtom}(3,distlist,order);
   end

   #DGP constructor from given dimension, distlist and order, for testing purposes
   function DGP(K::Int64,distlist::Dict{Tuple{T,T},Distance},order::Vector{T}) where {T}
       new{T}(K,distlist,order);
   end

   # overriding Base show function
   function Base.show(io::IO,dgp::DGP{T}) where {T}
      n = length(dgp.order);
      m = nd(dgp);
      print(io,"DGP{",nameof(T),"} (K = ",dgp.K," with ",n," elements and ",m," distances)");
   end
end

# showing the details of a DGP instance
function details(dgp::DGP{T}) where {T}
   println("DGP{",nameof(T),"} (K = ",dgp.K,") {");
   for u in dgp.order
      for v in dgp.order
         if u != v
            if haskey(dgp.distlist,(u,v))
               println("  (",u,",",v,") => ",dgp.distlist[(u,v)]);
            end
         end
      end
   end
   print("}");
end

# finding the references for a given element v of a DGP instance
function references(dgp::DGP{T},v::T) where {T}
   if !(v in dgp.order) throw(ArgumentError("The given element does not belong to the DGP instance")) end

   # looking for the position of the element in the order
   k = findfirst(x -> x == v,dgp.order);

   # looking for the references
   refs = Vector{T}();
   for i in 1:k - 1
      u = dgp.order[i];
      if haskey(dgp.distlist,(u,v)) || haskey(dgp.distlist,(v,u))
         push!(refs,u);
      end
   end

   # reference vector is ready
   return refs;
end

# printing in MDjeep format 
function MDjeep_format(dgp::DGP{T}) where {T}
   file = open("test_instance.nmr","w");
   len = length(dgp.order);
   for i in 1:len
      u = dgp.order[i];
      for j in 1:i - 1
         v = dgp.order[j];
         D = nothing;
         if (haskey(dgp.distlist,(u,v))) D = dgp.distlist[(u,v)] end
         if (haskey(dgp.distlist,(v,u))) D = dgp.distlist[(v,u)] end
         if D != nothing
            line = string(j) * " " * string(i) * " " * @sprintf("%20.16f",D.lb) * " " * @sprintf("%20.16f",D.ub);
	    if typeof(v) == myAtom && typeof(u) == myAtom
               line = line * " " * v.name * " " * u.name;
               line = line * " " * string(v.resnum) * " " * string(u.resnum);
               line = line * " " * v.resname * " " * u.resname;
            end
            line = line * "\n";
            write(file,line);
         end
      end
   end
   println("File 'test_instance.nmr' written in MDjeep format.");
   close(file);
end

# number of edges associated with one vertex within a dgp instance
function count_edges(dgp::DGP{T}, vertex::T) where {T}
    count = 0
    for (u, v) in keys(dgp.distlist)
        if u == vertex || v == vertex
            count += 1
        end
    end
    return count
end

#recursive function that returns a list of all possible label matching
function generateAllLabelsMatching(possible_labels_matching::Dict{Int64, Vector{Int64}}, current_match::Dict{Int64,Int64})

	current_candidate = -1

	# if all labels are matched then return current_match
	if length(current_match) == length(possible_labels_matching)
		return [current_match]
	end

	# Get the label that isn't already matched
	for i in 1:length(possible_labels_matching)
    	if !haskey(current_match, i)
        	current_candidate = i
        	break
    	end
	end

	all_labels_matching = Vector{Dict{Int64, Int64}}()

    # Call the function recursively for each possible matching of the current label, empty list if no possible matching
	for label in possible_labels_matching[current_candidate]
		current_match_copy = deepcopy(current_match)
		current_possible_matching_copy = deepcopy(possible_labels_matching)
		current_match_copy[current_candidate] = label
		#remove the label from the other labels
		for j in 1:length(possible_labels_matching)
			if current_candidate != j
				current_possible_matching_copy[j] = filter!(x -> x != label, current_possible_matching_copy[j])
				# If any label ends up with no possible matches, return an empty list
                if isempty(current_possible_matching_copy[j])
    				return []
				end
			end
		end
		# add the possible matching to the list of possible matching
		all_labels_matching = [all_labels_matching; generateAllLabelsMatching(current_possible_matching_copy, current_match_copy)]
	end

	return all_labels_matching
end

# the various equals operators of DGP

#override of the '==' operator
#two DGP instances satisfy the '==' operator if they have the same dimension, number of elements,
#the same structure and the same distances disregarding the order of the elements and what they represent. (isomorphic)
function Base.:(==)(dgp1::DGP{T},dgp2::DGP{T}) where {T}
   # checking the dimensions
   if dgp1.K != dgp2.K
        return false
   end

   # checking the number of vertices
   if length(dgp1.order) != length(dgp2.order)
        return false
   end

    numOfVertices = length(dgp1.order);

    # creating labels for the vertices
    labels1 = Dict{Int64,T}()
    labels2 = Dict{Int64,T}()
	#dgp1 to possible dgp2 labels
    possible_labels_matching = Dict{Int64, Vector{Int64}}()

    for i in 1:numOfVertices
        labels1[i] = dgp1.order[i]
        labels2[i] = dgp2.order[i]
    end

    # Check for each label if they can find at least one corresponding label in the other DGP
    # by checking how many edges they are included in
    for i in 1:numOfVertices
        referenceCount1 = count_edges(dgp1, labels1[i])
        possible_labels_matching[i] = Vector{Int64}()

        for j in 1:numOfVertices
            referenceCount2 = count_edges(dgp2, labels2[j])

            if referenceCount1 == referenceCount2
                push!(possible_labels_matching[i], j)
            end
        end

        # If no matching label was found for this vertex, the graphs are not equal
        if isempty(possible_labels_matching[i])
            return false
        end
    end

	fixed_matches = Dict{Int64, Int64}()

    # Remove the labels that have already been matched (meaning only one match) and then check if the graphs
	#can still be equal
	for i in 1:numOfVertices
		if length(possible_labels_matching[i]) == 1
			matched_Label = possible_labels_matching[i][1]
			fixed_matches[i] = matched_Label
			#remove the label from the other labels
			for j in 1:numOfVertices
				if i != j
					possible_labels_matching[j] = filter!(x -> x != matched_Label, possible_labels_matching[j])
					# If no matching label was found for this vertex after filter, the graphs are not equal
					if isempty(possible_labels_matching[i])
        				return false
					end
				end
			end
		end
	end

	# Get all different permutation of possible matching between the two lists of labels
	all_labels_matching = generateAllLabelsMatching(possible_labels_matching, fixed_matches)

	# If no possible matching was found, the graphs are not equal
	if isempty(all_labels_matching)
		return false
	end

	found = false
	# Check if the two graphs are equal for at least one of the possible label matching
	# by checking if the distances are equal
	for possible_match in all_labels_matching
		for i in 1:numOfVertices
			dgp1_object = labels1[i]
			dgp2_object = labels2[possible_match[i]]
			# put in a list all the distances of the first DGP
			distances1 = Vector{Distance}()
			for (u, v) in keys(dgp1.distlist)
        		if u == dgp1_object || v == dgp1_object
            		push!(distances1, dgp1.distlist[(u,v)])
				end
        	end
			# put in a list all the distances of the second DGP
			distances2 = Vector{Distance}()
			for (u, v) in keys(dgp2.distlist)
				if u == dgp2_object || v == dgp2_object
					push!(distances2, dgp2.distlist[(u,v)])
				end
			end
			# iterate over the distances of the first DGP and check if they are equal to the distances of the second DGP
			for distance1 in distances1
				index = findfirst(x -> x == distance1, distances2)
    			if index !== nothing
				deleteat!(distances2, index)
				end
			end
			# if distances2 is empty, the two DGP are equal
			if isempty(distances2)
				found = true
				break
			end
		end
	end
	return found
end

#override of the '!=' operator
#two DGP instances satisfy the '!=' operator if they do not satisfy the '==' operator
function Base.:(!=)(dgp1::DGP{T},dgp2::DGP{T}) where {T}
   return !(dgp1 == dgp2);
end

#definition of the '===' operator
#two DGP instances satisfy the '===' operator if satisfy the '==' operator and if the
#order of the elements is the same in both instances disregarding what the elements actually are.
function equalOrder(dgp1::DGP{T}, dgp2::DGP{T}) where {T}
   # checking the dimensions
   if dgp1.K != dgp2.K
        return false
   end

   # checking the number of vertices
   if length(dgp1.order) != length(dgp2.order)
        return false
   end

    numOfVertices = length(dgp1.order);

    # checking that the number of edges are the same for each vertex
	for i in 1:numOfVertices
		if count_edges(dgp1, dgp1.order[i]) != count_edges(dgp2, dgp2.order[i])
			return false
		end
	end

	# checking distances are the same for all matching vertices
	for i in 1:numOfVertices
		dgp1_object = dgp1.order[i]
		dgp2_object = dgp2.order[i]
		# put in a list all the distances of the first DGP
		distances1 = Vector{Distance}()
		for (u, v) in keys(dgp1.distlist)
			if u == dgp1_object || v == dgp1_object
				push!(distances1, dgp1.distlist[(u,v)])
			end
		end
		# put in a list all the distances of the second DGP
		distances2 = Vector{Distance}()
		for (u, v) in keys(dgp2.distlist)
			if u == dgp2_object || v == dgp2_object
				push!(distances2, dgp2.distlist[(u,v)])
			end
		end
		# iterate over the distances of the first DGP and check if they are equal to the distances of the second DGP
		for distance1 in distances1
			index = findfirst(x -> x == distance1, distances2)
			if index !== nothing
			deleteat!(distances2, index)
			end
		end
		# if distances2 is empty, the two DGP are equal
		if !isempty(distances2)
			return false
		end
	end

	return true
end

#definition of the '!==' operator
#two DGP instances satisfy the '!==' operator if they do not satisfy the '===' operator
function notEqualOrder(dgp1::DGP{T},dgp2::DGP{T}) where {T}
   return !equalOrder(dgp1,dgp2)
end

#definition of the '====' operator
#two DGP instances satisfy the '====' operator if they satisfy the '===' operator and if all the
#elements of the two instances are the exact same.
function equals(dgp1::DGP{T},dgp2::DGP{T}) where {T}
	#checking if the two instances satisfy the '===' operator (they are isomorphic)
	if (notEqualOrder(dgp1,dgp2))
    	return false
   end

   numOfVertices = length(dgp1.order)

   #checking if the two instances have the same elements
   for i in 1:numOfVertices
   		if dgp1.order[i] != dgp2.order[i]
   			return false
   		end
   end

   return true
end


#definition of the '!====' operator
#two DGP instances satisfy the '!====' operator if they do not satisfy the '====' operator
function notEquals(dgp1::DGP{T},dgp2::DGP{T}) where {T}
   return !(dgp1 ≡ dgp2);
end
