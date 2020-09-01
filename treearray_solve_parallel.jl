function check_tas_parallel(tas; post = "/scratch/smkatz/post2d/")
    updated_tas = Dict()
    tas_rc = Dict()
    # Do the tau = 0 trees for each pra
    println("τ: 0")
    for pra in actions
        ta = tas[(pra, 0.0)]
        post_prefix = "$(post)pra$(pra+1)tau0"
        tas_rc[pra] = remotecall(label_collisions, mod(pra + 1, nprocs() - 1) + 2, ta, post_prefix)
    end
    for pra in actions
        updated_tas[(pra, 0.0)] = fetch(tas_rc[pra])
    end
    # for the rest of the taus
    for τ = 1:40
        println("τ: $τ")
        τint = convert(Int64, τ)
        # do update for each pra
        for pra in actions
            ta = tas[(pra, τ)]
            post_prefix = "$(post)pra$(pra+1)tau$(τint)"
            tas_rc[pra] = remotecall(update, mod(pra + 1, nprocs() - 1) + 2, ta, updated_tas, pra, τ, post_prefix)
        end
        for pra in actions
            updated_tas[(pra, τ)] = fetch(tas_rc[pra])
        end
    end
    return updated_tas
end

function check_and_split_tas_parallel(tas; qthreshold = 0.5, threshold = 0.9,
    post = "/scratch/smkatz/post2d/", prefix = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra0")
    
    updated_tas = Dict()
    tas_rc = Dict()
    # Do the tau = 0 trees for each pra
    println("τ: 0")
    for pra in actions
        ta = tas[(pra, 0.0)]
        post_prefix = "$(post)pra$(pra+1)tau0"
        tas_rc[pra] = remotecall(label_collisions, mod(pra + 1, nprocs() - 1) + 2, ta, post_prefix)
    end
    for pra in actions
        updated_tas[(pra, 0.0)] = fetch(tas_rc[pra])
    end

    # allocate some memory to store temp tas
    # (I think this will help with memory)
    # temp_tas = Dict()
    # tot_size = 200000
    # for pra in actions
    #     temp_tas[pra] = treearray(tot_size)
    # end

    # for the rest of the taus
    for τ = 1:40
        start = time()
        println("τ: $τ")
        τint = convert(Int64, τ)
        # do update for each pra
        for pra in actions
            ta = tas[(pra, τ)]
            post_prefix = "$(post)pra$(pra+1)tau$(τint)"
            tas_rc[pra] = remotecall(update_or_split, mod(pra + 1, nprocs() - 1) + 2, 
                                        ta, updated_tas, pra, τ, post_prefix,
                                        qthreshold, threshold, prefix)
        end
        for pra in actions
            updated_tas[(pra, τ)] = fetch(tas_rc[pra])
        end
        println("Time: $(time() - start)")
    end
end

function update_parallel(pra, τ, post)   
    τint = convert(Int64, τ)

    prefix = "/scratch/smkatz/pre2d/pra$(pra+1)tau$(τint)"
    ta = mmap_tree(prefix)

    post_prefix = "$(post)pra$(pra+1)tau$(τint)"

    updated_ta = treearray_copy_and_map(ta, 40000, post_prefix)
    #updated_ta = treearray_copy_and_extend(ta, 40000)
    
    nodes = updated_ta.nodes
    splits = updated_ta.splits
    dims = updated_ta.dims
    leaf_data = updated_ta.leaf_data
    qvals = updated_ta.qvals
    
    cats = []
    lbs = []
    ubs = []
    
    lb_s = Stack{Vector{Float64}}()
    ub_s = Stack{Vector{Float64}}()
    s = Stack{Int32}()

    push!(lb_s, [-0.5, -0.5])
    push!(ub_s, [0.5, 0.5])
    push!(s, 1)

    while !isempty(s)
        curr = pop!(s)
        curr_lbs = pop!(lb_s)
        curr_ubs = pop!(ub_s)

        if nodes[curr] < 0 # leaf node
            # Get the possible actions
            cats = findall(leaf_data[:, -nodes[curr]]) .- 1
            for i = 1:length(cats)
                next_prefix = "$(post)pra$(cats[i]+1)tau$(τint)"
                next_ta = mmap_tree(next_prefix)
                _, vals, probs = transition_ta(next_ta, curr_lbs, curr_ubs, pra)
                for j = 1:length(vals)
                    qvals[cats[i] + 1, -nodes[curr]] += probs[j]*vals[j]
                end
            end
        else # interior node
            # Traverse tree and keep track of bounds
            dim = findall(dims[:, curr])[1]
            split = splits[curr]
            # Go left, upper bounds will change
            left_ubs = copy(curr_ubs)
            left_ubs[dim] = split

            push!(lb_s, curr_lbs)
            push!(ub_s, left_ubs)
            push!(s, nodes[curr])

            # Go right, lower bounds will change
            right_lbs = copy(curr_lbs)
            right_lbs[dim] = split
            
            push!(lb_s, right_lbs)
            push!(ub_s, curr_ubs)
            push!(s, nodes[curr] + 1)
        end
    end

    # Write the final tree array to a file and memory map it
    # write_to_files_and_truncate(temp_ta, post_prefix)
    # mapped_updated_ta = mmap_tree(post_prefix)
    # return mapped_updated_ta
    sync_treearray(updated_ta)
    return true
end