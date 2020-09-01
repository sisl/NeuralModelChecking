function check_tas(tas; post = "/scratch/smkatz/post2d/")
    updated_tas = Dict()
    # Do the tau = 0 trees for each pra
    println("τ: 0")
    for pra in actions
        ta = tas[(pra, 0.0)]
        post_prefix = "$(post)pra$(pra+1)tau0"
        updated_tas[(pra, 0.0)] = label_collisions(ta, post_prefix)
    end
    # for the rest of the taus
    for τ = 1:40
        println("τ: $τ")
        τint = convert(Int64, τ)
        # do update for each pra
        for pra in actions
            ta = tas[(pra, τ)]
            post_prefix = "$(post)pra$(pra+1)tau$(τint)"
            updated_tas[(pra, τ)] = update(ta, updated_tas, pra, τ, post_prefix)
        end
    end
end

function label_collisions(ta::TREEARRAY, post_prefix)
    collision_h = 100 / 16000
    
    updated_ta = treearray_copy(ta)
    
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
            # Label with probability 1 if it indicates a collision, 0 else
            collision = abs(curr_lbs[1]) < collision_h || abs(curr_ubs[1]) < collision_h
            leaf_ind = -nodes[curr]
            action_inds = findall(leaf_data[:, leaf_ind])

            q = zeros(9)
            if collision
                q[action_inds] .= 1.0
            end
            qvals[:, leaf_ind] = q
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
    write_to_files(updated_ta, post_prefix)
    mapped_updated_ta = mmap_tree(post_prefix)
    return mapped_updated_ta
end

function update(ta::TREEARRAY, updated_tas, pra, τ, post_prefix)   
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
                next_ta = updated_tas[cats[i], τ - 1]
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
    return updated_ta
end

"""
Check and split
"""

function check_and_split_tas(tas; qthreshold = 0.5, threshold = 0.9,
    post = "/scratch/smkatz/post2d/", prefix = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra0")
    
    updated_tas = Dict()
    # Do the tau = 0 trees for each pra
    println("τ: 0")
    for pra in actions
        ta = tas[(pra, 0.0)]
        post_prefix = "$(post)pra$(pra+1)tau0"
        updated_tas[(pra, 0.0)] = label_collisions(ta, post_prefix)
    end
    # for the rest of the taus
    for τ = 1:40
        start = time()
        println("τ: $τ")
        τint = convert(Int64, τ)
        # do update for each pra
        for pra in actions
            ta = tas[(pra, τ)]
            post_prefix = "$(post)pra$(pra+1)tau$(τint)"
            updated_tas[(pra, τ)] = update_or_split(ta, updated_tas, pra, τ, post_prefix,
                                        qthreshold, threshold, prefix)
        end
        println("Time: $(time() - start)")
    end
end

function update_or_split(ta::TREEARRAY, updated_tas, pra, τ, post_prefix, qthreshold, threshold, prefix)
    updated_ta = treearray_copy_and_map(ta, 400000, post_prefix)
    #updated_ta = treearray_copy_and_extend(ta, 80000)
    #copy_over!(ta, temp_ta)
    
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

        if updated_ta.nodes[curr] < 0 # leaf node
            # Get the possible actions
            update_leaf!(updated_ta, updated_tas, pra, τ, curr, curr_lbs, curr_ubs, s, lb_s, ub_s, 
                            qthreshold = qthreshold, threshold = threshold, prefix = prefix)
        else # interior node
            # Traverse tree and keep track of bounds
            dim = findall(updated_ta.dims[:, curr])[1]
            split = updated_ta.splits[curr]
            # Go left, upper bounds will change
            left_ubs = copy(curr_ubs)
            left_ubs[dim] = split

            push!(lb_s, curr_lbs)
            push!(ub_s, left_ubs)
            push!(s, updated_ta.nodes[curr])

            # Go right, lower bounds will change
            right_lbs = copy(curr_lbs)
            right_lbs[dim] = split
            
            push!(lb_s, right_lbs)
            push!(ub_s, curr_ubs)
            push!(s, updated_ta.nodes[curr] + 1)
        end
    end

    # Write the final tree array to a file and memory map it
    # write_to_files_and_truncate(temp_ta, post_prefix)
    # mapped_updated_ta = mmap_tree(post_prefix)
    # return mapped_updated_ta
    sync_treearray(updated_ta)
    return updated_ta
end

function update_leaf!(ta::Union{TREEARRAY, SHARED_TREEARRAY}, updated_tas, pra, τ, curr, curr_lbs, curr_ubs, s, lb_s, ub_s; qthreshold = 0.5, threshold = 0.9,
    prefix = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra0")
    
    nodes = ta.nodes
    splits = ta.splits
    dims = ta.dims
    leaf_data = ta.leaf_data
    qvals = ta.qvals

    cats = findall(leaf_data[:, -nodes[curr]]) .- 1
    for i = 1:length(cats)
        next_ta = updated_tas[cats[i], τ - 1]
        max_range, vals, probs = transition_ta(next_ta, curr_lbs, curr_ubs, pra)
        
        should_split = (max_range > threshold) && (curr_ubs[1] - curr_lbs[1] > 10 / 16000)
        if should_split
            split_and_reverify!(ta, curr, curr_lbs, curr_ubs, s, lb_s, ub_s, pra, τ, 1, prefix = prefix)
            break
        else
            for j = 1:length(vals)
                qvals[cats[i] + 1, -nodes[curr]] += probs[j]*vals[j]
            end

            # Check if need to split due to actions
            qval_diff = maximum(qvals[:,-nodes[curr]]) - minimum(qvals[:,-nodes[curr]])
            split_h = (qval_diff > qthreshold) && (curr_ubs[1] - curr_lbs[1] > 10 / 16000)
            split_ḣ₀ = (qval_diff > qthreshold) && (curr_ubs[2] - curr_lbs[2] > 1 / 200)

            if split_h && split_ḣ₀
                split_and_reverify!(ta, curr, curr_lbs, curr_ubs, s, lb_s, ub_s, pra, τ, 1, prefix = prefix)

                right = pop!(s)
                right_lbs = pop!(lb_s)
                right_ubs = pop!(ub_s)

                left = pop!(s)
                left_lbs = pop!(lb_s)
                left_ubs = pop!(ub_s)

                split_and_reverify!(ta, left, left_lbs, left_ubs, s, lb_s, ub_s, pra, τ, 2, prefix = prefix)
                split_and_reverify!(ta, right, right_lbs, right_ubs, s, lb_s, ub_s, pra, τ, 2, prefix = prefix)
                break
            elseif split_h
                split_and_reverify!(ta, curr, curr_lbs, curr_ubs, s, lb_s, ub_s, pra, τ, 1, prefix = prefix)
                break
            elseif split_ḣ₀
                split_and_reverify!(ta, curr, curr_lbs, curr_ubs, s, lb_s, ub_s, pra, τ, 2, prefix = prefix)
                break
            end
        end
    end
    #sync_treearray(ta)
end

function split_and_reverify!(ta::Union{TREEARRAY, SHARED_TREEARRAY}, curr, lbs, ubs, s, lb_s, ub_s, pra, τ, dim;
    prefix = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra0")
    
    cats = findall(ta.leaf_data[:, -ta.nodes[curr]]) .- 1 # subtract 1 because index at 1

    left_cats = cats
    right_cats = cats

    split = (lbs[dim] + ubs[dim]) ./ 2

    # Go left, upper bounds will change
    left_ubs = copy(ubs)
    left_ubs[dim] = split
    # Go right, lower bounds will change
    right_lbs = copy(lbs)
    right_lbs[dim] = split

    # Get categories for new leaves
    if length(cats) > 1
        left_cats = get_categories!(cats, lbs, left_ubs, pra, τ, prefix = prefix)
        if length(left_cats) == 0
            left_cats = get_categories!(actions, lbs, left_ubs, pra, τ, prefix = prefix)
        end

        right_cats = get_categories!(cats, right_lbs, ubs, pra, τ, prefix = prefix)
        if length(right_cats) == 0
            right_cats = get_categories!(actions, right_lbs, ubs, pra, τ, prefix = prefix)
        end
    end

    # Split the node

    # Reset values for curr node
    ta.splits[curr] = split
    ta.dims[dim, curr] = true
    ta.nodes[curr] = ta.node_index

    # Set values for leaf nodes
    ta.nodes[ta.node_index] = -ta.leaf_index
    ta.nodes[ta.node_index + 1] = -(ta.leaf_index + 1)

    ta.leaf_data[left_cats .+ 1, ta.leaf_index] .= true # +1 because index at 1
    ta.leaf_data[right_cats .+ 1, ta.leaf_index + 1] .= true # +1 because index at 1

    # Add everything to the stack
    push!(s, ta.node_index)
    push!(lb_s, lbs)
    push!(ub_s, left_ubs)

    push!(s, ta.node_index + 1)
    push!(lb_s, right_lbs)
    push!(ub_s, ubs)

    # Update indices
    ta.node_index += 2
    ta.leaf_index += 2
end