function check_stas(stas; post = "/scratch/smkatz/post3dshared/")
    updated_stas = Dict()
    tas_rc = Dict()
    # Do the tau = 0 trees for each pra
    println("τ: 0")
    for pra in actions
        sta = stas[(pra, 0.0)]
        post_prefix = "$(post)pra$(pra+1)tau0"
        updated_stas[(pra, 0.0)] = shared_treearray_copy(sta, post_prefix)
    end

    @sync @distributed for pra in actions
        label_collisions_distrib(updated_stas[(pra, 0.0)])
    end
    
    for τ = 1:40
        println("τ: $τ")
        τint = convert(Int64, τ)
        
        for pra in actions
            sta = stas[(pra, τ)]
            post_prefix = "$(post)pra$(pra+1)tau$(τint)"
            updated_stas[(pra, τ)] = shared_treearray_copy_and_extend(sta, 80000, post_prefix)
        end

        @sync @distributed for pra in actions
            update_distrib(updated_stas[(pra, τ)], updated_stas, pra, τ)
        end
    end

    return updated_stas
end

function label_collisions_distrib(updated_sta::SHARED_TREEARRAY)
    collision_h = 100 / 16000
    
    nodes = updated_sta.nodes
    splits = updated_sta.splits
    dims = updated_sta.dims
    leaf_data = updated_sta.leaf_data
    qvals = updated_sta.qvals
    
    cats = []
    lbs = []
    ubs = []
    
    lb_s = Stack{Vector{Float64}}()
    ub_s = Stack{Vector{Float64}}()
    s = Stack{Int32}()

    push!(lb_s, [-0.5, -0.5, -0.5])
    push!(ub_s, [0.5, 0.5, 0.5])
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

    return updated_sta
end

function update_distrib(updated_sta::SHARED_TREEARRAY, updated_stas, pra, τ)   
    
    nodes = updated_sta.nodes
    splits = updated_sta.splits
    dims = updated_sta.dims
    leaf_data = updated_sta.leaf_data
    qvals = updated_sta.qvals
    
    cats = []
    lbs = []
    ubs = []
    
    lb_s = Stack{Vector{Float64}}()
    ub_s = Stack{Vector{Float64}}()
    s = Stack{Int32}()

    push!(lb_s, [-0.5, -0.5, -0.5])
    push!(ub_s, [0.5, 0.5, 0.5])
    push!(s, 1)

    while !isempty(s)
        curr = pop!(s)
        curr_lbs = pop!(lb_s)
        curr_ubs = pop!(ub_s)

        if nodes[curr] < 0 # leaf node
            # Get the possible actions
            cats = findall(leaf_data[:, -nodes[curr]]) .- 1
            for i = 1:length(cats)
                next_ta = updated_stas[cats[i], τ - 1]
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
    
    return updated_sta
end

function check_and_split_stas(stas; post = "/scratch/smkatz/post3dshared/", qthreshold = 0.5, threshold = 0.9,
    prefix = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra0")

    updated_stas = Dict()
    tas_rc = Dict()
    # Do the tau = 0 trees for each pra
    println("τ: 0")
    for pra in actions
        sta = stas[(pra, 0.0)]
        post_prefix = "$(post)pra$(pra+1)tau0"
        updated_stas[(pra, 0.0)] = shared_treearray_copy(sta, post_prefix)
    end

    @sync @distributed for pra in actions
        label_collisions_distrib(updated_stas[(pra, 0.0)])
    end
    
    # for τ = 1:1
    #     println("τ: $τ")
    #     start = time()

    #     τint = convert(Int64, τ)
        
    #     for pra in actions
    #         sta = stas[(pra, τ)]
    #         post_prefix = "$(post)pra$(pra+1)tau$(τint)"
    #         updated_stas[(pra, τ)] = shared_treearray_copy_and_extend(sta, 3000000, post_prefix)
    #     end

    #     @sync @distributed for pra in actions
    #         update_or_split_distrib(updated_stas[(pra, τ)], updated_stas, pra, τ, 
    #                     qthreshold = qthreshold, threshold = threshold, prefix = prefix)
    #     end
    #     println("Time required: $(time() - start)")
    # end

    return updated_stas
end

function check_and_split_stas(stas, updated_stas, taumin, taumax; post = "/scratch/smkatz/post3dshared/", qthreshold = 0.5, threshold = 0.9,
    prefix = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra0")

    tas_rc = Dict()
    
    for τ = taumin:taumax
        println("τ: $τ")
        start = time()

        τint = convert(Int64, τ)
        
        for pra in actions
            sta = stas[(pra, τ)]
            post_prefix = "$(post)pra$(pra+1)tau$(τint)"
            updated_stas[(pra, τ)] = shared_treearray_copy_and_extend(sta, 7000000, post_prefix)
        end

        @sync @distributed for pra in actions
            update_or_split_distrib(updated_stas[(pra, τ)], updated_stas, pra, τ, 
                        qthreshold = qthreshold, threshold = threshold, prefix = prefix)
        end
        println("Time required: $(time() - start)")
    end

    return updated_stas
end

function update_or_split_distrib(updated_sta::SHARED_TREEARRAY, updated_stas, pra, τ; qthreshold = 0.5, threshold = 0.9,
    prefix = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra0")
    
    cats = []
    lbs = []
    ubs = []
    
    lb_s = Stack{Vector{Float64}}()
    ub_s = Stack{Vector{Float64}}()
    s = Stack{Int32}()

    push!(lb_s, [-0.5, -0.5, -0.5])
    push!(ub_s, [0.5, 0.5, 0.5])
    push!(s, 1)

    while !isempty(s)
        curr = pop!(s)
        curr_lbs = pop!(lb_s)
        curr_ubs = pop!(ub_s)

        if updated_sta.nodes[curr] < 0 # leaf node
            # Get the possible actions
            update_leaf!(updated_sta, updated_stas, pra, τ, curr, curr_lbs, curr_ubs, s, lb_s, ub_s, 
                            qthreshold = qthreshold, threshold = threshold, prefix = prefix)
        else # interior node
            # Traverse tree and keep track of bounds
            dim = findall(updated_sta.dims[:, curr])[1]
            split = updated_sta.splits[curr]
            # Go left, upper bounds will change
            left_ubs = copy(curr_ubs)
            left_ubs[dim] = split

            push!(lb_s, curr_lbs)
            push!(ub_s, left_ubs)
            push!(s, updated_sta.nodes[curr])

            # Go right, lower bounds will change
            right_lbs = copy(curr_lbs)
            right_lbs[dim] = split
            
            push!(lb_s, right_lbs)
            push!(ub_s, curr_ubs)
            push!(s, updated_sta.nodes[curr] + 1)
        end
    end

    return updated_sta
end