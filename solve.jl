function check(kdtrees; worst_case = false)
    # Do the tau = 0 trees for each pra
    println("τ: 0")
    for pra in actions
        tree = kdtrees[(pra, 0.0)]
        label_collisions(tree)
    end
    # for the rest of the taus
    for τ = 1:40
        println("τ: $τ")
        # do update for each pra
        for pra in actions
            tree = kdtrees[(pra, τ)]
            update(tree, kdtrees, pra, τ, worst_case = worst_case)
        end
    end
end

# NOTE: currently not working
function check_parallel(kdtrees)
    # Do the tau = 0 trees for each pra
    println("τ: 0")
    for pra in actions
        tree = kdtrees[(pra, 0.0)]
        remotecall(label_collisions, mod(pra, nprocs()-1)+2, tree)
    end
    # for the rest of the taus
    for τ = 1:40
        println("τ: $τ")
        # do update for each pra
        for pra in actions
            tree = kdtrees[(pra, τ)]
            remotecall(update, mod(pra, nprocs()-1)+2, tree)
        end
    end
end

function label_collisions(tree)
    lbs = [-0.5, -0.5]
    ubs = [0.5, 0.5]
    label_collisions_helper(tree, lbs, ubs)
end

function label_collisions_helper(tree, lbs, ubs)
    # Base case
    if typeof(tree) == LEAFNODE
        # Label with probability 1 if it indicates a collision, 0 else
        collision_h = 100 ./ 16000
        collision = abs(lbs[1]) < collision_h || abs(ubs[1]) < collision_h
        if collision
            #println("lbs: $lbs, ubs: $ubs")
        end
        nvals = length(tree.cats)
        qvals = collision ? ones(Float64, nvals) : zeros(Float64, nvals)
        tree.qvals = qvals
    else
        # Traverse tree and keep track of bounds
        dim = tree.dim
        split = tree.split
        # Go left, upper bounds will change
        left_ubs = copy(ubs)
        left_ubs[dim] = split
        label_collisions_helper(tree.left, lbs, left_ubs)
        # Go right, lower bounds will change
        right_lbs = copy(lbs)
        right_lbs[dim] = split
        label_collisions_helper(tree.right, right_lbs, ubs)
    end
end

function update(tree, kdtrees, pra, τ; worst_case = false)
    lbs = [-0.5, -0.5]
    ubs = [0.5, 0.5]
    update_helper(tree, kdtrees, pra, τ, lbs, ubs, worst_case = worst_case)
end

function update_helper(tree, kdtrees, pra, τ, lbs, ubs; worst_case = false)
    # Base case
    if typeof(tree) == LEAFNODE
        # Initialize the q values
        tree.qvals = zeros(length(tree.cats))

        # Iterate over possible actions
        for i = 1:length(tree.cats)
            next_tree = kdtrees[(tree.cats[i], τ - 1)]
            nodes, probs = transition(next_tree, lbs, ubs, pra, worst_case = worst_case)
            for j = 1:length(nodes)
                tree.qvals[i] += probs[j]*maximum(nodes[j].qvals)
            end
        end
    else
        # Traverse tree and keep track of bounds
        dim = tree.dim
        split = tree.split
        # Go left, upper bounds will change
        left_ubs = copy(ubs)
        left_ubs[dim] = split
        update_helper(tree.left, kdtrees, pra, τ, lbs, left_ubs, worst_case = worst_case)
        # Go right, lower bounds will change
        right_lbs = copy(lbs)
        right_lbs[dim] = split
        update_helper(tree.right, kdtrees, pra, τ, right_lbs, ubs, worst_case = worst_case)
    end
end