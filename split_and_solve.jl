function check_and_split(kdtrees; qthreshold = 0.5, threshold = 0.9, both = false, 
    prefix = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra0")
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
            update_or_split(tree, kdtrees, pra, τ, qthreshold = qthreshold, threshold = threshold, both = both)
        end
    end
end

function update_or_split(tree, kdtrees, pra, τ; qthreshold = 0.5, threshold = 0.9, both = false,
    prefix = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra0")
    lbs = [-0.5, -0.5]
    ubs = [0.5, 0.5]
    update_or_split_helper(tree, kdtrees, pra, τ, lbs, ubs, 
                             qthreshold = qthreshold, threshold = threshold, both = both, prefix = prefix)
end

function update_or_split_helper(tree, kdtrees, pra, τ, lbs, ubs; qthreshold = 0.5, threshold = 0.9, both = false,
    prefix = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra0")
    traverse_left = true
    traverse_right = true
    
    dim = tree.dim
    split = tree.split
    # Go left, upper bounds will change
    left_ubs = copy(ubs)
    left_ubs[dim] = split
    # Go right, lower bounds will change
    right_lbs = copy(lbs)
    right_lbs[dim] = split

    # Base case
    if typeof(tree.left) == LEAFNODE
        traverse_left = false
        # Figure out if need to split
        update_leaf!(tree, :left, kdtrees, pra, τ, lbs, left_ubs, 
                         qthreshold = qthreshold, threshold = threshold, both = both, prefix = prefix)
    end

    # Base case
    if typeof(tree.right) == LEAFNODE
        traverse_right = false
        # Figure out if need to split
        update_leaf!(tree, :right, kdtrees, pra, τ, right_lbs, ubs, 
                         qthreshold = qthreshold, threshold = threshold, both = both, prefix = prefix)
    end

    if traverse_left
        update_or_split_helper(tree.left, kdtrees, pra, τ, lbs, left_ubs, 
                                 qthreshold = qthreshold, threshold = threshold, prefix = prefix)
    end
    if traverse_right
        update_or_split_helper(tree.right, kdtrees, pra, τ, right_lbs, ubs, 
                                 qthreshold = qthreshold, threshold = threshold, prefix = prefix)
    end
end

function update_leaf!(tree, side, kdtrees, pra, τ, lbs, ubs; qthreshold = 0.5, threshold = 0.9, both = false,
    prefix = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra0")
    
    leaf = side == :left ? tree.left : tree.right
    # Figure out if need to split
    for i = 1:length(leaf.cats)
        next_tree = kdtrees[(leaf.cats[i], τ - 1)]
        nodes, probs, worst_nodes, worst_probs = transition_worst(next_tree, lbs, ubs, pra)
        next_vals = [maximum(nodes[i].qvals) for i = 1:length(nodes)]
        diff = maximum(next_vals) - minimum(next_vals)
        should_split = (diff > threshold) && (ubs[1] - lbs[1] > 10 / 16000)
        if should_split
            both ? split_both!(tree, lbs, ubs, side, prefix = prefix) : 
                   split!(tree, lbs, ubs, pra, τ, 1, side, prefix = prefix)
            if side == :left
                update_or_split_helper(tree.left, kdtrees, pra, τ, lbs, ubs, qthreshold = qthreshold, threshold = threshold, both = both)
            else
                update_or_split_helper(tree.right, kdtrees, pra, τ, lbs, ubs, qthreshold = qthreshold, threshold = threshold, both = both)
            end
            break
        else
            # Initialize the q values
            leaf.qvals = zeros(length(leaf.cats))
            for j = 1:length(worst_nodes)
                leaf.qvals[i] += worst_probs[j]*maximum(worst_nodes[j].qvals)
            end

            # Check if need to split due to actions
            qval_diff = maximum(leaf.qvals) - minimum(leaf.qvals)
            split_h = (qval_diff > qthreshold) && (ubs[1] - lbs[1] > 10 / 16000)
            split_ḣ₀ = (qval_diff > qthreshold) && (ubs[2] - lbs[2] > 1 / 200)
            if split_h && split_ḣ₀
                split!(tree, lbs, ubs, pra, τ, 1, side, prefix = prefix)

                split = (lbs[1] + ubs[1]) ./ 2
                # Go left, upper bounds will change
                left_ubs = copy(ubs)
                left_ubs[1] = split
                # Go right, lower bounds will change
                right_lbs = copy(lbs)
                right_lbs[1] = split

                if side == :left
                    split!(tree.left, lbs, left_ubs, pra, τ, 2, :left, prefix = prefix)
                    split!(tree.left, lbs, left_ubs, pra, τ, 2, :right, prefix = prefix)
                    update_or_split_helper(tree.left, kdtrees, pra, τ, lbs, ubs, qthreshold = qthreshold, threshold = threshold, both = both)
                else
                    split!(tree.right, right_lbs, ubs, pra, τ, 2, :left, prefix = prefix)
                    split!(tree.right, right_lbs, ubs, pra, τ, 2, :right, prefix = prefix)
                    update_or_split_helper(tree.right, kdtrees, pra, τ, lbs, ubs, qthreshold = qthreshold, threshold = threshold, both = both)
                end
                break
            elseif split_h
                split!(tree, lbs, ubs, pra, τ, 1, side, prefix = prefix)
                if side == :left
                    update_or_split_helper(tree.left, kdtrees, pra, τ, lbs, ubs, qthreshold = qthreshold, threshold = threshold, both = both)
                else
                    update_or_split_helper(tree.right, kdtrees, pra, τ, lbs, ubs, qthreshold = qthreshold, threshold = threshold, both = both)
                end
                break
            elseif split_ḣ₀
                split!(tree, lbs, ubs, pra, τ, 2, side, prefix = prefix)
                if side == :left
                    update_or_split_helper(tree.left, kdtrees, pra, τ, lbs, ubs, qthreshold = qthreshold, threshold = threshold, both = both)
                else
                    update_or_split_helper(tree.right, kdtrees, pra, τ, lbs, ubs, qthreshold = qthreshold, threshold = threshold, both = both)
                end
                break
            end
        end
    end
end

function split!(tree, lbs, ubs, pra, τ, dim, side;
     prefix = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra0")
    
    # side == :left ? println("before: $(typeof(tree.left))") : println("before: $(typeof(tree.right))")
    
    split = (lbs[dim] + ubs[dim]) ./ 2

    cats = side == :left ? tree.left.cats : tree.right.cats

    left_cats = cats
    right_cats = cats

    # New leaves
    if length(cats) > 1
        # Go left, upper bounds will change
        left_ubs = copy(ubs)
        left_ubs[dim] = split
        # Go right, lower bounds will change
        right_lbs = copy(lbs)
        right_lbs[dim] = split

        left_cats = get_categories!(cats, lbs, left_ubs, pra, τ, prefix = prefix)
        if length(left_cats) == 0
            left_cats = get_categories!(actions, lbs, left_ubs, pra, τ, prefix = prefix)
        end
        right_cats = get_categories!(cats, right_lbs, ubs, pra, τ, prefix = prefix)
        if length(right_cats) == 0
            right_cats = get_categories!(actions, right_lbs, ubs, pra, τ, prefix = prefix)
        end
    end

    if side == :left
        tree.left = kdnode(dim = dim, split = split)
        tree.left.left = leafnode(cats = left_cats)
        tree.left.right = leafnode(cats = right_cats)
    else
        tree.right = kdnode(dim = dim, split = split)
        tree.right.left = leafnode(cats = left_cats)
        tree.right.right = leafnode(cats = right_cats)
    end
    # side == :left ? println("after: $(typeof(tree.left))") : println("after: $(typeof(tree.right))")
end

function split_both!(tree, lbs, ubs, side;
    prefix = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra0")
    if side == :left
        cats = tree.left.cats
        split1 = (lbs[1] + ubs[1]) ./ 2
        split2 = (lbs[2] + ubs[2]) ./ 2

        tree.left = kdnode(dim = 1, split = split1)

        tree.left.left = kdnode(dim = 2, split = split2)
        tree.left.right = kdnode(dim = 2, split = split2)

        tree.left.left.left = leafnode(cats = cats)
        tree.left.left.right = leafnode(cats = cats)
        tree.left.right.left = leafnode(cats = cats)
        tree.left.right.right = leafnode(cats = cats)
    else
        cats = tree.right.cats
        split1 = (lbs[1] + ubs[1]) ./ 2
        split2 = (lbs[2] + ubs[2]) ./ 2

        tree.right = kdnode(dim = 1, split = split1)

        tree.right.left = kdnode(dim = 2, split = split2)
        tree.right.right = kdnode(dim = 2, split = split2)

        tree.right.left.left = leafnode(cats = cats)
        tree.right.left.right = leafnode(cats = cats)
        tree.right.right.left = leafnode(cats = cats)
        tree.right.right.right = leafnode(cats = cats)
    end
end