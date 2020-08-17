using DataStructures

"""
KDTrees
"""

mutable struct LEAFNODE
    cats::Vector{Int64}
    qvals::Vector{Float64}
end

mutable struct KDNODE
    dim::Int32
    split::Float64
    left::Union{KDNODE, LEAFNODE}
    right::Union{KDNODE, LEAFNODE}
end

function leafnode(;cats = collect(0:8), qvals = Vector{Float64}())
    return LEAFNODE(cats, qvals)
end

function kdnode(;dim = 1, split = 0.0, left = leafnode(), right = leafnode())
    return KDNODE(dim, split, left, right)
end

function get_leaf(root_node, point)
    curr_node = root_node
    while typeof(curr_node) != LEAFNODE
        val = point[curr_node.dim]
        split = curr_node.split
        curr_node = val < split ? curr_node.left : curr_node.right
    end
    return curr_node
end

# Note: this only works for 2D right now
function get_leaf_and_endpoints(root_node, point)
    bounds = [-0.5 0.5; -0.5 0.5]
    curr_node = root_node
    while typeof(curr_node) != LEAFNODE
        val = point[curr_node.dim]
        split = curr_node.split
        if val < split
            bounds[curr_node.dim, 2] = split
            curr_node = curr_node.left
        else
            bounds[curr_node.dim, 1] = split
            curr_node = curr_node.right
        end
    end
    endpoints = []
    for i = 1:2
        for j = 1:2
            push!(endpoints, [bounds[1,i], bounds[2,j]])
        end
    end
    return curr_node, endpoints
end

function get_bounds_and_cats(root_node)
    cats = []
    lbs = []
    ubs = []
    
    lb_s = Stack{Vector{Float64}}()
    ub_s = Stack{Vector{Float64}}()
    s = Stack{Union{LEAFNODE, KDNODE}}()

    push!(lb_s, [-0.5, -0.5])
    push!(ub_s, [0.5, 0.5])
    push!(s, root_node)

    while !isempty(s)
        curr = pop!(s)
        curr_lbs = pop!(lb_s)
        curr_ubs = pop!(ub_s)

        if typeof(curr) == LEAFNODE
            push!(cats, curr.cats)
            push!(lbs, curr_lbs)
            push!(ubs, curr_ubs)
        else
            # Traverse tree and keep track of bounds
            dim = curr.dim
            split = curr.split
            # Go left, upper bounds will change
            left_ubs = copy(curr_ubs)
            left_ubs[dim] = split

            push!(lb_s, curr_lbs)
            push!(ub_s, left_ubs)
            push!(s, curr.left)

            # Go right, lower bounds will change
            right_lbs = copy(curr_lbs)
            right_lbs[dim] = split
            
            push!(lb_s, right_lbs)
            push!(ub_s, curr_ubs)
            push!(s, curr.right)
        end
    end

    return lbs, ubs, cats
end

function get_bounds(root_node, point)
    lbs = [-0.5, -0.5]
    ubs = [0.5, 0.5]
    curr_node = root_node
    while typeof(curr_node) != LEAFNODE
        val = point[curr_node.dim]
        split = curr_node.split
        if val < split
            ubs[curr_node.dim] = split
            curr_node = curr_node.left
        else
            lbs[curr_node.dim] = split
            curr_node = curr_node.right
        end
    end
    return lbs, ubs
end

function get_overlapping_nodes(root_node, lbs, ubs)
    curr_node_list = Vector{LEAFNODE}()
    overlapping_nodes_helper(root_node, lbs, ubs, curr_node_list)
    return curr_node_list
end

function overlapping_nodes_helper(curr_node, lbs, ubs, curr_node_list)
    # Base case
    if typeof(curr_node) == LEAFNODE
        push!(curr_node_list, curr_node)
    else
        # Check if can prune either side
        dim = curr_node.dim
        split = curr_node.split
        if split > ubs[dim]
            # Can prune the right half
            overlapping_nodes_helper(curr_node.left, lbs, ubs, curr_node_list)
        elseif split < lbs[dim]
            # Can prune the left half
            overlapping_nodes_helper(curr_node.right, lbs, ubs, curr_node_list)
        else
            overlapping_nodes_helper(curr_node.left, lbs, ubs, curr_node_list)
            overlapping_nodes_helper(curr_node.right, lbs, ubs, curr_node_list)
        end
    end
end

""" Probabilities """

function get_max_prob(tree)
    # Base case
    if typeof(tree) == LEAFNODE
        return maximum(tree.qvals)
    else
        return max(get_max_prob(tree.left), get_max_prob(tree.right))
    end
end

function get_max_prob(tree, ḣ₀min, ḣ₀max)
    lbs = [-0.5, -0.5]
    ubs = [0.5, 0.5]
    return get_max_prob_helper(tree, lbs, ubs, ḣ₀min, ḣ₀max)
end

function get_max_prob_helper(tree, lbs, ubs, ḣ₀min, ḣ₀max)
    if typeof(tree) == LEAFNODE
        allowed = lbs[2] * 200.0 ≥ ḣ₀min && ubs[2] * 200.0 ≤ ḣ₀max
        val = allowed ? maximum(tree.qvals) : 0.0
        return lbs, ubs, val
    else
        # Traverse tree and keep track of bounds
        dim = tree.dim
        split = tree.split
        # Go left, upper bounds will change
        left_ubs = copy(ubs)
        left_ubs[dim] = split
        left_max_lbs, left_max_ubs, left_max = get_max_prob_helper(tree.left, lbs, left_ubs, ḣ₀min, ḣ₀max)
        # Go right, lower bounds will change
        right_lbs = copy(lbs)
        right_lbs[dim] = split
        right_max_lbs, right_max_ubs, right_max = get_max_prob_helper(tree.right, right_lbs, ubs, ḣ₀min, ḣ₀max)
        if left_max > right_max
            return left_max_lbs, left_max_ubs, left_max
        else
            return right_max_lbs, right_max_ubs, right_max
        end
    end
end

function get_max_prob_and_bounds(tree)
    if typeof(tree) == LEAFNODE
        return lbs, ubs, maximum(tree.qvals)
    else
        # Traverse tree and keep track of bounds
        dim = tree.dim
        split = tree.split
        # Go left, upper bounds will change
        left_ubs = copy(ubs)
        left_ubs[dim] = split
        left_max_lbs, left_max_ubs, left_max = get_max_prob_and_bounds_helper(tree.left, lbs, left_ubs)
        # Go right, lower bounds will change
        right_lbs = copy(lbs)
        right_lbs[dim] = split
        right_max_lbs, right_max_ubs, right_max = get_max_prob_and_bounds_helper(tree.right, right_lbs, ubs)
        if left_max > right_max
            return left_max_lbs, left_max_ubs, left_max
        else
            return right_max_lbs, right_max_ubs, right_max
        end
    end
end

function get_max_prob_and_bounds_helper(tree, lbs, ubs)
    # Base case
    if typeof(tree) == LEAFNODE
        return lbs, ubs, maximum(tree.qvals)
    else
        # Traverse tree and keep track of bounds
        dim = tree.dim
        split = tree.split
        # Go left, upper bounds will change
        left_ubs = copy(ubs)
        left_ubs[dim] = split
        left_max_lbs, left_max_ubs, left_max = get_max_prob_and_bounds_helper(tree.left, lbs, left_ubs)
        # Go right, lower bounds will change
        right_lbs = copy(lbs)
        right_lbs[dim] = split
        right_max_lbs, right_max_ubs, right_max = get_max_prob_and_bounds_helper(tree.right, right_lbs, ubs)
        if left_max > right_max
            return left_max_lbs, left_max_ubs, left_max
        else
            return right_max_lbs, right_max_ubs, right_max
        end
    end
end

function unnormalize_point(point::Vector{Float64})
    return point .* [16000.0, 200.0]
end

function normalize_point(point::Vector{Float64})
    return point ./ [16000.0, 200.0]
end