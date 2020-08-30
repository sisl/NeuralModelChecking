using DataStructures
using Mmap

"""
Array Trees

Interior nodes:
nodes (Int64): left index (right index is left index + 1 always!)
dims (BitArray): one hot encoding of dimension
splits (Float64): split values

Leaf nodes:
nodes (Int64): -leaf_data index (know it is leaf node if this is a negative value)
dims (BitArray): N/A will just keep as garbage value or zero
splits (Float64): N/A will just keep as garbage value or zero
leaf_data (BitArray): one hot encoding of possible advisories
"""
mutable struct TREEARRAY
    nodes::Vector{Int32}
    splits::Vector{Float64}
    dims::BitArray
    leaf_data::BitArray
    qvals::Array{Float64, 2}
    node_index::Int32
    leaf_index::Int32
end

# CREATION AND MEMORY MANAGEMENT

function treearray(nodes::Vector{Int32},
                   splits::Vector{Float64},
                   dims::BitArray,
                   leaf_data::BitArray,
                   node_index::Int32,
                   leaf_index::Int32)    
    return TREEARRAY(nodes, splits, dims, leaf_data, 
                     Array{Float64, 2}(undef, 1, 1), node_index, leaf_index)
end

function treearray(totsize::Int64)
    return TREEARRAY(zeros(Int32, totsize), zeros(Float64, totsize),
                           falses(2, totsize), falses(9, totsize), zeros(Float64, 9, totsize),
                           1, 1)
end

function treearray_copy(ta::TREEARRAY)
    return TREEARRAY(copy(ta.nodes), copy(ta.splits), copy(ta.dims), copy(ta.leaf_data),
                zeros(9, size(ta.leaf_data, 2)), ta.node_index, ta.leaf_index)
end

function treearray_copy_and_extend(ta::TREEARRAY, totsize)
    updated_ta = TREEARRAY(zeros(Int32, totsize), zeros(Float64, totsize),
                           falses(2, totsize), falses(9, totsize), zeros(Float64, 9, totsize),
                           length(ta.nodes) + 1, size(ta.leaf_data, 2) + 1)

    # Fill in the data
    node_index = length(ta.nodes)
    leaf_index = size(ta.leaf_data, 2)

    updated_ta.nodes[1:node_index] = ta.nodes
    updated_ta.splits[1:node_index] = ta.splits
    updated_ta.dims[:, 1:node_index] = ta.dims
    updated_ta.leaf_data[:, 1:leaf_index] = ta.leaf_data

    return updated_ta
end

function sync_treearray(ta::TREEARRAY)
    Mmap.sync!(ta.nodes)
    Mmap.sync!(ta.splits)
    Mmap.sync!(ta.dims)
    Mmap.sync!(ta.leaf_data)
    Mmap.sync!(ta.qvals)
end

function mmap_trees(filepath)
    tas = Dict()
    for pra in 1:9
        for τ = 0.0:40.0
            τint = convert(Int64, τ)
            prefix = "$(filepath)pra$(pra)tau$(τint)"
            tas[(pra - 1, τ)] = mmap_tree(prefix)
        end
    end
    return tas
end

function mmap_tree(prefix)
    # Nodes file
    n = open("$(prefix)_n.bin")
    a = read(n, Int)
    node_index = convert(Int32, a + 1)
    nodes = Mmap.mmap(n, Vector{Int32}, a)

    # Splits file
    s = open("$(prefix)_s.bin")
    a = read(s, Int)
    splits = Mmap.mmap(s, Vector{Float64}, a)

    # Dims file
    d = open("$(prefix)_d.bin")
    a = read(d, Int)
    b = read(d, Int)
    dims = Mmap.mmap(d, BitArray, (a,b))

    # Leaf_data file
    l = open("$(prefix)_l.bin")
    a = read(l, Int)
    b = read(l, Int)
    leaf_index = convert(Int32, b + 1)
    leaf_data = Mmap.mmap(l, BitArray, (a,b))

    if isfile("$(prefix)_q.bin")
        # Qvals file
        q = open("$(prefix)_q.bin")
        a = read(q, Int)
        b = read(q, Int)
        qvals = Mmap.mmap(q, Matrix{Float64}, (a,b))
        return TREEARRAY(nodes, splits, dims, leaf_data, qvals, node_index, leaf_index)
    end

    return treearray(nodes, splits, dims, leaf_data, node_index, leaf_index)
end

function mmap_tree_modify(prefix)
    # Nodes file
    n = open("$(prefix)_n.bin", "r+")
    a = read(n, Int)
    node_index = convert(Int32, a + 1)
    nodes = Mmap.mmap(n, Vector{Int32}, a)

    # Splits file
    s = open("$(prefix)_s.bin", "r+")
    a = read(s, Int)
    splits = Mmap.mmap(s, Vector{Float64}, a)

    # Dims file
    d = open("$(prefix)_d.bin", "r+")
    a = read(d, Int)
    b = read(d, Int)
    dims = Mmap.mmap(d, BitArray, (a,b))

    # Leaf_data file
    l = open("$(prefix)_l.bin", "r+")
    a = read(l, Int)
    b = read(l, Int)
    leaf_index = convert(Int32, b + 1)
    leaf_data = Mmap.mmap(l, BitArray, (a,b))

    if isfile("$(prefix)_q.bin")
        # Qvals file
        q = open("$(prefix)_q.bin", "r+")
        a = read(q, Int)
        b = read(q, Int)
        qvals = Mmap.mmap(q, Matrix{Float64}, (a,b))
        return TREEARRAY(nodes, splits, dims, leaf_data, qvals, node_index, leaf_index)
    end

    return treearray(nodes, splits, dims, leaf_data, node_index, leaf_index)
end

function write_to_files(ta::TREEARRAY, prefix)
    # Nodes file
    n = open("$(prefix)_n.bin", "w+")
    write(n, length(ta.nodes))
    write(n, ta.nodes)
    close(n)

    # Splits file
    s = open("$(prefix)_s.bin", "w+")
    write(s, length(ta.splits))
    write(s, ta.splits)
    close(s)

    # Dims file
    d = open("$(prefix)_d.bin", "w+")
    write(d, size(ta.dims,1))
    write(d, size(ta.dims,2))
    write(d, ta.dims)
    close(d)

    # Leaf_data file
    l = open("$(prefix)_l.bin", "w+")
    write(l, size(ta.leaf_data,1))
    write(l, size(ta.leaf_data,2))
    write(l, ta.leaf_data)
    close(l)

    # Qvals file
    q = open("$(prefix)_q.bin", "w+")
    write(q, size(ta.qvals,1))
    write(q, size(ta.qvals,2))
    write(q, ta.qvals)
    close(q)
end

function write_to_files_and_truncate(ta::TREEARRAY, prefix)
    node_end = ta.node_index - 1
    leaf_end = ta.leaf_index - 1
    
    # Nodes file
    n = open("$(prefix)_n.bin", "w+")
    write(n, node_end)
    write(n, ta.nodes[1:node_end])
    close(n)

    # Splits file
    s = open("$(prefix)_s.bin", "w+")
    write(s, node_end)
    write(s, ta.splits[1:node_end])
    close(s)

    # Dims file
    d = open("$(prefix)_d.bin", "w+")
    write(d, size(ta.dims,1))
    write(d, node_end)
    write(d, ta.dims[:,1:node_end])
    close(d)

    # Leaf_data file
    l = open("$(prefix)_l.bin", "w+")
    write(l, size(ta.leaf_data,1))
    write(l, leaf_end)
    write(l, ta.leaf_data[:,1:leaf_end])
    close(l)

    # Qvals file
    q = open("$(prefix)_q.bin", "w+")
    write(q, size(ta.qvals,1))
    write(q, leaf_end)
    write(q, ta.qvals[:,1:leaf_end])
    close(q)
end


# TRAVERSAL AND EVALUATION

function get_bounds_and_cats(ta::TREEARRAY)
    nodes = ta.nodes
    splits = ta.splits
    dims = ta.dims
    leaf_data = ta.leaf_data
    
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
            one_hot = leaf_data[:, -nodes[curr]]
            push!(cats, findall(one_hot))
            push!(lbs, curr_lbs)
            push!(ubs, curr_ubs)
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

    return lbs, ubs, cats
end

# Will return positive indices of the corresponding leaf data
function get_overlapping_nodes(ta::TREEARRAY, lbs, ubs)
    nodes = ta.nodes
    splits = ta.splits
    dims = ta.dims
    leaf_data = ta.leaf_data
    
    overlapping_node_data_inds = []
    
    s = Stack{Int32}()

    push!(s, 1)

    while !isempty(s)
        curr = pop!(s)

        if nodes[curr] < 0 # leaf node
            push!(overlapping_node_data_inds, -nodes[curr])
        else # interior node
            # Check if can prune either side
            dim = findfirst(dims[:, curr])
            split = splits[curr]
            if split > ubs[dim]
                # Can prune the right half
                push!(s, nodes[curr])
            elseif split < lbs[dim]
                # Can prune the left half
                push!(s, nodes[curr] + 1)
            else
                push!(s, nodes[curr])
                push!(s, nodes[curr] + 1)
            end
        end
    end

    return overlapping_node_data_inds
end

function get_qvals(ta::TREEARRAY, point)
    curr = 1 # root of tree
    while ta.nodes[curr] > 0 # not at leaf
        split = ta.splits[curr]
        dim = findfirst(ta.dims[:, curr])

        # Go left or right
        curr = point[dim] < split ? ta.nodes[curr] : ta.nodes[curr] + 1
    end
    return ta.qvals[:, -ta.nodes[curr]]
end

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