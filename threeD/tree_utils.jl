using DataStructures
using Mmap
using SharedArrays

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

mutable struct SHARED_TREEARRAY
    nodes::SharedArray{Int32, 1}
    splits::SharedArray{Float64, 1}
    dims::SharedArray{Bool, 2}
    leaf_data::SharedArray{Bool, 2}
    qvals::SharedArray{Float64, 2}
    node_index::Int32
    leaf_index::Int32
end

# SHARED ARRAYS

# CREATION AND MEMORY MANAGEMENT

function create_files_for_sa(old_prefix, new_prefix)
    n = open("$(old_prefix)_n.bin")
    a = read(n, Int)
    node_index = convert(Int32, a + 1)
    nodes = Mmap.mmap(n, Vector{Int32}, a)
    close(n)

    n_new = open("$(new_prefix)_n.bin", "w+")
    write(n_new, nodes)

    close(n)
    close(n_new)

    # Splits file
    s = open("$(old_prefix)_s.bin")
    a = read(s, Int)
    splits = Mmap.mmap(s, Vector{Float64}, a)

    s_new = open("$(new_prefix)_s.bin", "w+")
    write(s_new, splits)

    close(s)
    close(s_new)

    # Dims file
    d = open("$(old_prefix)_d.bin")
    a = read(d, Int)
    b = read(d, Int)
    dims = Mmap.mmap(d, BitArray, (a,b))
    dims_bool = convert(Array{Bool}, copy(dims))

    d_new = open("$(new_prefix)_d.bin", "w+")
    write(d_new, dims_bool)

    close(d)
    close(d_new)

    # Leaf_data file
    l = open("$(old_prefix)_l.bin")
    a = read(l, Int)
    b = read(l, Int)
    leaf_index = convert(Int32, b + 1)
    leaf_data = Mmap.mmap(l, BitArray, (a,b))
    leaf_data_bool = convert(Array{Bool}, copy(leaf_data))

    l_new = open("$(new_prefix)_l.bin", "w+")
    write(l_new, leaf_data_bool)

    close(l)
    close(l_new)

    i = open("$(new_prefix)_i.bin", "w+")
    write(i, node_index - 1)
    write(i, leaf_index - 1)

    close(i)
end

function truncate_files(old_prefix, new_prefix, totsize)
    n = open("$(old_prefix)_n.bin")
    nodes = Mmap.mmap(n, Vector{Int32}, totsize)
    node_index = findfirst(nodes .== 0) - 1
    leaf_index = convert(Int64, -minimum(nodes))
    close(n)

    n_new = open("$(new_prefix)_n.bin", "w+")
    write(n_new, nodes[1:node_index])

    close(n)
    close(n_new)

    # Splits file
    s = open("$(old_prefix)_s.bin")
    splits = Mmap.mmap(s, Vector{Float64}, totsize)

    s_new = open("$(new_prefix)_s.bin", "w+")
    write(s_new, splits[1:node_index])

    close(s)
    close(s_new)

    # Dims file
    d = open("$(old_prefix)_d.bin")
    dims = Mmap.mmap(d, Array{Bool, 2}, (2, totsize))

    d_new = open("$(new_prefix)_d.bin", "w+")
    write(d_new, dims[:, 1:node_index])

    close(d)
    close(d_new)

    # Leaf_data file
    l = open("$(old_prefix)_l.bin")
    leaf_data = Mmap.mmap(l, Array{Bool, 2}, (9, totsize))

    l_new = open("$(new_prefix)_l.bin", "w+")
    write(l_new, leaf_data[:, 1:leaf_index])

    close(l)
    close(l_new)

    # Qvals file
    q = open("$(old_prefix)_q.bin")
    qvals = Mmap.mmap(q, Array{Float64, 2}, (9, totsize))

    q_new = open("$(new_prefix)_q.bin", "w+")
    write(q_new, qvals[:, 1:leaf_index])

    close(q)
    close(q_new)

    i = open("$(new_prefix)_i.bin", "w+")
    write(i, node_index)
    write(i, leaf_index)

    close(i)
end

function shared_treearray(prefix)
    # Get the dimensions
    i = open("$(prefix)_i.bin")
    a = read(i, Int)
    b = read(i, Int)
    close(i)

    # Read everything into shared arrays
    nodes = SharedArray{Int32, 1}("$(prefix)_n.bin", (a,))
    splits = SharedArray{Float64, 1}("$(prefix)_s.bin", (a,))
    dims = SharedArray{Bool, 2}("$(prefix)_d.bin", (3, a))
    leaf_data = SharedArray{Bool, 2}("$(prefix)_l.bin", (9, b))
    if isfile("$(prefix)_q.bin")
        qvals = SharedArray{Float64, 2}("$(prefix)_q.bin", (9, b))
    else
        qvals = SharedArray{Float64, 2}((1,1))
    end

    node_index = convert(Int32, a + 1)
    leaf_index = convert(Int32, b + 1)

    return SHARED_TREEARRAY(nodes, splits, dims, leaf_data, qvals, node_index, leaf_index)
end

function shared_treearray(prefix, a, b)
    # Read everything into shared arrays
    nodes = SharedArray{Int32, 1}("$(prefix)_n.bin", (a,))
    splits = SharedArray{Float64, 1}("$(prefix)_s.bin", (a,))
    dims = SharedArray{Bool, 2}("$(prefix)_d.bin", (3, a))
    leaf_data = SharedArray{Bool, 2}("$(prefix)_l.bin", (9, b))
    if isfile("$(prefix)_q.bin")
        qvals = SharedArray{Float64, 2}("$(prefix)_q.bin", (9, b))
    else
        qvals = SharedArray{Float64, 2}((1,1))
    end

    node_index = convert(Int32, a + 1)
    leaf_index = convert(Int32, b + 1)

    return SHARED_TREEARRAY(nodes, splits, dims, leaf_data, qvals, node_index, leaf_index)
end

function shared_treearrays(filepath)
    stas = Dict()
    for pra in 1:9
        for τ = 0.0:40.0
            τint = convert(Int64, τ)
            prefix = "$(filepath)pra$(pra)tau$(τint)"
            stas[(pra - 1, τ)] = shared_treearray(prefix)
        end
    end
    return stas
end

function shared_treearrays(filepath, a, b)
    stas = Dict()
    for pra in 1:9
        for τ = 0.0:40.0
            τint = convert(Int64, τ)
            prefix = "$(filepath)pra$(pra)tau$(τint)"
            stas[(pra - 1, τ)] = shared_treearray(prefix, a, b)
        end
    end
    return stas
end

function shared_treearrays(filepath, a, b, taumax)
    stas = Dict()
    for pra in 1:9
        for τ = 0.0:taumax
            τint = convert(Int64, τ)
            prefix = "$(filepath)pra$(pra)tau$(τint)"
            stas[(pra - 1, τ)] = shared_treearray(prefix, a, b)
        end
    end
    return stas
end

function shared_treearrays(filepath, a, b, taumin, taumax)
    stas = Dict()
    for pra in 1:9
        for τ = taumin:taumax
            τint = convert(Int64, τ)
            prefix = "$(filepath)pra$(pra)tau$(τint)"
            stas[(pra - 1, τ)] = shared_treearray(prefix, a, b)
        end
    end
    return stas
end

# CAREFUL, ONLY LOAD A PREFIX YOU ARE COOL WITH MODIFYING THE FILES IN
function shared_treearray_copy(sta::SHARED_TREEARRAY, post_prefix)
    a = sta.node_index - 1
    b = sta.leaf_index - 1

    nodes = SharedArray{Int32, 1}("$(post_prefix)_n.bin", (a,), mode="w+")
    nodes[:] = sta.nodes
    splits = SharedArray{Float64, 1}("$(post_prefix)_s.bin", (a,), mode="w+")
    splits[:] = sta.splits
    dims = SharedArray{Bool, 2}("$(post_prefix)_d.bin", (3, a), mode="w+")
    dims[:,:] = sta.dims
    leaf_data = SharedArray{Bool, 2}("$(post_prefix)_l.bin", (9, b), mode="w+")
    leaf_data[:,:] = sta.leaf_data
    qvals = SharedArray{Float64, 2}("$(post_prefix)_q.bin", (9, b), mode="w+")
    qvals[:,:] = zeros(9,b)

    return SHARED_TREEARRAY(nodes, splits, dims, leaf_data, qvals, sta.node_index, sta.leaf_index)
end

# CAREFUL, ONLY LOAD A PREFIX YOU ARE COOL WITH MODIFYING THE FILES IN
function shared_treearray_copy_and_extend(sta::SHARED_TREEARRAY, totsize, post_prefix)
    a = sta.node_index - 1
    b = sta.leaf_index - 1

    nodes = SharedArray{Int32, 1}("$(post_prefix)_n.bin", (totsize,), mode="w+")
    nodes[1:a] = sta.nodes
    splits = SharedArray{Float64, 1}("$(post_prefix)_s.bin", (totsize,), mode="w+")
    splits[1:a] = sta.splits
    dims = SharedArray{Bool, 2}("$(post_prefix)_d.bin", (3, totsize), mode="w+")
    dims[:,1:a] = sta.dims
    leaf_data = SharedArray{Bool, 2}("$(post_prefix)_l.bin", (9, totsize), mode="w+")
    leaf_data[:,1:b] = sta.leaf_data
    qvals = SharedArray{Float64, 2}("$(post_prefix)_q.bin", (9, totsize), mode="w+")
    qvals[:,:] = zeros(9,totsize)

    return SHARED_TREEARRAY(nodes, splits, dims, leaf_data, qvals, sta.node_index, sta.leaf_index)
end

# TRAVERSAL AND EVALUATION

function get_bounds_and_cats(sta::SHARED_TREEARRAY)
    nodes = sta.nodes
    splits = sta.splits
    dims = sta.dims
    leaf_data = sta.leaf_data
    
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

function get_max_prob(sta::SHARED_TREEARRAY)
    nodes = sta.nodes
    splits = sta.splits
    dims = sta.dims
    qvals = sta.qvals

    s = Stack{Int32}()
    push!(s, 1)

    max_prob = 0.0

    while !isempty(s)
        curr = pop!(s)

        if nodes[curr] < 0 # leaf node
            prob = maximum(qvals[:, -nodes[curr]])
            max_prob = prob > max_prob ? prob : max_prob
        else # interior node
            # Traverse tree
            dim = findall(dims[:, curr])[1]
            split = splits[curr]

            # Go left
            push!(s, nodes[curr])
            # Go right
            push!(s, nodes[curr] + 1)
        end
    end
    return max_prob
end

# Will return positive indices of the corresponding leaf data
function get_overlapping_nodes(ta::Union{TREEARRAY, SHARED_TREEARRAY}, lbs, ubs)
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

function get_qvals(ta::Union{TREEARRAY, SHARED_TREEARRAY}, point)
    curr = 1 # root of tree
    while ta.nodes[curr] > 0 # not at leaf
        split = ta.splits[curr]
        dim = findfirst(ta.dims[:, curr])

        # Go left or right
        curr = point[dim] < split ? ta.nodes[curr] : ta.nodes[curr] + 1
    end
    return ta.qvals[:, -ta.nodes[curr]]
end

# MMAP

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

function mmapped_treearray(totsize::Int64, prefix)
    # Nodes file
    n = open("$(prefix)_n.bin", "w+")
    write(n, totsize)
    nodes = Mmap.mmap(n, Vector{Int32}, totsize, shared = true)
    close(n)

    # Splits file
    s = open("$(prefix)_s.bin", "w+")
    write(s, totsize)
    splits = Mmap.mmap(s, Vector{Int32}, totsize, shared = true)
    close(s)

    # Dims file
    d = open("$(prefix)_d.bin", "w+")
    write(d, 2)
    write(d, totsize)
    dims = Mmap.mmap(d, BitArray, (2, totsize), shared = true)
    close(d)

    # Leaf_data file
    l = open("$(prefix)_l.bin", "w+")
    write(l, 9)
    write(l, totsize)
    leaf_data = Mmap.mmap(l, BitArray, (9, totsize), shared = true)
    close(l)

    # Qvals file
    q = open("$(prefix)_q.bin", "w+")
    write(q, 9)
    write(q, totsize)
    qvals = Mmap.mmap(q, Array{Float64,2}, (9, totsize), shared = true)
    close(q)

    return mmap_tree_modify(prefix)
end

function copy_over!(from::TREEARRAY, to::TREEARRAY)
    node_index = length(from.nodes)
    leaf_index = size(from.leaf_data, 2)

    to.nodes = zeros(Int32, length(to.nodes))
    to.nodes[1:node_index] = from.nodes

    to.splits = zeros(Float64, length(to.splits))
    to.splits[1:node_index] = from.splits

    to.dims = falses(size(to.dims,1), size(to.dims,2))
    to.dims[:, 1:node_index] = from.dims

    to.leaf_data = falses(size(to.leaf_data,1), size(to.leaf_data,2))
    to.leaf_data[:, 1:leaf_index] = from.leaf_data

    to.node_index = length(from.nodes) + 1
    to.leaf_index = size(from.leaf_data, 2) + 1
end

function treearray_copy_and_map(ta::TREEARRAY, totsize, prefix)
    updated_ta = mmapped_treearray(totsize, prefix)

    # Fill in the data
    node_index = length(ta.nodes)
    leaf_index = size(ta.leaf_data, 2)

    updated_ta.nodes[1:node_index] = ta.nodes
    updated_ta.splits[1:node_index] = ta.splits
    updated_ta.dims[:, 1:node_index] = ta.dims
    updated_ta.leaf_data[:, 1:leaf_index] = ta.leaf_data

    updated_ta.node_index = node_index + 1
    updated_ta.leaf_index = leaf_index + 1

    sync_treearray(updated_ta)
    return updated_ta
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
    nodes = Mmap.mmap(n, Vector{Int32}, a, shared = true)

    # Splits file
    s = open("$(prefix)_s.bin")
    a = read(s, Int)
    splits = Mmap.mmap(s, Vector{Float64}, a, shared = true)

    # Dims file
    d = open("$(prefix)_d.bin")
    a = read(d, Int)
    b = read(d, Int)
    dims = Mmap.mmap(d, BitArray, (a,b), shared = true)

    # Leaf_data file
    l = open("$(prefix)_l.bin")
    a = read(l, Int)
    b = read(l, Int)
    leaf_index = convert(Int32, b + 1)
    leaf_data = Mmap.mmap(l, BitArray, (a,b), shared = true)

    if isfile("$(prefix)_q.bin")
        # Qvals file
        q = open("$(prefix)_q.bin")
        a = read(q, Int)
        b = read(q, Int)
        qvals = Mmap.mmap(q, Matrix{Float64}, (a,b), shared = true)
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

function unnormalize_point_3d(point::Vector{Float64})
    return point .* [16000.0, 200.0, 200.0]
end

function normalize_point_3d(point::Vector{Float64})
    return point ./ [16000.0, 200.0, 200.0]
end