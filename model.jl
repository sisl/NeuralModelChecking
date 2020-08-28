"""
Constants
"""

ḣ₁ = -90
τs = collect(0:40)

"""
Advisories
"""

COC = 0
DNC = 1
DND = 2
DES1500 = 3
CL1500 = 4
SDES1500 = 5
SCL1500 = 6
SDES2500 = 7
SCL2500 = 8

actions = [COC, DNC, DND, DES1500, CL1500, SDES1500, SCL1500, SDES2500, SCL2500]

"""
Dictionaries to define transitions
"""

# Tuple of probabilities and corresponding acceleration in ft/s^2
probs = [0.5,0.3,0.2]
g = 32.2
accels = Dict(COC      => ([0.34,0.33,0.33], [0.0,-g/3,g/3]),
              DNC      =>(probs, [-g/3,-g/2,g/3]),
              DND      => (probs, [g/3,g/2,-g/3]),
              DES1500  => (probs, [-g/3,-g/2,g/3]),
              CL1500   => (probs, [g/3,g/2,-g/3]),
              SDES1500 => (probs, [-g/2.5,-g/2,g/3]),
              SCL1500  => (probs, [g/2.5,g/2,-g/3]),
              SDES2500 => (probs, [-g/2.5,-g/2,g/3]),
              SCL2500  => (probs, [g/2.5,g/2,-g/3]),
              -1       => ([0.34,0.33,0.33], [0,-g/8,g/8]))

# Velocity range where aircraft is NON-compliant with advisory (ft/s)
velRanges = Dict(COC      => (-Inf,Inf),
                 DNC      => (0.0,Inf),
                 DND      => (-Inf,0.0),
                 DES1500  => (-25.0,Inf),
                 CL1500   => (-Inf,25.0),
                 SDES1500 => (-25.0,Inf),
                 SCL1500  => (-Inf,25.0),
                 SDES2500 => (-41.67,Inf),
                 SCL2500  => (-Inf,41.67))

# Bounds should be normalized
function transition(next_tree, lbs, ubs, pra; worst_case = false)
    nodes = Vector{LEAFNODE}()
    probs = Vector{Float64}()

    ownProbs, ownAccels = accels[pra]

    for i = 1:length(ownAccels)
        next_lbs, next_ubs = dynamics(lbs, ubs, ownAccels[i], pra)

        overlapping_nodes = get_overlapping_nodes(next_tree, next_lbs, next_ubs)
        if worst_case
            worst_node = overlapping_nodes[1]
            worst_val = maximum(overlapping_nodes[1].qvals)
            for j = 2:length(overlapping_nodes)
                val = maximum(overlapping_nodes[j].qvals)
                if val > worst_val
                    worst_node = overlapping_nodes[j]
                    worst_val = val
                end
            end
            push!(nodes, worst_node)
            push!(probs, ownProbs[i])
        else
            append!(nodes, overlapping_nodes)

            num_nodes = length(overlapping_nodes)
            node_probs = ones(num_nodes) .* ownProbs[i] ./ num_nodes
            append!(probs, node_probs)
        end
    end
    return nodes, probs
end

function transition_worst(next_tree, lbs, ubs, pra)
    nodes = Vector{LEAFNODE}()
    probs = Vector{Float64}()
    worst_nodes = Vector{LEAFNODE}()
    worst_probs = Vector{Float64}()

    ownProbs, ownAccels = accels[pra]

    for i = 1:length(ownAccels)
        next_lbs, next_ubs = dynamics(lbs, ubs, ownAccels[i], pra)

        overlapping_nodes = get_overlapping_nodes(next_tree, next_lbs, next_ubs)

        worst_node = overlapping_nodes[1]
        worst_val = maximum(overlapping_nodes[1].qvals)
        for j = 2:length(overlapping_nodes)
            val = maximum(overlapping_nodes[j].qvals)
            if val > worst_val
                worst_node = overlapping_nodes[j]
                worst_val = val
            end
        end
        push!(worst_nodes, worst_node)
        push!(worst_probs, ownProbs[i])

        append!(nodes, overlapping_nodes)

        num_nodes = length(overlapping_nodes)
        node_probs = ones(num_nodes) .* ownProbs[i] ./ num_nodes
        append!(probs, node_probs)
    end
    return nodes, probs, worst_nodes, worst_probs
end

function transition_range(next_tree, lbs, ubs, pra)
    worst_nodes = Vector{LEAFNODE}()
    worst_probs = Vector{Float64}()

    max_range = 0.0

    ownProbs, ownAccels = accels[pra]

    for i = 1:length(ownAccels)
        next_lbs, next_ubs = dynamics(lbs, ubs, ownAccels[i], pra)

        overlapping_nodes = get_overlapping_nodes(next_tree, next_lbs, next_ubs)
        vals = [maximum(overlapping_nodes[i].qvals) for i = 1:length(overlapping_nodes)]

        worst_node = overlapping_nodes[1]
        worst_val = maximum(overlapping_nodes[1].qvals)
        for j = 2:length(overlapping_nodes)
            if vals[j] > worst_val
                worst_node = overlapping_nodes[j]
                worst_val = vals[j]
            end
        end
        push!(worst_nodes, worst_node)
        push!(worst_probs, ownProbs[i])

        range = maximum(vals) - minimum(vals)
        max_range = range > max_range ? range : max_range
    end
    return max_range, worst_nodes, worst_probs    
end

function transition_ta(ta::TREEARRAY, lbs, ubs, pra)
    worst_vals = Vector{Float64}()
    worst_probs = Vector{Float64}()

    max_range = 0.0

    ownProbs, ownAccels = accels[pra]

    for i = 1:length(ownAccels)
        next_lbs, next_ubs = dynamics(lbs, ubs, ownAccels[i], pra)

        overlapping_data_inds = get_overlapping_nodes(ta, next_lbs, next_ubs)
        vals = [maximum(ta.qvals[:, overlapping_data_inds[i]]) for i = 1:length(overlapping_data_inds)]

        push!(worst_vals, maximum(vals))
        push!(worst_probs, ownProbs[i])

        range = maximum(vals) - minimum(vals)
        max_range = range > max_range ? range : max_range
    end
    return max_range, worst_vals, worst_probs 
end

# Bounds should be normalized
function dynamics(lbs, ubs, ḧ₀, a)
    true_lbs = unnormalize_point(lbs)
    true_ubs = unnormalize_point(ubs)

    min_h = true_lbs[1]
    max_h = true_ubs[1]
    min_ḣ₀ = true_lbs[2]
    max_ḣ₀ = true_ubs[2]

    vLow, vHigh = velRanges[a]
    if (vLow ≥ max_ḣ₀) .| (vHigh ≤ min_ḣ₀)
        ḧ₀ = 0
    elseif vLow > max_ḣ₀ + ḧ₀
        ḧ₀ = vLow - max_ḣ₀
    elseif vHigh < min_ḣ₀ + ḧ₀
        ḧ₀ = vHigh - min_ḣ₀
    end

    min_next_ḣ₀ = min_ḣ₀ + ḧ₀
    max_next_ḣ₀ = max_ḣ₀ + ḧ₀

    min_next_h = min_h - max_ḣ₀ - 0.5*ḧ₀ + ḣ₁
    max_next_h = max_h - min_ḣ₀ - 0.5*ḧ₀ + ḣ₁

    min_point = normalize_point([min_next_h, min_next_ḣ₀])
    max_point = normalize_point([max_next_h, max_next_ḣ₀])
    return min_point, max_point
end