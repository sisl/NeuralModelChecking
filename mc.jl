using Distributions

function load_nnets(;network_path = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra0")
    nnets = Dict()
    for pra in actions
        nnets[pra] = read_network("$(network_path)$(pra+1)_v5_25HU_1000.nnet")
    end
    return nnets
end

function mc_sim(s₀, num_sims, nnets)
    nmacs = 0
    for i = 1:num_sims
        s = rollout(s₀, nnets)
        abs(s[1]) < 100 ? nmacs += 1 : nothing
    end
    return nmacs/num_sims
end

function rollout(s₀, nnets)
    s = deepcopy(s₀)
    for τ = 40:-1:1
        qvals = evaluate_network(nnets[s[4]], [s; τ])
        a = argmax(qvals) - 1
        next_s = mc_transition(s, a)
        s = deepcopy(next_s)
    end
    return s
end

function mc_transition(s, a)
    pra = s[4]
    probs₀, accels₀ = accels[pra]
    dist₀ = Categorical(probs₀)
    ḧ₀ = accels₀[rand(dist₀)]

    next_h, next_ḣ₀ = mc_dynamics(s, pra, ḧ₀)

    return [next_h, next_ḣ₀, ḣ₁, a]
end

function mc_dynamics(s, a, ḧ₀)
	h, ḣ₀, ḣ₁ = s[1], s[2], s[3]
    
    vLow, vHigh = velRanges[a]
    if (vLow ≥ ḣ₀) .| (vHigh ≤ ḣ₀)
        ḧ₀ = 0
    elseif vLow > ḣ₀ + ḧ₀
        ḧ₀ = vLow - ḣ₀
    elseif vHigh < ḣ₀ + ḧ₀
        ḧ₀ = vHigh - ḣ₀
    end

    next_h = h - ḣ₀ - 0.5*ḧ₀ + ḣ₁
    next_ḣ₀ = ḣ₀ + ḧ₀

    return next_h, clamp(next_ḣ₀, -100, 100), ḣ₁
end