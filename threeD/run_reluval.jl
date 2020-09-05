function verify_region(lbs, ubs; advs = collect(0:8), network_path = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra01_v5_25HU_1000.nnet")
    poss_advs = []
    advs_string = "$(advs[1])"
    nadvs = length(advs)
    for i = 2:nadvs
        advs_string = "$advs_string $(advs[i])"
    end
    cmd = `./network_test_VC $(network_path) $(lbs[1]) $(lbs[2]) $(lbs[3]) $(lbs[4]) $(ubs[1]) $(ubs[2]) $(ubs[3]) $(ubs[4]) $advs`
    output = read(cmd, String)
    vals = parse.(Float64, split(output, ",")[1:end-1])
    times = vals[2:2:2*nadvs]
    advs = vals[1:2:2*nadvs-1]
    return times, advs
end

function get_categories!(cats, lbs, ubs, pra, τ; 
        prefix = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra0")
    
    # + 1 because networks labeled from 1
    network_path = "$(prefix)$(pra + 1)_v5_25HU_1000.nnet" 
    lbs, ubs = unnormalize_point_3d(lbs), unnormalize_point_3d(ubs)
    lbs = vcat(lbs, [τ])
    ubs = vcat(ubs, [τ])
    old_cats = deepcopy(cats)
    times, poss_cats = verify_region(lbs, ubs, advs = old_cats, network_path = network_path)
    poss_inds = findall(poss_cats .> 0)
    return old_cats[poss_inds]
end