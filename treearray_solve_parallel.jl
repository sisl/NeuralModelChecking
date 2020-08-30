function check_tas_parallel(tas; post = "/scratch/smkatz/post2d/")
    updated_tas = Dict()
    tas_rc = Dict()
    # Do the tau = 0 trees for each pra
    println("τ: 0")
    for pra in actions
        ta = tas[(pra, 0.0)]
        post_prefix = "$(post)pra$(pra+1)tau0"
        tas_rc[pra] = remotecall(label_collisions, mod(pra + 1, nprocs() - 1) + 2, ta, post_prefix)
    end
    for pra in actions
        updated_tas[(pra, 0.0)] = fetch(tas_rc[pra])
    end
    # for the rest of the taus
    for τ = 1:40
        println("τ: $τ")
        τint = convert(Int64, τ)
        # do update for each pra
        for pra in actions
            ta = tas[(pra, τ)]
            post_prefix = "$(post)pra$(pra+1)tau$(τint)"
            tas_rc[pra] = remotecall(update, mod(pra + 1, nprocs() - 1) + 2, ta, updated_tas, pra, τ, post_prefix)
        end
        for pra in actions
            updated_tas[(pra, τ)] = fetch(tas_rc[pra])
        end
    end
end

function check_and_split_tas_parallel(tas; qthreshold = 0.5, threshold = 0.9,
    post = "/scratch/smkatz/post2d/", prefix = "/scratch/smkatz/VerticalCAS/networks/bugfix_pra0")
    
    updated_tas = Dict()
    tas_rc = Dict()
    # Do the tau = 0 trees for each pra
    println("τ: 0")
    for pra in actions
        ta = tas[(pra, 0.0)]
        post_prefix = "$(post)pra$(pra+1)tau0"
        tas_rc[pra] = remotecall(label_collisions, mod(pra + 1, nprocs() - 1) + 2, ta, post_prefix)
    end
    for pra in actions
        updated_tas[(pra, 0.0)] = fetch(tas_rc[pra])
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
            tas_rc[pra] = remotecall(update_or_split, mod(pra + 1, nprocs() - 1) + 2, 
                                        ta, updated_tas, pra, τ, post_prefix,
                                        qthreshold, threshold, prefix)
        end
        for pra in actions
            updated_tas[(pra, τ)] = fetch(tas_rc[pra])
        end
        println("Time: $(time() - start)")
    end
end