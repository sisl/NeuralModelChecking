using PGFPlots
using Interact
using Colors
using ColorBrewer
using DataStructures

ḣ₀s = vcat(LinRange(-100,-60,5),
           LinRange(-50,-35,4),
           LinRange(-30,30,21),
           LinRange(35,50,4),
           LinRange(60,100,5))

"""
Tree Arrays
"""
function plot_nadvs_3d(sta::SHARED_TREEARRAY, ḣ₁)
    lbs_u, ubs_u, cats = get_bounds_and_cats(sta)

    # Unnormalize everying
    lbs = [unnormalize_point_3d(lbs_u[i]) for i = 1:length(lbs_u)]
    ubs = [unnormalize_point_3d(ubs_u[i]) for i = 1:length(ubs_u)]

    ymin = -8000
    ymax = 8000
    xmin = -100
    xmax = 100

    ax = Axis(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, width="10cm", height="10cm", 
    xlabel=L"$h$", ylabel=L"$\dot{h}_0$", title="Number of Possible Advisories")

    for i = 1:length(lbs)
        # Determine if it contains ḣ₁
        if lbs[i][3] < ḣ₁ < ubs[i][3]
            if length(cats[i]) == 1
                color = "blue"
            elseif length(cats[i]) == 2
                color = "red"
            else
                color = "yellow"
            end
            push!(ax, Plots.Command(get_filled_rectangle([lbs[i][2], lbs[i][1]],
                                                        [ubs[i][2], ubs[i][1]], color)))
        end
    end

    return ax
end

function viz_prob_shared(stas)
    currSavePlot = 0
    first_call = true

    @manipulate for fileName in textbox(value="myFile.pdf",label="File Name") |> onchange,
        savePlot in button("Save Plot"), 
        vl in pad(0.3em, nothing),
		refreshPlot in button("Refresh"),
    	nbin = 100,
        xmin = 0.0,
        xmax = 40.0,
        ymin = -2000.0,
        ymax = 2000.0,
        logprob = [true, false],
        pra = collect(0:8),
        ḣ₀ = ḣ₀s

        function get_heat(x, y)
            npoint = normalize_point([y, ḣ₀])
            τ = round(x)
            sta = stas[(pra, τ)]
            qvals = get_qvals(sta, npoint)
            if logprob
                return -log(maximum(qvals) + 1e-10)
            else
                return -maximum(qvals)
            end
        end

        ax = Axis([
            Plots.Image(get_heat, (xmin, xmax), (ymin, ymax),
            xbins = nbin, ybins = nbin, colormap = pasteljet, colorbar=false),
            ], xmin=xmin, xmax=xmax, ymin=ymin,ymax=ymax, width="10cm", height="8cm", 
               xlabel="Tau (s)", ylabel="Relative Alt (ft)", title="Probabilities")

        if savePlot > currSavePlot
        	ax2 = Axis([
                Plots.Image(get_heat, (xmin, xmax), (ymin, ymax),
                xbins = nbin, ybins = nbin, colormap = pasteljet, colorbar=false),
                ], xmin=xmin, xmax=xmax, ymin=ymin,ymax=ymax, width="10cm", height="8cm", 
                   xlabel="Tau (s)", ylabel="Relative Alt (ft)", title="Probabilities")
            PGFPlots.save(fileName, ax2, include_preamble=:false)
            currSavePlot += 1
        end

        if first_call
			first_call = false
			PGFPlots.cleanup(ax)
			return nothing
		else
			return ax
		end
    end
end

"""
Kdtrees
"""

function viz_probability(kdtrees)
    currSavePlot = 0
    first_call = true

    @manipulate for fileName in textbox(value="myFile.pdf",label="File Name") |> onchange,
        savePlot in button("Save Plot"), 
        vl in pad(0.3em, nothing),
		refreshPlot in button("Refresh"),
    	nbin = 100,
        xmin = 0.0,
        xmax = 40.0,
        ymin = -2000.0,
        ymax = 2000.0,
        logprob = [true, false],
        pra = collect(0:8),
        ḣ₀ = ḣ₀s

        function get_heat(x, y)
            npoint = normalize_point([y, ḣ₀])
            τ = round(x)
            tree = kdtrees[(pra, τ)]
            leaf = get_leaf(tree, npoint)
            if logprob
                return -log(maximum(leaf.qvals) + 1e-10)
            else
                return -maximum(leaf.qvals)
            end
        end

        ax = Axis([
            Plots.Image(get_heat, (xmin, xmax), (ymin, ymax),
            xbins = nbin, ybins = nbin, colormap = pasteljet, colorbar=false),
            ], xmin=xmin, xmax=xmax, ymin=ymin,ymax=ymax, width="10cm", height="8cm", 
               xlabel="Tau (s)", ylabel="Relative Alt (ft)", title="Probabilities")

        if savePlot > currSavePlot
        	ax2 = Axis([
                Plots.Image(get_heat, (xmin, xmax), (ymin, ymax),
                xbins = nbin, ybins = nbin, colormap = pasteljet, colorbar=false),
                ], xmin=xmin, xmax=xmax, ymin=ymin,ymax=ymax, width="10cm", height="8cm", 
                   xlabel="Tau (s)", ylabel="Relative Alt (ft)", title="Probabilities")
            PGFPlots.save(fileName, ax2, include_preamble=:false)
            currSavePlot += 1
        end

        if first_call
			first_call = false
			PGFPlots.cleanup(ax)
			return nothing
		else
			return ax
		end
    end
end

function viz_probability_tree(kdtrees)
    currSavePlot = 0

    @manipulate for fileName in textbox(value="myFile.pdf",label="File Name") |> onchange,
		savePlot in button("Save Plot"), 
    	nbin = 100,
        xmin = -100.0,
        xmax = 100.0,
        ymin = -2000.0,
        ymax = 2000.0,
        logprob = [true, false],
        pra = collect(0:8),
        τ = collect(0.0:40.0)

        function get_heat(x, y)
            npoint = normalize_point([y, x])
            tree = kdtrees[(pra, τ)]
            leaf = get_leaf(tree, npoint)
            if logprob
                return -log(maximum(leaf.qvals) + 1e-10)
            else
                return -maximum(leaf.qvals)
            end
        end

        zmin = logprob ? 0.0 : -1.0
        zmax = logprob ? -log(1e-10) : 0.0

        ax = Axis([
            Plots.Image(get_heat, (xmin, xmax), (ymin, ymax), zmin = zmin, zmax = zmax,
            xbins = nbin, ybins = nbin, colormap = pasteljet, colorbar=false),
            ], xmin=xmin, xmax=xmax, ymin=ymin,ymax=ymax, width="10cm", height="8cm", 
               xlabel="Tau (s)", ylabel="Relative Alt (ft)", title="Probabilities")

        if savePlot > currSavePlot
        	ax2 = Axis([
                Plots.Image(get_heat, (xmin, xmax), (ymin, ymax), zmin = zmin, zmax = zmax,
                xbins = nbin, ybins = nbin, colormap = pasteljet, colorbar=false),
                ], xmin=xmin, xmax=xmax, ymin=ymin,ymax=ymax, width="10cm", height="8cm", 
                   xlabel="Tau (s)", ylabel="Relative Alt (ft)", title="Probabilities")
            PGFPlots.save(fileName, ax2, include_preamble=:false)
            currSavePlot += 1
        end

        return ax
    end
end

function viz_probability_comparison(left_kdtrees, right_kdtrees)
    currSavePlot = 0

    @manipulate for fileName in textbox(value="myFile.pdf",label="File Name") |> onchange,
		savePlot in button("Save Plot"), 
    	nbin = 100,
        xmin = 0.0,
        xmax = 40.0,
        ymin = -2000.0,
        ymax = 2000.0,
        logprob = [true, false],
        pra = collect(0:8),
        ḣ₀ = ḣ₀s

        function get_heat_left(x, y)
            npoint = normalize_point([y, ḣ₀])
            τ = round(x)
            tree = left_kdtrees[(pra, τ)]
            leaf = get_leaf(tree, npoint)
            if logprob
                return -log(maximum(leaf.qvals) + 1e-10)
            else
                return -maximum(leaf.qvals)
            end
        end

        function get_heat_right(x, y)
            npoint = normalize_point([y, ḣ₀])
            τ = round(x)
            tree = right_kdtrees[(pra, τ)]
            leaf = get_leaf(tree, npoint)
            if logprob
                return -log(maximum(leaf.qvals) + 1e-10)
            else
                return -maximum(leaf.qvals)
            end
        end

        zmin = logprob ? 0.0 : -1.0
        zmax = logprob ? -log(1e-10) : 0.0

        ax_left = Axis([
            Plots.Image(get_heat_left, (xmin, xmax), (ymin, ymax), zmin = zmin, zmax = zmax,
            xbins = nbin, ybins = nbin, colormap = pasteljet, colorbar=false),
            ], xmin=xmin, xmax=xmax, ymin=ymin,ymax=ymax, width="10cm", height="8cm", 
               xlabel="Tau (s)", ylabel="Relative Alt (ft)", title="Probabilities")

        ax_right = Axis([
            Plots.Image(get_heat_right, (xmin, xmax), (ymin, ymax), zmin = zmin, zmax = zmax,
            xbins = nbin, ybins = nbin, colormap = pasteljet, colorbar=false),
            ], xmin=xmin, xmax=xmax, ymin=ymin,ymax=ymax, width="10cm", height="8cm", 
               xlabel="Tau (s)", ylabel="Relative Alt (ft)", title="Probabilities")

        g = GroupPlot(2, 1, groupStyle = "horizontal sep=3cm")
        push!(g, ax_left)
        push!(g, ax_right)

        if savePlot > currSavePlot
        	ax2_left = Axis([
                Plots.Image(get_heat_left, (xmin, xmax), (ymin, ymax),  zmin = zmin, zmax = zmax,
                xbins = nbin, ybins = nbin, colormap = pasteljet, colorbar=false),
                ], xmin=xmin, xmax=xmax, ymin=ymin,ymax=ymax, width="10cm", height="8cm", 
                   xlabel="Tau (s)", ylabel="Relative Alt (ft)", title="Probabilities")

            ax2_right = Axis([
                Plots.Image(get_heat_right, (xmin, xmax), (ymin, ymax),  zmin = zmin, zmax = zmax,
                xbins = nbin, ybins = nbin, colormap = pasteljet, colorbar=false),
                ], xmin=xmin, xmax=xmax, ymin=ymin,ymax=ymax, width="10cm", height="8cm", 
                   xlabel="Tau (s)", ylabel="Relative Alt (ft)", title="Probabilities")

            g2 = GroupPlot(2, 1, groupStyle = "horizontal sep=3cm")
            push!(g2, ax_left)
            push!(g2, ax_right)

            PGFPlots.save(fileName, g2, include_preamble=:false)
            currSavePlot += 1
        end

        return g
    end
end

function plot_slice(tree, τ; logprob = false)
    nbin = 100
    xmin = -100.0
    xmax = 100.0
    ymin = -8000.0
    ymax = 8000.0
    
    function get_heat(x, y)
        npoint = normalize_point([y, x])
        leaf = get_leaf(tree, npoint)
        if logprob
            return -log(maximum(leaf.qvals) + 1e-10)
        else
            return -maximum(leaf.qvals)
        end
    end

    zmin = logprob ? 0.0 : -1.0
    zmax = logprob ? -log(1e-10) : 0.0

    ax = Axis([
            Plots.Image(get_heat, (xmin, xmax), (ymin, ymax), zmin = zmin, zmax = zmax,
            xbins = nbin, ybins = nbin, colormap = pasteljet, colorbar=false),
            ], xmin=xmin, xmax=xmax, ymin=ymin,ymax=ymax, width="10cm", height="8cm", 
               xlabel="Tau (s)", ylabel="Relative Alt (ft)", title="Tau = $τ")
    return ax
end

function plot_nadvs(tree)
    lbs_u, ubs_u, cats = get_bounds_and_cats(tree)

    # Unnormalize everying
    lbs = [unnormalize_point(lbs_u[i]) for i = 1:length(lbs_u)]
    ubs = [unnormalize_point(ubs_u[i]) for i = 1:length(ubs_u)]

    ymin = -8000
    ymax = 8000
    xmin = -100
    xmax = 100

    ax = Axis(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, width="10cm", height="10cm", 
    xlabel=L"$h$", ylabel=L"$\dot{h}_0$", title="Number of Possible Advisories")

    for i = 1:length(lbs)
        if length(cats[i]) == 1
            color = "blue"
        elseif length(cats[i]) == 2
            color = "red"
        else
            color = "yellow"
        end
        push!(ax, Plots.Command(get_filled_rectangle([lbs[i][2], lbs[i][1]],
                                                     [ubs[i][2], ubs[i][1]], color)))
    end

    return ax
end

function plot_regions(kdtrees, pra, ḣ₀)
    ymin = -8000
    ymax = 8000
    xmin = 0
    xmax = 41

    ax = Axis(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, width="10cm", height="10cm", 
    xlabel=L"$\tau$", ylabel=L"$h$", title="Regions")
    
    for τ = 0:40
        tree = kdtrees[(pra, τ)]
        lbs, ubs, _ = get_bounds_and_cats(tree)
        normalized_ḣ₀ = ḣ₀ / 200
        for i = 1:length(lbs)
            if lbs[i][2] ≤ normalized_ḣ₀ ≤ ubs[i][2]
                push!(ax, Plots.Command(get_rectangle([τ, lbs[i][1] * 16000], 
                                                      [τ+1, ubs[i][1] * 16000])))
            end
        end
    end

    return ax
end

function plot_nadvs_side(kdtrees, pra, ḣ₀; ymin = -8000, ymax = 8000)
    xmin = 0
    xmax = 41

    ax = Axis(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, width="10cm", height="10cm", 
    xlabel=L"$\tau$", ylabel=L"$h$", title="Regions")
    
    for τ = 0:40
        tree = kdtrees[(pra, τ)]
        lbs, ubs, cats = get_bounds_and_cats(tree)
        normalized_ḣ₀ = ḣ₀ / 200
        for i = 1:length(lbs)
            contains_ḣ₀ = lbs[i][2] ≤ normalized_ḣ₀ ≤ ubs[i][2]
            in_plot = lbs[i][1] ≥ ymin && ubs[i][1] ≤ ymax
            if contains_ḣ₀ && in_plot
                if length(cats[i]) == 1
                    color = "blue"
                elseif length(cats[i]) == 2
                    color = "red"
                else
                    color = "yellow"
                end
                push!(ax, Plots.Command(get_filled_rectangle([τ, lbs[i][1] * 16000], 
                                                      [τ+1, ubs[i][1] * 16000], color)))
            end
        end
    end

    return ax
end

function plot_nadvs_overlay(kdtrees1, kdtrees2, pra, ḣ₀, max_τ; 
                                                ymin = 0, ymax = 2000)
    
    xmin = 0
    xmax = 40

    ax = Axis(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, width="10cm", height="10cm", 
    xlabel=L"$\tau$", ylabel=L"$h$", title="Regions")

    # Rectangles 1
    for τ = 0:max_τ
        tree = kdtrees1[(pra, τ)]
        lbs, ubs, cats = get_bounds_and_cats(tree)
        normalized_ḣ₀ = ḣ₀ / 200
        for i = 1:length(lbs)
            contains_ḣ₀ = lbs[i][2] ≤ normalized_ḣ₀ ≤ ubs[i][2]
            in_plot = lbs[i][1] ≥ ymin && ubs[i][1] ≤ ymax
            if contains_ḣ₀ && in_plot
                if length(cats[i]) == 1
                    color = "blue"
                elseif length(cats[i]) == 2
                    color = "red"
                else
                    color = "yellow"
                end
                push!(ax, Plots.Command(get_filled_rectangle([max(0.0, τ-0.5), lbs[i][1] * 16000], 
                                                      [min(40.0, τ+0.5), ubs[i][1] * 16000],
                                                      color)))
            end
        end
    end

    # Rectangles 2
    for τ = max_τ+1:40
        tree = kdtrees2[(pra, τ)]
        lbs, ubs, cats = get_bounds_and_cats(tree)
        normalized_ḣ₀ = ḣ₀ / 200
        for i = 1:length(lbs)
            contains_ḣ₀ = lbs[i][2] ≤ normalized_ḣ₀ ≤ ubs[i][2]
            in_plot = lbs[i][1] ≥ ymin && ubs[i][1] ≤ ymax
            if contains_ḣ₀ && in_plot
                if length(cats[i]) == 1
                    color = "blue"
                elseif length(cats[i]) == 2
                    color = "red"
                else
                    color = "yellow"
                end
                push!(ax, Plots.Command(get_filled_rectangle([max(0.0, τ-0.5), lbs[i][1] * 16000], 
                                                      [min(40.0, τ+0.5), ubs[i][1] * 16000],
                                                      color)))
            end
        end
    end

    return ax
end

function plot_probability(kdtrees, pra, ḣ₀, max_τ; ymin = 0, ymax = 2000, logprob =false)
    xmin = 0
    xmax = 40

    ax = Axis(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, width="10cm", height="10cm", 
    xlabel=L"$\tau$", ylabel=L"$h$", title="Probability")

    # Probability
    nbin = 500

    function get_heat(x, y)
        npoint = normalize_point([y, ḣ₀])
        τ = round(x)
        tree = kdtrees[(pra, τ)]
        leaf = get_leaf(tree, npoint)
        if logprob
            return -log(maximum(leaf.qvals) + 1e-10)
        else
            return -maximum(leaf.qvals)
        end
    end

    push!(ax, Plots.Image(get_heat, (xmin, min(40.0, max_τ+0.5)), (ymin, ymax),
        xbins = nbin, ybins = nbin, colormap = pasteljet, colorbar=false))

    return ax
end

function plot_probability_with_regions(kdtrees, pra, ḣ₀, max_τ; logprob =false, plot_middle=false)
    ymin = -8000
    ymax = 8000
    xmin = 0
    xmax = 40

    ax = Axis(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, width="10cm", height="10cm", 
    xlabel=L"$\tau$", ylabel=L"$h$")

    # Probability
    nbin = 500

    function get_heat(x, y)
        npoint = normalize_point([y, ḣ₀])
        τ = round(x)
        tree = kdtrees[(pra, τ)]
        leaf = get_leaf(tree, npoint)
        if logprob
            return -log(maximum(leaf.qvals) + 1e-10)
        else
            return -maximum(leaf.qvals)
        end
    end

    push!(ax, Plots.Image(get_heat, (xmin, min(40.0, max_τ+0.5)), (ymin, ymax),
        xbins = nbin, ybins = nbin, colormap = pasteljet, colorbar=false))

    # Rectangles
    for τ = 0:max_τ
        tree = kdtrees[(pra, τ)]
        lbs, ubs, _ = get_bounds_and_cats(tree)
        normalized_ḣ₀ = ḣ₀ / 200
        for i = 1:length(lbs)
            should_plot = plot_middle || lbs[i][1] * 16000 > 1000 || (ubs[i][1] * 16000 < -1000)
            if (lbs[i][2] ≤ normalized_ḣ₀ ≤ ubs[i][2]) && should_plot
                push!(ax, Plots.Command(get_rectangle([max(0.0, τ-0.5), lbs[i][1] * 16000], 
                                                      [min(40.0, τ+0.5), ubs[i][1] * 16000],
                                                      color = "white", linewidth = 0)))
            end
        end
    end

    return ax
end

function plot_probability_with_regions_both(kdtrees1, kdtrees2, pra, ḣ₀, max_τ; logprob =false, plot_middle=false)
    ax1 = plot_probability_with_regions(kdtrees1, pra, ḣ₀, max_τ, logprob = logprob, plot_middle = plot_middle)
    ax2 = plot_probability_with_regions(kdtrees2, pra, ḣ₀, max_τ, logprob = logprob, plot_middle = plot_middle)
    g = GroupPlot(2, 1, groupStyle = "horizontal sep=3cm")
    push!(g, ax2)
    push!(g, ax1)
    return g
end

function plot_probability_with_regions_overlay(kdtrees1, kdtrees2, pra, ḣ₀, max_τ; 
                                                logprob =false, plot_middle = true)
    ymin = -8000
    ymax = 8000
    xmin = 0
    xmax = 40

    ax = Axis(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, width="15cm", height="15cm", 
    xlabel=L"$\tau$", ylabel=L"$h$", title="Regions")

    # Probability 1
    nbin = 100

    function get_heat1(x, y)
        npoint = normalize_point([y, ḣ₀])
        τ = round(x)
        tree = kdtrees1[(pra, τ)]
        leaf = get_leaf(tree, npoint)
        if logprob
            return -log(maximum(leaf.qvals) + 1e-10)
        else
            return -maximum(leaf.qvals)
        end
    end

    push!(ax, Plots.Image(get_heat1, (xmin, min(40.0, max_τ+0.5)), (ymin, ymax),
        xbins = nbin, ybins = nbin, colormap = pasteljet, colorbar=false))

    # Probability 2
    function get_heat2(x, y)
        npoint = normalize_point([y, ḣ₀])
        τ = round(x)
        tree = kdtrees2[(pra, τ)]
        leaf = get_leaf(tree, npoint)
        if logprob
            return -log(maximum(leaf.qvals) + 1e-10)
        else
            return -maximum(leaf.qvals)
        end
    end

    push!(ax, Plots.Image(get_heat2, (max(0.0, max_τ-0.5), 40.0), (ymin, ymax),
        xbins = nbin, ybins = nbin, colormap = pasteljet, colorbar=false))

    # Rectangles 1
    for τ = 0:max_τ
        tree = kdtrees1[(pra, τ)]
        lbs, ubs, _ = get_bounds_and_cats(tree)
        normalized_ḣ₀ = ḣ₀ / 200
        for i = 1:length(lbs)
            should_plot = plot_middle || lbs[i][1] * 16000 > 1000 || (ubs[i][1] * 16000 < -1000)
            if (lbs[i][2] ≤ normalized_ḣ₀ ≤ ubs[i][2]) && should_plot
                push!(ax, Plots.Command(get_rectangle([max(0.0, τ-0.5), lbs[i][1] * 16000], 
                                                      [min(40.0, τ+0.5), ubs[i][1] * 16000],
                                                      color = "white", linewidth = 0)))
            end
        end
    end

    # Rectangles 2
    for τ = max_τ+1:40
        tree = kdtrees2[(pra, τ)]
        lbs, ubs, _ = get_bounds_and_cats(tree)
        normalized_ḣ₀ = ḣ₀ / 200
        for i = 1:length(lbs)
            should_plot = plot_middle || lbs[i][1] * 16000 > 1000 || (ubs[i][1] * 16000 < -1000)
            if (lbs[i][2] ≤ normalized_ḣ₀ ≤ ubs[i][2]) && should_plot
                push!(ax, Plots.Command(get_rectangle([max(0.0, τ-0.5), lbs[i][1] * 16000], 
                                                      [min(40.0, τ+0.5), ubs[i][1] * 16000],
                                                      color = "white", linewidth = 0)))
            end
        end
    end

    return ax
end

function plot_probability_with_regions_overlay_reg(kdtrees1, kdtrees2, pra, ḣ₀, max_τ; 
                                                logprob =false, plot_middle = true)
    ymin = -8000
    ymax = 8000
    xmin = 0
    xmax = 40

    ax = Axis(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, width="10cm", height="10cm", 
    xlabel=L"$\tau$", ylabel=L"$h$")

    # Probability 1
    nbin = 500

    function get_heat1(x, y)
        npoint = normalize_point([y, ḣ₀])
        τ = round(x)
        tree = kdtrees1[(pra, τ)]
        leaf = get_leaf(tree, npoint)
        if logprob
            return -log(maximum(leaf.qvals) + 1e-10)
        else
            return -maximum(leaf.qvals)
        end
    end

    push!(ax, Plots.Image(get_heat1, (xmin, min(40.0, max_τ+0.5)), (ymin, ymax),
        xbins = nbin, ybins = nbin, colormap = pasteljet, colorbar=false))

    # Rectangles 1
    for τ = 0:max_τ
        tree = kdtrees1[(pra, τ)]
        lbs, ubs, _ = get_bounds_and_cats(tree)
        normalized_ḣ₀ = ḣ₀ / 200
        for i = 1:length(lbs)
            should_plot = plot_middle || lbs[i][1] * 16000 > 1000 || (ubs[i][1] * 16000 < -1000)
            if (lbs[i][2] ≤ normalized_ḣ₀ ≤ ubs[i][2]) && should_plot
                push!(ax, Plots.Command(get_rectangle([max(0.0, τ-0.5), lbs[i][1] * 16000], 
                                                      [min(40.0, τ+0.5), ubs[i][1] * 16000],
                                                      color = "white", linewidth = 0)))
            end
        end
    end

    # Rectangles 2
    for τ = max_τ+1:40
        tree = kdtrees2[(pra, τ)]
        lbs, ubs, _ = get_bounds_and_cats(tree)
        normalized_ḣ₀ = ḣ₀ / 200
        for i = 1:length(lbs)
            should_plot = plot_middle || lbs[i][1] * 16000 > 1000 || (ubs[i][1] * 16000 < -1000)
            if (lbs[i][2] ≤ normalized_ḣ₀ ≤ ubs[i][2]) && should_plot
                push!(ax, Plots.Command(get_rectangle([max(0.0, τ-0.5), lbs[i][1] * 16000], 
                                                      [min(40.0, τ+0.5), ubs[i][1] * 16000],
                                                      color = "black", linewidth = 0)))
            end
        end
    end

    return ax
end

function gif_view(kdtrees1, kdtrees2, pra, ḣ₀, max_τ; logprob =false, plot_middle=false)
    ax1 = plot_probability_with_regions_overlay_reg(kdtrees1, kdtrees2, pra, ḣ₀, max_τ, logprob = logprob, plot_middle = plot_middle)
    ax2 = plot_probability_with_regions_overlay_reg(kdtrees2, kdtrees2, pra, ḣ₀, max_τ, logprob = logprob, plot_middle = plot_middle)
    g = GroupPlot(2, 1, groupStyle = "horizontal sep=3cm")
    push!(g, ax2)
    push!(g, ax1)
    return g
end

function gif_view_nadvs(kdtrees1, kdtrees2, pra, ḣ₀, max_τ; ymin = 0, ymax = 2000, logprob =false)
    ax1 = plot_probability(kdtrees1, pra, ḣ₀, max_τ, ymin = ymin, ymax = ymax, logprob = logprob)
    ax2 = plot_probability(kdtrees2, pra, ḣ₀, max_τ, ymin = ymin, ymax = ymax, logprob = logprob)
    ax3 = plot_nadvs_overlay(kdtrees1, kdtrees2, pra, ḣ₀, max_τ, ymin = ymin, ymax = ymax)
    ax4 = plot_nadvs_side(kdtrees2, pra, ḣ₀, ymin = ymin, ymax = ymax)
    g = GroupPlot(2, 2, groupStyle = "horizontal sep=3cm, vertical sep=2cm")
    push!(g, ax2)
    push!(g, ax1)
    push!(g, ax4)
    push!(g, ax3)
    return g
end

function get_rectangle(lb, ub; color = "black", linewidth = 0.25)
    return "\\draw[$(color),line width=$(string(linewidth))mm] ($(string(lb[1])),$(string(lb[2]))) rectangle ($(string(ub[1])),$(string(ub[2])));"
end

function get_filled_rectangle(lb, ub, color)
    return "\\filldraw[fill=$(color), draw=black] ($(string(lb[1])),$(string(lb[2]))) rectangle ($(string(ub[1])),$(string(ub[2])));"
end