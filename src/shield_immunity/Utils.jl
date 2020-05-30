function population_wide_states(simulation; scaling=1.0)
    ((pop, _), solution) = simulation
    solution_m = hcat(solution.u...) .* scaling 
    
    states = zeros(pop.n_states, length(solution.t))
    
    state_idxs = idxs_per_state(pop)    
    for (i, idxs)=enumerate(state_idxs)
        states[i, :] = sum(solution_m[idxs, :], dims=1)
    end
    
    # Deaths per day.
    D_per_day = [states[8, 2:end] - states[8, 1:(end-1)]; 0] 
    
    # Deaths per class.
    D_per_class = solution_m[state_idxs[end], end]
    
    # Hospitalised.
    JJ = findmax(states[6, :])[2]    
    Ihs = solution_m[state_idxs[5], JJ]
    Ihc = solution_m[state_idxs[6], JJ]
    day_peak = solution.t[JJ]
    ICU = sum(states[6, JJ])
    
    return states, D_per_day, D_per_class, Ihs, Ihc, day_peak, ICU
end;

function plot_dynamics(simulations; scaling=1, bw=false, kwargs...)
    
    low_results = map(x->population_wide_states(x; scaling=scaling), simulations["low"])
    high_results = map(x->population_wide_states(x; scaling=scaling), simulations["high"])
    γ_low, γ_high = simulations["low"][1][1][6:7], simulations["high"][1][1][6:7]

    # Peak ranges.
    peak_l_min, peak_l_max = Int(ceil(minimum(map(x->x[6], low_results)))), Int(ceil(maximum(map(x->x[6], low_results))))
    peak_h_min, peak_h_max = Int(ceil(minimum(map(x->x[6], high_results)))), Int(ceil(maximum(map(x->x[6], high_results))))
    
    # ICU ranges.
    ICU_l_min, ICU_l_max = Int(ceil(minimum(map(x->x[7], low_results)))), Int(ceil(maximum(map(x->x[7], low_results))))
    ICU_h_min, ICU_h_max = Int(ceil(minimum(map(x->x[7], high_results)))), Int(ceil(maximum(map(x->x[7], high_results))))
    
    classes = 1:length(low_results[1][3])
    
    if bw
        color_list = reshape(range(colorant"black", stop=colorant"gray",length=length(low_results)), 1, length(low_results))
        line_kwargs = Dict(:lc=>color_list)
        point_kwargs = Dict(:seriescolor=>color_list, :markershape=>:circle, :markerstrokecolor=>:white, :xticks=>classes)
    else
        line_kwargs = Dict(:lc=>:auto)
        point_kwargs = Dict(:seriescolor=>:auto, :markershape=>:circle, :markerstrokecolor=>:white, :xticks=>classes)
    end
    
    # Ploting Section.
    ## Titles...
    main_title = plot(title = "COVID-19 Epidemic\n Scale per $scaling", grid = false, showaxis = false, bottom_margin = -50Plots.px)
    low_title = "Low R0\nDay of peak: $peak_l_max - $peak_l_min\n ICU: $ICU_l_max - $ICU_l_min"
    high_title = "High R0\n Day of peak: $peak_h_max - $peak_h_min\n ICU: $ICU_h_max - $ICU_h_min"
    
    ## Low plots.
    pl1 = plot(reduce(hcat, map(x->x[2], low_results)), ylabel="Deaths per day", title=low_title; line_kwargs...)
    pl2 = plot(reduce(hcat, map(x->x[1][6, :], low_results)), xlabel="Time (days)", ylabel="ICU bed demand"; line_kwargs...)
    pl3 = plot(reduce(hcat, map(x->x[3], low_results)), ylabel="Cumulative Deaths"; point_kwargs...)
    pl4 = plot(reduce(hcat, map(x->x[4], low_results)), ylabel="Hospitilised subacute"; point_kwargs...)
    pl5 = plot(reduce(hcat, map(x->x[5], low_results)), xlabel="Class", ylabel="Hospitilised ICU"; point_kwargs...)
    
    ## High plots.
    ph1 = plot(reduce(hcat, map(x->x[2], high_results)), title=high_title; line_kwargs...)
    ph2 = plot(reduce(hcat, map(x->x[1][6, :], high_results)), xlabel="Time (days)"; line_kwargs...)
    ph3 = plot(reduce(hcat, map(x->x[3], high_results)); point_kwargs...)
    ph4 = plot(reduce(hcat, map(x->x[4], high_results)); point_kwargs...)
    ph5 = plot(reduce(hcat, map(x->x[5], high_results)), xlabel="Class"; point_kwargs...)

    ll = @layout[title{0.05h}; A1 A2; B1 B2; C1 C2; D1 D2; E1 E2]
    figure = plot(main_title, pl1, ph1, pl2, ph2, pl3, ph3, pl4, ph4, pl5, ph5, layout=ll, legend=false; kwargs...)
    return figure
end;