function SerologicalShields(du, u, par, t)
    pop, β_a, β_s, p, γ_e, γ_a, γ_s, γ_h, h, ξ, μ , α, sero_idxs = par
    S_idx, E_idx, Ia_idx, Is_idx, Ihsub_idx, Ihcrit_idx, R_idx, D_indx = idxs_per_state(pop)
    
    # Calculate totals.
    total_Ia = sum(u[Ia_idx]) # Asymptomatic infectious individuals.
    total_Is = sum(u[Is_idx]) # Symptomatic infectious individuals.
    total_shieldR = sum(u[R_idx][sero_idxs]) # Classes as potential serological shields.
    
    not_D_set = Set(reduce(vcat, [Ihsub_idx, Ihcrit_idx, D_indx]))
    total_not_D = sum(u[collect(setdiff(Set(idxs(pop)), not_D_set))]) # All individuals not in the hospitalised (Ihsub,Ihcrit) or dead (D) states.
    
    for (i, idxs) in enumerate(idxs_per_class(pop))
        S, E, Ia, Is, Ihsub, Ihcrit, R, D = u[idxs]
        
        du[idxs[1]] = dS = -β_a*S*total_Ia/(total_not_D+α*total_shieldR) - β_s*S*total_Is/(total_not_D+α*total_shieldR)
        du[idxs[2]] = dE = +β_a*S*total_Ia/(total_not_D+α*total_shieldR) + β_s*S*total_Is/(total_not_D+α*total_shieldR) - γ_e*E
        du[idxs[3]] = dIa = p[i]*γ_e*E - γ_a*Ia
        du[idxs[4]] = dIs = (1-p[i])*γ_e*E - γ_s*Is
        du[idxs[5]] = dIhsub = h[i]*(1-ξ[i])*γ_s*Is - γ_h*Ihsub
        du[idxs[6]] = dIhcrit = h[i]*ξ[i]*γ_s*Is - γ_h*Ihcrit
        du[idxs[7]] = dR = γ_a*Ia + (1-h[i])*γ_s*Is + γ_h*Ihsub + (1-μ[i])*γ_h*Ihcrit
        du[idxs[8]] = dD = μ[i]*γ_h*Ihcrit
    end
end;

function run_simulations(u0, tspan, p, stop_at)
    
    simulations = Dict()
    pop, γ_a, γ_s, αs, p = p[1], p[6], p[7], p[12], deepcopy(p)
    case_threshold = stop_at / pop.size
    
    # Set up callback to stop integration when number of cases reaches a threshold.
    # This must be equal to zero at stopping condition.
    condition(u, t, integrator) = sum(u[1:n_classes(pop)]) - (1 - case_threshold)
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!)
    
    # Solver kwargs.
    algoritm = Tsit5()
    kwargs = Dict(:abstol=>1e-8, :reltol=>1e-8, :saveat=>1, :isoutofdomain=>(y, p, t)->any(x-> x < 0, y))
    
    for scenario in ["low", "high"]
        p[6:7] .= scenario == "low" ? (γ_a[1], γ_s[1]) : (γ_a[2], γ_s[2])
        p[12] = 0
        
        # Initial spin-up.
        prob_init = ODEProblem(SerologicalShields, u0, tspan, p)
        
        # Solve the initial spin-up, checking for positivity and using callback definition.
        sol_init = solve(prob_init, algoritm, callback=cb; kwargs...)

        # Assign end spin-up states as new initial condition.
        u_init = sol_init.u[end]
        
        # Run one simulation per each α.
        solutions = []
        for α=αs
            pα = deepcopy(p)
            pα[12] = α
            prob = ODEProblem(SerologicalShields, u_init, tspan, pα)
            sol = solve(prob, algoritm; kwargs...)
            push!(solutions, (pα, sol))
        end
        simulations[scenario] = solutions
    end
    
    return simulations
end;