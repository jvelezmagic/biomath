const steps_per_day = 24

function sir_initiation(;
    # Agents.
    N = 1000,
    interaction_radius = 0.012,
    speed = 0.002,
    βmax = 0.8,
    βmin = 0.4,
    
    # Pandemic.
    death_rate = 0.044,
    detection_time = 14 * steps_per_day,
    infection_period = 30 * steps_per_day,
    initial_infected = 1,
    reinfection_probability = 0.05,
    
    # Interventions.
    start_lockdown = Inf,
    start_vaccination = Inf,
    stop_lockdown_after = Inf,
    periodic_lockdown = false,
    interval_between_locks = Inf,
    fraction_isolated = 0.0, # In percentage.
    vaccination_probability = 0.001,
    vaccine_effectiveness = 1.0,
        
    # Extra.
    dt = 1.0,
    seed = 42
)
    # Assign access properties in the model.
    stop_lockdown = start_lockdown + stop_lockdown_after
    ind_isolated = nothing
    active_lockdown, active_vaccination, n_step = false, false, 0
    properties = @dict(
        N,
        active_lockdown,
        active_vaccination,
        death_rate,
        detection_time,
        dt,
        infection_period,
        interaction_radius,
        interval_between_locks,
        ind_isolated,
        fraction_isolated,
        n_step,
        periodic_lockdown,
        reinfection_probability,
        seed,
        speed,
        start_lockdown,
        start_vaccination,
        stop_lockdown,
        stop_lockdown_after,
        vaccination_probability,
        vaccine_effectiveness,
    )
    
    space = ContinuousSpace(2)
    model = ABM(PoorSoul, space, properties = properties)
    
    # Add initial individual.
    Random.seed!(seed) # Ensure reproducibility.
    for ind in 1:N
        pos = Tuple(rand(2))
        status = ind ≤ N - initial_infected ? :S : :I
        mass = 1.0
        vel = sincos(2π * rand()) .* speed
        
        # Very high transmission probability.
        # We are modelling close encounters after all.
        β = (βmax - βmin) * rand() + βmin
        add_agent!(pos, model, vel, mass, 0, status, β)
    end
    
    Agents.index!(model)
    return model
end

function transmit!(a1, a2, rp, ve)
    # For transmission, only 1 can have the disease (otherwise nothing happends). 
    count(a.status == :I for a in (a1, a2)) ≠ 1 && return
    infected, healthy = a1.status == :I ? (a1, a2) : (a2, a1)
    
    rand() > infected.β && return
    
    if healthy.status == :R
        rand() > rp && return
    elseif healthy.status == :V
        ve > rand() && return
    end
    
    healthy.status = :I
end

function sir_model_step!(model)
    r = model.interaction_radius
    activate_lockdown_or_vaccination!(model)
    for (a1, a2) in interacting_pairs(model, r, :nearest)
        transmit!(a1, a2, model.reinfection_probability, model.vaccine_effectiveness)
        elastic_collision!(a1, a2, :mass)
    end
    model.n_step += 1
end

function sir_agent_step!(agent, model)
    move_agent!(agent, model, model.dt)
    update!(agent)
    vaccinate_recover_or_die!(agent, model)
end

update!(agent) = agent.status == :I && (agent.days_infected += 1)

function vaccinate_recover_or_die!(agent, model)
    
    if agent.status in [:S, :R] && model.active_vaccination
        if rand() ≤ model.vaccination_probability
           agent.status = :V 
        end
    elseif agent.days_infected ≥ model.infection_period
        if rand() ≤ model.death_rate
            kill_agent!(agent, model)
        else
            agent.status = :R
            agent.days_infected = 0
        end
    end
end

function activate_lockdown_or_vaccination!(model)
    
    if !model.active_lockdown &&
        model.n_step ≥ model.start_lockdown &&
        model.n_step ≤ model.stop_lockdown
        
        model.active_lockdown = true
        
        n_agents = nagents(model)
        n_isolated = Int(ceil(n_agents * model.fraction_isolated))
        if n_isolated > n_agents n_isolated = n_agents end
        
        Random.seed!(model.seed)
        ind_isolated = sample(1:n_agents, n_isolated; replace=false)
        
        for (i, agent) in enumerate(allagents(model))
            if i in ind_isolated
                agent.mass = Inf
                agent.vel = (0.0, 0.0)
            end
        end
        model.ind_isolated = ind_isolated

    elseif model.active_lockdown && model.n_step == model.stop_lockdown
        
        if model.periodic_lockdown
            model.start_lockdown = model.n_step + model.interval_between_locks
            model.stop_lockdown = model.start_lockdown + model.stop_lockdown_after 
        end
        
        model.active_lockdown = false
        
        Random.seed!(model.seed)
        for(i, agent) in enumerate(allagents(model))
            if i in model.ind_isolated
                agent.mass = 1.0
                agent.vel = sincos(2π * rand()) .* model.speed
            end
        end
        model.seed += 1
    end
    
    if !model.active_vaccination && model.n_step ≥ model.start_vaccination ≥ 0
        model.active_vaccination = true
    end
end