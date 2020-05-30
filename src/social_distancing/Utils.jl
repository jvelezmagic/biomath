function agent_color(a)
    
    color_scheme = Dict(
        :S => :black,#"#2b2b33",
        :I => :red, #"#e58c8a",
        :R => :green,#"#80ff72",
        :V => :blue#"#7ee8fa"
    )
    return color_scheme[a.status]
end;

function gif_SIR(model, n_days, check_each; out="tmp.gif", fps=25)
    model = deepcopy(model)
    e = model.space.extend
    animation = @animate for i in 0:check_each:((n_days) * steps_per_day)
        p1 = plotabm(
            model;
            ac = agent_color,
            as = 4,
            msc=:auto,
            showaxis = false,
            grid = false,
            xlims = (0, e[1]),
            ylims = (0, e[2]),
        )

        day = div(i, steps_per_day)
        
        title = string("Day: ", day, " Step: ", i, "\nLockdown: ", model.active_lockdown, " Vaccination: ", model.active_vaccination)
        title!(p1, title)
        step!(model, sir_agent_step!, sir_model_step!, check_each)
    end
    
    return gif(animation, out, fps = fps)
end

suceptible(x) = count(i == :S for i in x)
infected(x) = count(i == :I for i in x)
recovered(x) = count(i == :R for i in x)
vaccinated(x) = count(i == :V for i in x)

function simulate(dict; out_dir = nothing, return_results = false)
        
    if out_dir == nothing && !return_results
        return nothing
    end
    
    if out_dir ≠ nothing
        model_label = savename(dict, "csv")
        out_file = joinpath(out_dir, model_label)
        
        if isfile(out_file)        
            return return_results ? CSV.read(out_file) : nothing 
        end
    end
    
    adata = [
        (:status, suceptible),
        (:status, infected),
        (:status, recovered),
        (:status, vaccinated)
    ]
        
    dict = deepcopy(dict)
    n_steps = pop!(dict, :nsteps)
    model = sir_initiation(;dict...)
    
    results, _ = run!(model, sir_agent_step!, sir_model_step!, n_steps; adata = adata)
    
    if return_results && out_dir ≠ nothing
        CSV.write(out_file, results)
        return results
    elseif return_results
        return results
    else
        CSV.write(out_file, results)
        return nothing
    end
end