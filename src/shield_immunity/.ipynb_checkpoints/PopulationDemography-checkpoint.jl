@with_kw struct PopulationDemography
    class_structure # Fraction of population in each class.
    size # Number of people in populations.
    n_states # Epidemic States.
end;

n_classes(pop::PopulationDemography) = length(pop.class_structure)
n_variables(pop::PopulationDemography) = n_classes(pop) * pop.n_states
idxs(pop::PopulationDemography) = collect(1:n_variables(pop))

function idx_matrix(pop::PopulationDemography)
    pop_idxs = idxs(pop)
    return reshape(pop_idxs, n_classes(pop), pop.n_states)
end;

function idxs_per_state(pop::PopulationDemography)
    m = idx_matrix(pop)
    return [m[:, i] for i=1:size(m, 2)]
end;

function idxs_per_class(pop::PopulationDemography)
    m = idx_matrix(pop)
    return [m[i, :] for i=1:size(m, 1)]
end;