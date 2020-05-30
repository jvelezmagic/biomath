mutable struct PoorSoul <: AbstractAgent
    id::Int
    pos::NTuple{2, Float64}
    vel::NTuple{2, Float64}
    mass::Float64
    days_infected::Int # Number of days since is infected.
    status::Symbol # :S, :I , :R or :V
    β::Float64 # Transmission probability. It reflects the level of hygiene of an individual.
end