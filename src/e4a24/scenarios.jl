function tltl(; kwargs...)
    @named e4a = earth4all(; kwargs...)
    systems = [e4a]
    con_eqs::Vector{Equation} = []
    return WorldDynamics.compose(systems, con_eqs)
end
