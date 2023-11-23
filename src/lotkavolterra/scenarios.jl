function lk_sys(; kwargs...)
    @named lk = lotka_volterra(; kwargs...)
    systems = [lk]
    con_eqs::Vector{Equation} = []
    return WorldDynamics.compose(systems, con_eqs)
end
