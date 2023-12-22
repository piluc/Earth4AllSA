include("functions.jl")

@variables t
D = Differential(t)

function lotka_volterra(; name, params=_params, inits=_inits, tables=_tables, ranges=_ranges)
    eqs = []

    # @parameters P0 = params[:P0] [description = "Parameter 0"]
    @parameters P1 = params[:P1] [description = "Parameter 1"]
    @parameters P2 = params[:P2] [description = "Parameter 2"]
    @parameters P3 = params[:P3] [description = "Parameter 3"]
    @parameters P4 = params[:P4] [description = "Parameter 4"]

    @variables U1(t) = inits[:U1] [description = "Prey"]
    @variables U2(t) = inits[:U2] [description = "Predator"]
    # @variables U3(t) [description = "Sum of predator and prey"]
    # @variables U4(t) = 0.0 [description = "Equation for fake parameter"]

    # add_equation!(eqs, D(U4) ~ P0 - log10(10^P0))
    add_equation!(eqs, D(U1) ~ P1 * U1 - P2 * U1 * U2)
    add_equation!(eqs, D(U2) ~ -P3 * U2 + P4 * U1 * U2)
    # add_equation!(eqs, U3 ~ U1 + U2)

    return ODESystem(eqs; name=name)
end
