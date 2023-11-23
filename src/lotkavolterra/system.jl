include("functions.jl")

@variables t
D = Differential(t)

function lotka_volterra(; name, params=_params, inits=_inits, tables=_tables, ranges=_ranges)
    eqs = []

    @parameters P1 = params[:P1] [description = "Parameter 1"]
    @parameters P2 = params[:P2] [description = "Parameter 2"]
    @parameters P3 = params[:P3] [description = "Parameter 3"]
    @parameters P4 = params[:P4] [description = "Parameter 4"]

    @variables U1(t) = inits[:U1] [description = "U1"]
    @variables U2(t) = inits[:U2] [description = "U2"]

    add_equation!(eqs, D(U1) ~ P1 * U1 - P2 * U1 * U2)
    add_equation!(eqs, D(U2) ~ -P3 * U2 + P4 * U1 * U2)

    return ODESystem(eqs; name=name)
end
