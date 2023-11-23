_params = Dict{Symbol,Float64}(
    :P1 => 1.5,
    :P2 => 1.0,
    :P3 => 3.0,
    :P4 => 1.0,
)

getparameters() = copy(_params)
