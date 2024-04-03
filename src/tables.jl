_tables = Dict{Symbol,Tuple{Vararg{Float64}}}(
    :ROCWSO => (
        0.06,
        0.02,
        0,
        -0.007,
        -0.01,
    ),
    :IEST => (
        1,
        1,
        0,
    ),
    :PSESTR => (
        0,
        1,
    ),
)
_ranges = Dict{Symbol,Tuple{Float64,Float64}}(
    :IEST => (0, 2),
    :PSESTR => (0, 1),
    :ROCWSO => (0, 2),
)

gettables() = copy(_tables)
getranges() = copy(_ranges)
