using ModelingToolkit
using DifferentialEquations

function tltl_solution()
    return WorldDynamics.solve(tltl(), (1980, 2100), solver=Euler(), dt=0.015625, dtmax=0.015625)
end
