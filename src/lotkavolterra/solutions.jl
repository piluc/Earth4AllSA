using ModelingToolkit
using DifferentialEquations

function lk_solution()
    return WorldDynamics.solve(lk_sys(), (0.0, 10.0), solver=Euler(), dt=0.015625, dtmax=0.015625)
end
