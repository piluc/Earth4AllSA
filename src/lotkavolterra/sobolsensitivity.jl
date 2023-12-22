using GlobalSensitivity, Statistics, ModelingToolkit, OrdinaryDiffEq, QuasiMonteCarlo, Plots

include("LotkaVolterra.jl")

global tmp = 0

function lk_ode_problem()
    @named lk = LotkaVolterra.lotka_volterra()
    sys = structural_simplify(lk)
    prob = ODEProblem(sys, [], (0.0, 10.0))
    return prob
end

function modify(p_from_gsa)
    global tmp = tmp + 1
    prob = lk_ode_problem()
    @named lk = LotkaVolterra.lotka_volterra()
    t = collect(range(0, stop=10, length=200))
    prob1 = remake(prob; p=p_from_gsa)
    sol = solve(prob1, Tsit5(); saveat=t)
    return [sol[lk.U1][100], sol[lk.U2][100]]
end

function upper_lower_bounds(prob)
    y = []
    push!(y, [prob.p[1] - 0.5, prob.p[1] + 0.5])
    for i in 2:lastindex(prob.p)
        push!(y, [prob.p[i] - 0.5, prob.p[i] + 0.5])
        # push!(y, [prob.p[i], prob.p[i]])
    end
    y = [[0.02, 0.044], [0.44, 0.68], [0.0226, 0.0354], [0.71, 1.15]]
    return y
end

function execute_sobol(ns)
    prob = lk_ode_problem()
    ub_lb = upper_lower_bounds(prob)
    return gsa(modify, Sobol(), ub_lb, samples=ns)
end

# execute_sobol(10);
