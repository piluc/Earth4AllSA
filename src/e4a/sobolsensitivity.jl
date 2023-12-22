using GlobalSensitivity, Statistics, ModelingToolkit, OrdinaryDiffEq, QuasiMonteCarlo, Plots

include("Earth4All.jl")

function e4a_ode_problem()
    @named e4a = Earth4All.earth4all()
    e4a_sys = structural_simplify(e4a)
    e4a_prob = ODEProblem(e4a_sys, [], (1980.0, 2100.0))
    return e4a_prob
end

function modify(p_from_gsa)
    e4a_prob = e4a_ode_problem()
    @named e4a = Earth4All.earth4all()
    t = collect(range(1980.0, stop=2100.0, length=2400))
    e4a_prob_modified = remake(e4a_prob; p=p_from_gsa)
    e4a_sol = solve(e4a_prob_modified, Tsit5(); saveat=t)
    return [mean(e4a_sol[e4a.AWBI])]
end

function upper_lower_bounds(prob)
    y = []

    push!(y, [prob.p[1] - 0.5, prob.p[1] + 0.5])
    for i in 2:lastindex(prob.p)
        push!(y, [prob.p[i], prob.p[i]])
    end

    y[48] = [0, 0.2]
    y[211] = [0, 0.02]
    y[212] = [0, 0.02]
    y[258] = [0, 0.01]

    y[50] = [0, 0.02]
    y[68] = [0.5, 1.6]
    y[54] = [0, 0.4]
    y[60] = [0, 0.04]

    y[237] = [0, 0.4]
    y[49] = [0, 0.04]
    y[51] = [0, 0.04]

    y[154] = [0.02, 0.02]
    y[141] = [0.05, 0.4]
    y[152] = [0.1, 1.0]
    y[150] = [0.1, 1.0]

    y[81] = [0.002, 0.008]
    y[83] = [0.5, 2.0]
    y[89] = [0.5, 2.0]
    y[78] = [0.2, 1.8]
    y[21] = [0.0, 16.0]

    y[10] = [0.0, 0.02]
    y[4] = [0.0, 0.02]
    y[129] = [5.0, 5.0]
    y[52] = [5.0, 5.0]

    return y
end

function execute_sobol(ns)
    prob = e4a_ode_problem()
    ub_lb = upper_lower_bounds(prob)
    return gsa(modify, Sobol(), ub_lb, samples=ns)
end

# execute_sobol(10);