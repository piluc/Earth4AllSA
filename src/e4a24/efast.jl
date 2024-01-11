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

    for i in 1:lastindex(prob.p)
        push!(y, [prob.p[i], prob.p[i]])
    end

    y[4] = [0, 0.2]
    y[19] = [0, 0.02]
    y[20] = [0, 0.02]
    y[22] = [0, 0.01]

    y[6] = [0, 0.02]
    y[10] = [0.5, 1.6]
    y[8] = [0, 0.4]
    y[9] = [0, 0.04]

    y[21] = [0, 0.4]
    y[5] = [0, 0.04]
    y[7] = [0, 0.04]

    y[18] = [0.02, 0.04]
    y[15] = [0.05, 0.4]
    y[17] = [0.1, 1.0]
    y[16] = [0.1, 1.0]

    y[12] = [0.002, 0.008]
    y[13] = [0.5, 2.0]
    y[14] = [0.5, 2.0]
    y[11] = [0.2, 1.8]
    y[3] = [0.0, 16.0]

    y[2] = [0.0, 0.02]
    y[1] = [0.0, 0.02]

    return y
end

function execute_efast(ns)
    prob = e4a_ode_problem()
    ub_lb = upper_lower_bounds(prob)
    return gsa(modify, eFAST(num_harmonics=8), ub_lb, samples=ns)
end

# execute_efast(10);