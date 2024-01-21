include("Earth4All.jl")

using Dates
using ModelingToolkit
using DifferentialEquations
using Plots
using Serialization
using Statistics

global p = [3, 9, 11, 13, 14, 19, 20, 21, 22]
global δ = [8.0, 0.02, 0.7, 0.5, 0.5, 0.01, 0.01, 0.2, 0.005]
global α_values = [0.0, 0.5, 1, 1.5, 2.0]
global α = [1, 1, 1, 1, 1, 1, 1, 1, 0]


function modify_pars(x)
    i = length(α)
    while (α[i] == length(α_values))
        α[i] = 1
        i = i - 1
    end
    if (i == 1)
        println("Finished α[1]=", α[i], " at time ", string(Dates.Time(Dates.now())))
    end
    α[i] = α[i] + 1
    y = copy(x)
    for i in 1:lastindex(p)
        y[p[i]] += δ[i] * α_values[α[i]]
    end
    # First group
    # y[4] = 0.1
    # y[19] = 0.01
    # y[20] = 0.01
    # y[22] = 0.005
    # # Second group
    # y[6] = 0.01
    # y[10] = 0.8
    # y[8] = 0.2
    # y[9] = 0.02
    # # Third group
    # y[21] = 0.2
    # y[5] = 0.02
    # y[7] = 0.02
    # # Fourth group
    # y[18] = 0.02
    # y[15] = 0.2
    # y[17] = 0.5
    # y[16] = 0.5
    # # Fifth group
    # y[12] = 0.004
    # y[13] = 1.0
    # y[14] = 1.0
    # y[11] = 0.9
    # y[3] = 8.0
    # # Sixth group
    # y[2] = 0.01
    # y[1] = 0.01
    return y
end

function par_func(prob, i, repeat)
    remake(prob, p=modify_pars(prob.p))
end

function out_func(sol, i)
    @named e4a = Earth4All.earth4all()
    li = lastindex(sol.t)
    return ([sol[e4a.AWBI][li], sol[e4a.GDPP][li], sol[e4a.INEQ][li], sol[e4a.OW][li], sol[e4a.POP][li], sol[e4a.STE][li]], false)
end

function ensemble(nt)
    @named e4a = Earth4All.earth4all()
    e4a_sys = structural_simplify(e4a)
    prob = ODEProblem(e4a_sys, [], (1980, 2100))
    ensemble_prob = EnsembleProblem(prob, output_func=out_func, prob_func=par_func)
    sol = solve(ensemble_prob, Euler(), EnsembleThreads(), trajectories=nt; dt=0.015625, dtmax=0.015625)
    # avg_s = sum(sol[i] for i in 1:nt) / nt
    # min_s = minimum(sol[i] for i in 1:nt)
    # max_s = maximum(sol[i] for i in 1:nt)
    # errlb_s = avg_s - min_s
    # errub_s = max_s - avg_s
    # plot(avg_s; ribbon=(errlb_s, errub_s))
    serialize("sols/sol_5_9.dat", sol)
end

# ensemble(10)
