include("Earth4All.jl")

using ModelingToolkit
using DifferentialEquations
using Plots
using Statistics

global p1 = 3
global p2 = 9
global p3 = 11
global p4 = 13
global α_values = [0.0, 0.5, 1, 1.5, 2.0]
global α1 = 1
global α2 = 1
global α3 = 1
global α4 = 0

function modify_pars(x)
    global α1, α2, α3, α4
    # println(α1, α2, α3, α4)
    if (α4 == 5)
        α4 = 1
        if (α3 == 5)
            α3 = 1
            if (α2 == 5)
                α2 = 1
                α1 = α1 + 1
            else
                α2 = α2 + 1
            end
        else
            α3 = α3 + 1
        end
    else
        α4 = α4 + 1
    end
    y = copy(x)
    y[p1] += 8.0 * α_values[α1]
    y[p2] += 0.02 * α_values[α2]
    y[p3] += 0.7 * α_values[α3]
    y[p4] += 0.5 * α_values[α4]
    # # First group
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
    return sol
end

# ensemble(10)