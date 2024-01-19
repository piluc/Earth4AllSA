include("Earth4All.jl")
##
using ModelingToolkit
using DifferentialEquations
using Plots
##
# @named e4a = Earth4All.earth4all()
##
# e4a_sys = structural_simplify(e4a)
# prob = ODEProblem(e4a_sys, [], (1980, 2100))

# function modify_inits(x)
#     return x .* (1 .+ 0.1 * randn(length(x)))
# end

function modify_pars(x)
    y = copy(x)
    y[4] = 0.1
    y[19] = 0.01
    y[20] = 0.01
    y[22] = 0.005

    y[6] = 0.01
    y[10] = 0.8
    y[8] = 0.2
    y[9] = 0.02

    y[21] = rand() * 0.4
    y[5] = 0.02
    y[7] = 0.02

    y[18] = 0.01
    y[15] = 0.2
    y[17] = 0.5
    y[16] = 0.5

    y[12] = 0.004
    y[13] = 1.0
    y[14] = 1.0
    y[11] = 0.9
    y[3] = 8.0

    y[2] = 0.01
    y[1] = 0.01
    return y
end

##
# This is the function that will be called to generate a new initial condition. 
# In this case, we are adding a small amount of normally distributed noise to the initial condition.
# function prob_func(prob, i, repeat)
#     remake(prob, u0=modify_inits(prob.u0))
# end

function par_func(prob, i, repeat)
    remake(prob, p=modify_pars(prob.p))
end

function out_func(sol, i)
    @named e4a = Earth4All.earth4all()
    return (sol[e4a.POP], false)
end

function ensemble(nt)
    @named e4a = Earth4All.earth4all()
    e4a_sys = structural_simplify(e4a)
    prob = ODEProblem(e4a_sys, [], (1980, 2100))
    ensemble_prob = EnsembleProblem(prob, output_func=out_func, prob_func=par_func)
    sol = solve(ensemble_prob, Euler(), EnsembleThreads(), trajectories=nt; dt=0.015625, dtmax=0.015625)
    avg_s = sum(sol[i] for i in 1:nt) / nt
    min_s = minimum(sol[i] for i in 1:nt)
    max_s = maximum(sol[i] for i in 1:nt)
    errlb_s = avg_s - min_s
    errub_s = max_s - avg_s
    plot(avg_s; ribbon=(errlb_s, errub_s))
end

# ensemble(10)