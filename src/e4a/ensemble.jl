include("Earth4All.jl")
##
using ModelingToolkit
using DifferentialEquations
using Plots
##
tltl = Earth4All.tltl()
##
sys = structural_simplify(tltl)
prob = ODEProblem(sys, [], (1980, 2100))

function modify_inits(x)
    return x .* (1 .+ 0.1 * randn(length(x)))
end

function modify_pars(x)
    println(x[1], " ", x[2])
    y = copy(x)
    y[1] = 0.0
    y[2] = 0.0
    return y
end

##
# This is the function that will be called to generate a new initial condition. 
# In this case, we are adding a small amount of normally distributed noise to the initial condition.
function prob_func(prob, i, repeat)
    remake(prob, u0=modify_inits(prob.u0))
end

function par_func(prob, i, repeat)
    remake(prob, p=modify_pars(prob.p))
end

##
ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
ensemble_prob = EnsembleProblem(prob, prob_func=par_func)
##
sol = solve(ensemble_prob, Euler(), EnsembleThreads(), trajectories=10; dt=0.01, dtmax=0.01)
##
summ = EnsembleSummary(sol)
##
# you can get the trajectory index with 
# ```
# println.(enumerate(states(sys))," -> ",getdescription.(states(sys)))
# ```
i = 26
fig = plot(
    summ,
    fillalpha=0.5,
    trajectories=i,
    title=(getdescription.(states(sys)))[i],
    titlefontsize=11
)
savefig(fig, "ensemble_example.png")
##
