include("LotkaVolterra.jl")
##
using ModelingToolkit
using DifferentialEquations
using Plots
##
@named model = LotkaVolterra.lotka_volterra()
##
sys = structural_simplify(model)
prob = ODEProblem(sys, [], (0, 10))

# function modify_inits(x)
#     return x .* (1 .+ 0.1 * randn(length(x)))
# end

function modify_pars(x)
    y = copy(x)
    y[1] = rand()
    y[2] = rand()
    return y
end

##
# This is the function that will be called to generate a new initial condition. 
# In this case, we are adding a small amount of normally distributed noise to the initial condition.
# function prob_func(prob, i, repeat)
#     remake(prob, u0=modify_inits(prob.u0))
# end

function modify_prob(prob, i, repeat)
    remake(prob, p=modify_pars(prob.p))
end

##
ensemble_prob = EnsembleProblem(prob, prob_func=modify_prob)
# ensemble_prob = EnsembleProblem(prob, prob_func=par_func)
##
sol = solve(ensemble_prob, Euler(), EnsembleThreads(), trajectories=10; dt=0.01, dtmax=0.01)
##
summ = EnsembleSummary(sol)
##
# you can get the trajectory index with 
# ```
# println.(enumerate(states(sys))," -> ",getdescription.(states(sys)))
# ```
i = 2
fig = plot(
    summ,
    fillalpha=0.5,
    trajectories=i,
    title=(getdescription.(states(sys)))[i],
    titlefontsize=11
)
# savefig(fig, "ensemble_example.png")
##
