using GlobalSensitivity, Statistics, ModelingToolkit, OrdinaryDiffEq, QuasiMonteCarlo, Plots

include("Earth4All.jl")

function e4a_ode_problem()
    @named e4a = Earth4All.earth4all()
    e4a_sys = structural_simplify(e4a)
    e4a_prob = ODEProblem(e4a_sys, [], (1980.0, 2100.0))
    return e4a_prob
end

function modify(p)
    e4a_prob = e4a_ode_problem()
    @named e4a = Earth4All.earth4all()
    t = collect(range(1980.0, stop=2100.0, length=2400))
    e4a_prob_modified = remake(e4a_prob; p=p)
    e4a_sol = solve(e4a_prob_modified, Tsit5(); saveat=t)
    return [mean(e4a_sol[e4a.AWBI])]
end

function upper_lower_bounds(prob)
    y = []
    push!(y, [prob.p[1] - 0.5, prob.p[1] + 0.5])
    for i in 2:lastindex(prob.p)
        push!(y, [prob.p[i], prob.p[i]])
    end
    y[47] = [0, 0.2]
    y[210] = [0, 0.02]
    y[211] = [0, 0.02]
    y[257] = [0, 0.01]
    return y
end

function execute_gsa(tnt, nt, pl)
    e4a_prob = e4a_ode_problem()
    ub_lb = upper_lower_bounds(e4a_prob)
    return gsa(modify, Morris(total_num_trajectory=tnt, num_trajectory=nt, len_design_mat=pl), ub_lb)
end
