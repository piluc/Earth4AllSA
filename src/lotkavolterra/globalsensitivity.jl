using GlobalSensitivity, Statistics, ModelingToolkit, OrdinaryDiffEq, QuasiMonteCarlo, Plots

include("LotkaVolterra.jl")

function lk_ode_problem()
    @named lk = LotkaVolterra.lotka_volterra()
    sys = structural_simplify(lk)
    prob = ODEProblem(sys, [], (0.0, 10.0))
    return prob
end

function modify(p_from_gsa)
    prob = lk_ode_problem()
    t = collect(range(0, stop=10, length=200))
    prob1 = remake(prob; p=p_from_gsa)
    sol = solve(prob1, Tsit5(); saveat=t)
    return [mean(sol[1, :]), mean(sol[2, :])]
end

function upper_lower_bounds(prob)
    y = []
    for i in 1:lastindex(prob.p)
        push!(y, [prob.p[i] - 0.5, prob.p[i] + 0.5])
    end
    return y
end

function execute_gsa(tnt, nt)
    prob = lk_ode_problem()
    ub_lb = upper_lower_bounds(prob)
    return gsa(modify, Morris(total_num_trajectory=tnt, num_trajectory=nt), ub_lb)
end
