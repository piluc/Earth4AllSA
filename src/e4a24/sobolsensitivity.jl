using GlobalSensitivity, Statistics, ModelingToolkit, OrdinaryDiffEq, QuasiMonteCarlo, PlotlyJS, Serialization

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
    return [mean(e4a_sol[e4a.AWBI][end-120:end]), mean(e4a_sol[e4a.AWBI][end-240:end]), mean(e4a_sol[e4a.GDPP][end-120:end]), mean(e4a_sol[e4a.GDPP][end-240:end]), mean(e4a_sol[e4a.INEQ][end-120:end]), mean(e4a_sol[e4a.INEQ][end-240:end]), mean(e4a_sol[e4a.OW][end-120:end]), mean(e4a_sol[e4a.OW][end-240:end]), mean(e4a_sol[e4a.POP][end-120:end]), mean(e4a_sol[e4a.POP][end-240:end]), mean(e4a_sol[e4a.STE][end-120:end]), mean(e4a_sol[e4a.STE][end-240:end])]
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

    y[18] = [0.02, 0.02]
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

function execute_sobol(ns)
    prob = e4a_ode_problem()
    ub_lb = upper_lower_bounds(prob)
    return gsa(modify, Sobol(), ub_lb, samples=ns)
end

function bar_plots(res)
    @named e4a = Earth4All.earth4all()
    fp = "/Users/piluc/Downloads/sobol/plots/"
    mkpath(fp)
    v = ["Average wellbeing index", "GDP per person", "Inequality", "Observed warming", "Population", "Social tension"]
    a = ["awbi", "gdpp", "ineq", "ow", "pop", "ste"]
    for i in 1:lastindex(v)
        p = plot([bar(x=getdescription.(parameters(e4a)), y=res.S1[2*i-1, :], name="Sobol first index"), bar(x=getdescription.(parameters(e4a)), y=res.ST[2*i-1, :], name="Sobol total index")], Layout(title=v[i] * " (averaged over 6 years)"))
        savefig(p, fp * a[i] * "_06.html")
        p = plot([bar(x=getdescription.(parameters(e4a)), y=res.S1[2*i, :], name="Sobol first index"), bar(x=getdescription.(parameters(e4a)), y=res.ST[2*i, :], name="Sobol total index")], Layout(title=v[i] * " (averaged over 12 years)"))
        savefig(p, fp * a[i] * "_12.html")
    end
end

# start_time = time()
# res = execute_sobol(6)
# println("Sobol executed in ", time() - start_time, " seconds")
# serialize("sobol.dat", res)
