include("Earth4All.jl")

using ModelingToolkit
using DifferentialEquations
using PlotlyJS
using Serialization
using Statistics
using WorldDynamics

global p_name = ["ERDN2OKF2022", "ERDCH4KC2022", "DACCO22100", "FGDC2022", "EETF2022", "EGTRF2022", "EPTF2022", "ETGBW", "GEIC", "FETPO", "GFCO2SCCS", "EROCEPA2022", "GFNE", "GREF", "GCWR", "GFNRM", "GFRA", "EROCFSP", "USPIS2022", "USPUS2022", "GEFR", "MIROTA2022"]
global tltl_alpha = fill(0.0, length(p_name))
global gl_alpha = fill(1.0, length(p_name))
global p_desc = ["Extra rate of decline in N2O per kg fertilizer from 2022", "Extra rate of decline in CH4 per kg crop after 2022 1/y", "Direct Air Capture of CO2 in 2100 GtCO2/y", "Fraction of Govmnt Debt Cancelled in 2022 1/y", "Extra Empowerment Tax From 2022 (share of NI)", "Extra General Tax Rate From 2022", "Extra Pension Tax From 2022 (share of NI)", "Extra Transfer of Govmnt Budget to Workers", "Goal for Extra Income from Commons (share of NI)", "Fraction of Extra Taxes Paid by Owners", "Goal for fraction of CO2-sources with CCS", "Extra ROC in energy productivity after 2022 1/y", "Goal for fraction new electrification", "Goal for renewable el fraction", "Goal for Crop Waste Reduction", "Goal for Fraction New Red Meat", "Goal for fraction regenerative agriculture", "Extra ROC in Food Sector Productivity from 2022 1/y", "Unconventional Stimulus in PIS from 2022 (share of GDP)", "Unconventional Stimulus in PUS from 2022 (share of GDP)", "Goal for Extra Fertility Reduction", "Max Imported ROTA from 2022 1/y"]

function _variables(e4a)
    variables = [
        (e4a.POP, 0, 10000, "Population"),
        (e4a.AWBI, 0, 2.4, "Average wellbeing"),
        (e4a.GDPP, 0, 60, "GDP per person"),
        (e4a.STE, 0, 2, "Social tension"),
        (e4a.INEQ, 0, 1.6, "Inequality"),
        (e4a.OW, 0, 4, "Global warming"),]
    return variables
end

function default_gl_pars(x)
    y = copy(x)
    # First group
    y[4] = 0.1
    y[19] = 0.01
    y[20] = 0.01
    y[22] = 0.005
    # Second group
    y[6] = 0.01
    y[10] = 0.8
    y[8] = 0.2
    y[9] = 0.02
    # Third group
    y[21] = 0.2
    y[5] = 0.02
    y[7] = 0.02
    # Fourth group
    y[18] = 0.02
    y[15] = 0.2
    y[17] = 0.5
    y[16] = 0.5
    # Fifth group
    y[12] = 0.004
    y[13] = 1.0
    y[14] = 1.0
    y[11] = 0.9
    y[3] = 8.0
    # Sixth group
    y[2] = 0.01
    y[1] = 0.01
    return y
end

function default_tltl_pars(x)
    y = copy(x)
    # First group
    y[4] = 0.0
    y[19] = 0.0
    y[20] = 0.0
    y[22] = 0.0
    # Second group
    y[6] = 0.0
    y[10] = 0.5
    y[8] = 0.0
    y[9] = 0.0
    # Third group
    y[21] = 0.0
    y[5] = 0.0
    y[7] = 0.00
    # Fourth group
    y[18] = 0.02
    y[15] = 0.05
    y[17] = 0.1
    y[16] = 0.1
    # Fifth group
    y[12] = 0.002
    y[13] = 0.5
    y[14] = 0.5
    y[11] = 0.2
    y[3] = 0.0
    # Sixth group
    y[2] = 0.0
    y[1] = 0.0
    return y
end

function e4a_sys_prob()
    @named e4a = Earth4All.earth4all()
    e4a_sys = structural_simplify(e4a)
    prob = ODEProblem(e4a_sys, [], (1980, 2100))
    return e4a, e4a_sys, prob
end

function compute_gl_sol(prob)
    prob = remake(prob, p=default_gl_pars(prob.p))
    sol = DifferentialEquations.solve(prob, Euler(); dt=0.015625, dtmax=0.015625)
    return sol
end

function compute_tltl_sol(prob)
    prob = remake(prob, p=default_tltl_pars(prob.p))
    sol = DifferentialEquations.solve(prob, Euler(); dt=0.015625, dtmax=0.015625)
    return sol
end

function compute_alpha_sol(prob, fixed_p, α, p)
    tltl_p = default_tltl_pars(prob.p)
    gl_p = default_gl_pars(prob.p)
    δ = gl_p - tltl_p
    for i in 1:lastindex(p)
        tltl_p[p[i]] = tltl_p[p[i]] + δ[p[i]] * α[i]
    end
    for i in 1:lastindex(fixed_p)
        tltl_p[fixed_p[i][1]] = fixed_p[i][2]
    end
    prob = remake(prob, p=tltl_p)
    sol = DifferentialEquations.solve(prob, Euler(); dt=0.015625, dtmax=0.015625)
    return sol
end

function compute_all_sols(fixed_p, α_values, p, suffix)
    @named e4a = Earth4All.earth4all()
    e4a_sys = structural_simplify(e4a)
    prob = ODEProblem(e4a_sys, [], (1980, 2100))
    sol_array = compute_all_sols(fixed_p, α_values, p, prob, e4a)
    mkpath("sols/")
    serialize("sols/sols_" * string(length(α_values)) * "_" * string(length(p)) * suffix * ".dat", sol_array)
end

function compute_all_sols(fixed_p, α_values, p, prob, e4a)
    α = fill(0, length(p))
    gl_p = default_gl_pars(prob.p)
    tltl_p = default_tltl_pars(prob.p)
    δ = gl_p - tltl_p
    sol_array = []
    compute_all_sols(length(α_values), α, 0, fixed_p, α_values, δ, prob, e4a, p, sol_array)
    return sol_array
end

function compute_all_sols(n, α, i, fixed_p, α_values, δ, prob, e4a, p, sol_array)
    if (i == length(α))
        tltl_p = default_tltl_pars(prob.p)
        for pi in 1:lastindex(p)
            tltl_p[p[pi]] = tltl_p[p[pi]] + δ[p[pi]] * α_values[α[pi]]
        end
        for i in 1:lastindex(fixed_p)
            tltl_p[fixed_p[i][1]] = fixed_p[i][2]
        end
        prob = remake(prob, p=tltl_p)
        sol = DifferentialEquations.solve(prob, Euler(); dt=0.015625, dtmax=0.015625)
        li = lastindex(sol.t)
        sol_v = [sol[e4a.AWBI][li], sol[e4a.GDPP][li], sol[e4a.INEQ][li], sol[e4a.OW][li], sol[e4a.POP][li], sol[e4a.STE][li]]
        push!(sol_array, sol_v)
        println(α)
    else
        for a in 1:n
            α[i+1] = a
            compute_all_sols(n, α, i + 1, fixed_p, α_values, δ, prob, e4a, p, sol_array)
        end
    end
end

function compute_all_subset_sols(prob, e4a)
    sol_array = []
    compute_all_subset_sols(fill(1, length(p_name)), length(p_name) + 1, prob, e4a, sol_array)
    mkpath("sols/")
    serialize("sols/sols_all_subsets.dat", sol_array)
end

function compute_all_subset_sols(a, b, prob, e4a, sol_array)
    if (b == 1)
        tltl_p = default_tltl_pars(prob.p)
        gl_p = default_gl_pars(prob.p)
        for pi in 1:lastindex(tltl_p)
            if (a[pi] == 1)
                tltl_p[pi] = gl_p[pi]
            end
        end
        prob = remake(prob, p=tltl_p)
        sol = DifferentialEquations.solve(prob, Euler(); dt=0.015625, dtmax=0.015625)
        li = lastindex(sol.t)
        sol_v = [sol[e4a.AWBI][li], sol[e4a.GDPP][li], sol[e4a.INEQ][li], sol[e4a.OW][li], sol[e4a.POP][li], sol[e4a.STE][li]]
        push!(sol_array, sol_v)
        println(a)
    else
        a[b-1] = 0
        compute_all_subset_sols(a, b - 1, prob, e4a, sol_array)
        a[b-1] = 1
        compute_all_subset_sols(a, b - 1, prob, e4a, sol_array)
    end
end

function check_all_subsets(sol_fn, dom_fn, apx_fn, e4a, gl_sol, abs_tol, rel_tol)
    sol_array = deserialize(sol_fn)
    li = lastindex(gl_sol.t)
    df = open(dom_fn, "w")
    af = open(apx_fn, "w")
    for r in lastindex(sol_array):-1:1
        if (sol_array[r][1] >= gl_sol[e4a.AWBI][li] && sol_array[r][2] >= gl_sol[e4a.GDPP][li] && sol_array[r][3] <= gl_sol[e4a.INEQ][li] && sol_array[r][4] <= gl_sol[e4a.OW][li] && sol_array[r][6] <= gl_sol[e4a.STE][li]) # && sol_array[r][5] >= gl_sol[e4a.POP][li])
            println("          ", r, " is a dominator of GL")
            write(df, string(r) * "\n")
        end
        if (isapprox(sol_array[r][1], gl_sol[e4a.AWBI][li]; atol=abs_tol, rtol=rel_tol) && isapprox(sol_array[r][2], gl_sol[e4a.GDPP][li]; atol=abs_tol, rtol=rel_tol) && isapprox(sol_array[r][3], gl_sol[e4a.INEQ][li]; atol=abs_tol, rtol=rel_tol) && isapprox(sol_array[r][4], gl_sol[e4a.OW][li]; atol=abs_tol, rtol=rel_tol) && isapprox(sol_array[r][6], gl_sol[e4a.STE][li]; atol=abs_tol, rtol=rel_tol)) # && isapprox(sol_array[r][5]; gl_sol[e4a.POP][li]; atol=abs_tol, rtol=rel_tol))
            println("          ", r, " is approximately the same as GL")
            write(af, string(r) * "\n")
        end
        if (r % 1000 == 0)
            println("Processed ", r, " sunsets")
        end
    end
    close(df)
    close(af)
end

function find_dominators(sol_fn, dom_fn, α_values, np)
    println("Computing ODE system and giant leap solution")
    e4a, _, prob = e4a_sys_prob()
    gl_sol = compute_gl_sol(prob)
    println("Reading solutions")
    sol_array = deserialize(sol_fn)
    println("Looking for dominators")
    f = open(dom_fn, "w")
    find_dominators(e4a, gl_sol, sol_array, f, α_values, np)
    close(f)
    println("Finished")
end

function find_dominators(e4a, gl_sol, sol_array, f, α_values, np)
    α = fill(0, np)
    find_dominators(length(α_values), α, 0, 0, e4a, gl_sol, sol_array, f, α_values)
end

function find_dominators(n, α, i, r, e4a, gl_sol, sol_array, f, α_values)
    if (i == length(α))
        r = r + 1
        if (r % 1000 == 0)
            println(r, " solutions have been analyzed")
        end
        li = lastindex(gl_sol.t)
        if (sol_array[r][1] >= gl_sol[e4a.AWBI][li] && sol_array[r][2] >= gl_sol[e4a.GDPP][li] && sol_array[r][3] <= gl_sol[e4a.INEQ][li] && sol_array[r][4] <= gl_sol[e4a.OW][li] && sol_array[r][6] <= gl_sol[e4a.STE][li]) # && sol_array[r][5] >= gl_sol[e4a.POP][li])
            p = fill(0.0, length(α))
            for j in 1:lastindex(p)
                p[j] = α_values[α[j]]
            end
            println(p)
            write(f, string(p) * "\n")
        end
        return r
    end
    for a in 1:n
        α[i+1] = a
        r = find_dominators(n, α, i + 1, r, e4a, gl_sol, sol_array, f, α_values)
    end
    return r
end

function find_dominator_indices(sol_fn)
    println("Computing ODE system and giant leap solution")
    e4a, _, prob = e4a_sys_prob()
    gl_sol = compute_gl_sol(prob)
    println("Reading solutions")
    sol_array = deserialize(sol_fn)
    println("Looking for dominators")
    di = []
    li = lastindex(gl_sol.t)
    for s in 1:lastindex(sol_array)
        sol = sol_array[s]
        if (sol[1] >= gl_sol[e4a.AWBI][li] && sol[2] >= gl_sol[e4a.GDPP][li] && sol[3] <= gl_sol[e4a.INEQ][li] && sol[4] <= gl_sol[e4a.OW][li] && sol[6] <= gl_sol[e4a.STE][li]) # && sol_array[r][5] >= gl_sol[e4a.POP][li])
            println(s)
            push!(di, s)
        end
    end
    println("Finished")
    return di
end

function plot_two_sols(scen1, sol1, scen2, sol2, vars)
    x = range(1, 7681, length=7681)
    traces = GenericTrace[]
    for v in 1:lastindex(vars)
        desc = getdescription(vars[v])
        trace1 = PlotlyJS.scatter(x=x, y=sol1[vars[v]], name=desc * "-" * scen1, line=attr(color="royalblue", dash="dash"))
        trace2 = PlotlyJS.scatter(x=x, y=sol2[vars[v]], name=desc * "-" * scen2, line=attr(color="firebrick", dash="dot"))
        push!(traces, trace1)
        push!(traces, trace2)
    end
    return PlotlyJS.plot(traces, Layout(title="GL versus dominator"))
end

function plot_sol(e4a, _sol, _title, sa, sl; kwargs...)
    @variables t
    return WorldDynamics.plotvariables(_sol, (t, 1980, 2100), _variables(e4a); title=_title, showaxis=sa, showlegend=sl, kwargs...)
end

function find_α(index, nα, np)
    α_values = [0.0, 0.5, 1, 1.5, 2.0]
    α = fill(0, np)
    α, _ = find_α(nα, α, 0, index, 0)
    p = fill(0.0, np)
    for i in 1:np
        p[i] = α_values[α[i]]
    end
    return p
end

function find_α(n, α, i, index, r)
    if (i == length(α))
        r = r + 1
        return α, r
    end
    β = fill(-1, length(α))
    for a in 1:n
        α[i+1] = a
        β, r = find_α(n, α, i + 1, index, r)
        if (r == index)
            break
        end
    end
    return β, r
end

function find_index(β, nα, np)
    α = fill(0, np)
    return find_index(nα, α, 0, β, 0)
end

function find_index(n, α, i, β, index)
    if (i == length(α))
        index = index + 1
        if (α == β)
            return index, true
        else
            return index, false
        end
    end
    f = false
    for a in 1:n
        α[i+1] = a
        index, f = find_index(n, α, i + 1, β, index)
        if (f)
            break
        end
    end
    return index, f
end

function read_alphas(α_fn)
    f = open(α_fn, "r")
    lines = readlines(f)
    close(f)
    return lines
end

function find_dominators_from_files(sol_fn, α_fn, doms_fn)
    println("Computing ODE system and giant leap solution")
    e4a, _, prob = e4a_sys_prob()
    gl_sol = compute_gl_sol(prob)
    println("Reading solutions")
    sol_array = deserialize(sol_fn)
    find_dominators_from_files(e4a, gl_sol, sol_array, α_fn, doms_fn)
end

function find_dominators_from_files(e4a, gl_sol, sol_array, α_fn, doms_fn)
    println("Reading α combinations")
    α_array = read_alphas(α_fn)
    println("Looking for dominators")
    li = lastindex(gl_sol.t)
    vv = [gl_sol[e4a.AWBI][li], gl_sol[e4a.GDPP][li], gl_sol[e4a.INEQ][li], gl_sol[e4a.OW][li], gl_sol[e4a.POP][li], gl_sol[e4a.STE][li]]
    f = open(doms_fn, "w")
    for r in 1:lastindex(sol_array)
        if (sol_array[r][1] >= vv[1] && sol_array[r][2] >= vv[2] && sol_array[r][3] <= vv[3] && sol_array[r][4] <= vv[4] && sol_array[r][6] <= vv[6]) # && sol_array[r][2][5] >= vv[5])
            println(α_array[r])
            write(f, α_array[r] * "\n")
        end
    end
    close(f)
    println("Finished")
end
