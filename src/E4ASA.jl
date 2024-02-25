include("Earth4All.jl")

using ColorTypes
using ColorSchemes
using ModelingToolkit
using DifferentialEquations
using PlotlyJS
using Serialization
using Statistics
using WorldDynamics

global p_name = ["ERDN2OKF2022", "ERDCH4KC2022", "DACCO22100", "FGDC2022", "EETF2022", "EGTRF2022", "EPTF2022", "ETGBW", "GEIC", "FETPO", "GFCO2SCCS", "EROCEPA2022", "GFNE", "GREF", "GCWR", "GFNRM", "GFRA", "EROCFSP", "USPIS2022", "USPUS2022", "GEFR", "MIROTA2022"]
global p_desc = ["Extra rate of decline in N2O per kg fertilizer from 2022", "Extra rate of decline in CH4 per kg crop after 2022 1/y", "Direct Air Capture of CO2 in 2100 GtCO2/y", "Fraction of Govmnt Debt Cancelled in 2022 1/y", "Extra Empowerment Tax From 2022 (share of NI)", "Extra General Tax Rate From 2022", "Extra Pension Tax From 2022 (share of NI)", "Extra Transfer of Govmnt Budget to Workers", "Goal for Extra Income from Commons (share of NI)", "Fraction of Extra Taxes Paid by Owners", "Goal for fraction of CO2-sources with CCS", "Extra ROC in energy productivity after 2022 1/y", "Goal for fraction new electrification", "Goal for renewable el fraction", "Goal for Crop Waste Reduction", "Goal for Fraction New Red Meat", "Goal for fraction regenerative agriculture", "Extra ROC in Food Sector Productivity from 2022 1/y", "Unconventional Stimulus in PIS from 2022 (share of GDP)", "Unconventional Stimulus in PUS from 2022 (share of GDP)", "Goal for Extra Fertility Reduction", "Max Imported ROTA from 2022 1/y"]
global p_desc_sorted = ["Direct Air Capture of CO2 in 2100 GtCO2/y", "Extra Empowerment Tax From 2022 (share of NI)", "Extra General Tax Rate From 2022", "Extra Pension Tax From 2022 (share of NI)", "Extra rate of decline in CH4 per kg crop after 2022 1/y", "Extra rate of decline in N2O per kg fertilizer from 2022", "Extra ROC in energy productivity after 2022 1/y", "Extra Transfer of Govmnt Budget to Workers", "Fraction of Extra Taxes Paid by Owners", "Fraction of Govmnt Debt Cancelled in 2022 1/y", "Goal for Crop Waste Reduction", "Goal for Extra Fertility Reduction", "Goal for Extra Income from Commons (share of NI)", "Goal for fraction of CO2-sources with CCS", "Goal for fraction new electrification", "Goal for Fraction New Red Meat", "Goal for fraction regenerative agriculture", "Goal for renewable el fraction", "Max Imported ROTA from 2022 1/y", "Unconventional Stimulus in PIS from 2022 (share of GDP)", "Unconventional Stimulus in PUS from 2022 (share of GDP)"]
global p_map = [3, 5, 6, 7, 2, 1, 12, 8, 10, 4, 15, 21, 9, 11, 13, 16, 17, 14, 22, 19, 20]
global tltl_alpha = fill(0.0, length(p_name))
global gl_alpha = fill(1.0, length(p_name))
global v_name = ["Average wellbeing", "GDP per person", "Inequality", "Observed warming", "Population", "Social tension"]

function _variables(e4a)
    variables = [
        (e4a.POP, 0, 10000, "Population"),
        (e4a.AWBI, 0, 5.0, "Average wellbeing"),
        (e4a.GDPP, 0, 65, "GDP per person"),
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
    y[18] = 0.0
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
    y[18] = 0.0
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

function local_sensitivity(e4a, prob)
    α_values = [0.0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.0]
    gl_p = default_gl_pars(prob.p)
    tltl_p = default_tltl_pars(prob.p)
    prob = remake(prob, p=tltl_p)
    tltl_sol = DifferentialEquations.solve(prob, Euler(); dt=0.015625, dtmax=0.015625)
    perc = []
    for p in 1:lastindex(gl_p)
        println("Computing percentages for parameter ", p)
        δ = gl_p[p] - tltl_p[p]
        sol = []
        for α in 1:lastindex(α_values)
            tltl_p = default_tltl_pars(prob.p)
            tltl_p[p] = tltl_p[p] + δ * α_values[α]
            if (tltl_p[p] <= 1.0 || !(p in [4, 10, 11, 13, 14, 16, 17, 21]))
                prob = remake(prob, p=tltl_p)
                push!(sol, DifferentialEquations.solve(prob, Euler(); dt=0.015625, dtmax=0.015625))
            else
                tltl_p = default_tltl_pars(prob.p)
                if (α > 1 && tltl_p[p] + δ * α_values[α-1] < 1.0)
                    tltl_p[p] = 0.99999999
                    prob = remake(prob, p=tltl_p)
                    push!(sol, DifferentialEquations.solve(prob, Euler(); dt=0.015625, dtmax=0.015625))
                end
            end
        end
        perc_p = []
        for v in [e4a.AWBI, e4a.GDPP, e4a.INEQ, e4a.OW, e4a.POP, e4a.STE]
            perc_v = []
            for i in 1:lastindex(sol)
                li = lastindex(sol[i][v])
                push!(perc_v, (sol[i][v][li] - tltl_sol[v][li]) / tltl_sol[v][li])
            end
            push!(perc_p, perc_v)
        end
        push!(perc, perc_p)
    end
    return perc
end

function plot_perc(perc, p, sl)
    gl_p = default_gl_pars(prob.p)
    tltl_p = default_tltl_pars(prob.p)
    δ = gl_p[p] - tltl_p[p]
    println(δ)
    α_values = [0.0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.0]
    x_labels = []
    for α in 1:lastindex(α_values)
        if (tltl_p[p] + δ * α_values[α] <= 1.0 || !(p in [4, 10, 11, 13, 14, 16, 17, 21]))
            push!(x_labels, tltl_p[p] + δ * α_values[α])
        else
            if (α > 1 && tltl_p[p] + δ * α_values[α-1] < 1.0)
                push!(x_labels, 1.0)
            end
        end
    end
    println(x_labels)
    data = GenericTrace[]
    for v in 1:lastindex(perc[p])
        push!(data, scatter(; x=x_labels, y=perc[p][v], name=v_name[v], mode="lines+markers"))
    end
    layout = Layout(; title=p_desc[p], yaxis_range=[-1.0, 2.0], showlegend=sl)
    PlotlyJS.plot(data, layout)
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

function plot_variables(solution, xrange, variables::Vector{<:NTuple{4,Any}}; title="", showaxis=true, showlegend=true, linetype="lines", colored=true, save=false)
    numvars = length(variables)

    @assert 1 ≤ numvars
    @assert (1 == length(xrange)) || (3 == length(xrange))
    @assert 4 == length(variables[1])


    colors = colored ? ColorSchemes.tab10.colors : fill(RGB(0.2, 0.2, 0.2), numvars)

    x_offset = 0.05
    x_domain = showaxis ? x_offset * numvars - 0.04 : 0.0


    traces = GenericTrace[]

    (xvalue, xmin, xmax) = (3 == length(xrange)) ? xrange : (xrange[1], -Inf, Inf)
    (var, varmin, varmax, varname) = variables[1]

    layout = Dict([
        ("title", attr(text=title, x=0.5)),
        ("showlegend", showlegend),
        ("plot_bgcolor", "#EEE"),
        ("xaxis", attr(
            domain=[x_domain + 0.02, 1.0],
            position=0.0,
            range=[xmin, xmax])),
        ("yaxis", attr(
            color=colors[1],
            visible=showaxis,
            name="",
            position=0.0,
            showgrid=false,
            range=[varmin, varmax],
            domain=[0.05, 1.0],
        ))
    ])

    push!(traces, scatter(
        x=solution[xvalue],
        y=solution[var],
        marker_color=colors[1],
        name=varname,
        mode=linetype, yaxis="y1"),
    )


    for i ∈ 2:numvars
        (var, varmin, varmax, varname) = variables[i]

        layout[string("yaxis", i)] = attr(
            color=colors[i],
            overlaying="y",
            visible=showaxis,
            name="",
            position=(i - 1) * x_offset,
            showgrid=false,
            range=[varmin, varmax],
        )

        push!(traces, scatter(
            x=solution[xvalue],
            y=solution[var],
            marker_color=colors[i],
            name=varname,
            mode=linetype,
            yaxis=string("y", i)),
        )
    end

    p = plot(traces, Layout(layout))
    save && savefig(p, "./" * title * ".svg")

    return p
end

function plot_two_sol_variables(scen1, sol1, scen2, sol2, xrange, variables::Vector{<:NTuple{4,Any}}; title="", showaxis=true, showlegend=true, linetype="lines", colored=true, save=false)
    numvars = length(variables)

    @assert 1 ≤ numvars
    @assert (1 == length(xrange)) || (3 == length(xrange))
    @assert 4 == length(variables[1])


    colors = colored ? ColorSchemes.tab10.colors : fill(RGB(0.2, 0.2, 0.2), numvars)

    x_offset = 0.05
    x_domain = showaxis ? x_offset * numvars - 0.04 : 0.0


    traces = GenericTrace[]

    (xvalue, xmin, xmax) = (3 == length(xrange)) ? xrange : (xrange[1], -Inf, Inf)
    (var, varmin, varmax, varname) = variables[1]

    layout = Dict([
        ("title", attr(text=title, x=0.5)),
        ("showlegend", showlegend),
        ("plot_bgcolor", "#EEE"),
        ("xaxis", attr(
            domain=[x_domain + 0.02, 1.0],
            position=0.0,
            range=[xmin, xmax])),
        ("yaxis", attr(
            color=colors[1],
            visible=showaxis,
            name="",
            position=0.0,
            showgrid=false,
            range=[varmin, varmax],
            domain=[0.05, 1.0],
        ))
    ])

    push!(traces, scatter(
        x=sol1[xvalue],
        y=sol1[var],
        marker_color=colors[1],
        name=varname * " (" * scen1 * ")",
        mode=linetype, yaxis="y1"),
    )

    push!(traces, scatter(
        x=sol2[xvalue],
        y=sol2[var],
        marker_color=colors[1],
        name=varname * " (" * scen2 * ")",
        mode=linetype, yaxis="y1",
        line=attr(dash="dash"))
    )

    for i ∈ 2:numvars
        (var, varmin, varmax, varname) = variables[i]

        layout[string("yaxis", i)] = attr(
            color=colors[i],
            overlaying="y",
            visible=showaxis,
            name="",
            position=(i - 1) * x_offset,
            showgrid=false,
            range=[varmin, varmax],
        )

        push!(traces, scatter(
            x=sol1[xvalue],
            y=sol1[var],
            marker_color=colors[i],
            name=varname * " (" * scen1 * ")",
            mode=linetype,
            yaxis=string("y", i)),
        )

        push!(traces, scatter(
            x=sol2[xvalue],
            y=sol2[var],
            marker_color=colors[i],
            name=varname * " (" * scen2 * ")",
            mode=linetype,
            yaxis=string("y", i),
            line=attr(dash="dash"))
        )
    end

    p = plot(traces, Layout(layout))
    save && savefig(p, "./" * title * ".svg")

    return p
end

function plot_sol_var(e4a, _sol, var_ind, _title, sa, sl; kwargs...)
    @variables t
    return plot_variables(_sol, (t, 1980, 2100), getindex(_variables(e4a), var_ind); title=_title, showaxis=sa, showlegend=sl, kwargs...)
end

function plot_two_sol_var(e4a, _scen1, _sol1, _scen2, _sol2, var_ind, _title, sa, sl; kwargs...)
    @variables t
    return plot_two_sol_variables(_scen1, _sol1, _scen2, _sol2, (t, 1980, 2100), getindex(_variables(e4a), var_ind); title=_title, showaxis=sa, showlegend=sl, kwargs...)
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

function read_vensim_dataset(fn, to_be_removed)
    f::IOStream = open(fn, "r")
    ds = Dict{String,Array{Float64}}()
    for line in eachline(f)
        split_line::Vector{String} = split(line, "\t")
        s = replace(split_line[1], "\"" => "")
        s = replace(s, to_be_removed => "")
        s = replace(s, " (Year)" => "")
        s = replace(s, "\$" => "dollar")
        v::Array{Float64} = Array{Float64}(undef, length(split_line) - 1)
        for i in 2:lastindex(split_line)
            v[i-1] = parse(Float64, split_line[i])
        end
        ds[lowercase(s)] = v
    end
    close(f)
    return ds
end

function compare(a, b, pepsi)
    max_re = -1
    max_re_a = 0
    max_re_b = 0
    max_re_i = 0
    for i in 1:lastindex(a)
        if (pepsi >= 0)
            re = abs(a[i] - b[i]) / (abs(b[i]) + pepsi)
        else
            re = 0
            if (a[i] != 0 || b[i] != 0)
                re = (2 * abs(a[i] - b[i])) / (abs(a[i]) + abs(b[i]))
                # re = abs(a[i] - b[i]) / max(abs(a[i]), abs(b[i]))
            end
        end
        if (re > max_re)
            max_re = max(max_re, re)
            max_re_a = a[i]
            max_re_b = b[i]
            max_re_i = i
        end
    end
    return max_re, max_re_a, max_re_b, max_re_i
end

function mre_sys(sol, sys, vs_ds, pepsi, nt, verbose)
    max_re = 0.0
    for v in states(sys)
        d = getdescription(v)
        if (d != "" && d != "Time instants" && !startswith(d, "LV functions") && !startswith(d, "RT functions"))
            if (get(vs_ds, lowercase(d), "") != "")
                re, _, _, _ = Earth4All.compare(sol[v][1:nt], vs_ds[lowercase(d)], pepsi)
                max_re = max(max_re, re)
                if (verbose)
                    println(d, "\t", re)
                end
            end
        end
    end
    return max_re
end

function all_mre(scen, sector, sys, sol)
    max_re = 0
    # sn = ["climate", "demand", "energy", "finance", "foodland", "inventory", "labourmarket", "other", "output", "population", "public", "wellbeing"]
    println("====" * uppercase(sector) * "====")
    vs_ds = read_vensim_dataset("vensim/" * lowercase(scen) * "/" * sector * ".txt", " : E4A-220501 " * scen)
    re = mre_sys(sol, sys, vs_ds, 1, 7681, true)
    max_re = max(max_re, re)
    println("====MAXIMUM ERROR===" * string(max_re))
end

# function sobol_bar_plots(e4a, res)
#     fp = "plots/sobol/"
#     mkpath(fp)
#     v = ["Average wellbeing index", "GDP per person", "Inequality", "Observed warming", "Population", "Social tension"]
#     a = ["awbi", "gdpp", "ineq", "ow", "pop", "ste"]
#     for i in 1:lastindex(v)
#         p = plot([bar(x=getdescription.(parameters(e4a)), y=res.S1[2*i-1, :], name="Sobol first index"), bar(x=getdescription.(parameters(e4a)), y=res.ST[2*i-1, :], name="Sobol total index")], Layout(title=v[i] * " (averaged over 6 years)"))
#         savefig(p, fp * a[i] * "_06.html")
#         savefig(p, fp * a[i] * "_06.png")
#         p = plot([bar(x=getdescription.(parameters(e4a)), y=res.S1[2*i, :], name="Sobol first index"), bar(x=getdescription.(parameters(e4a)), y=res.ST[2*i, :], name="Sobol total index")], Layout(title=v[i] * " (averaged over 12 years)"))
#         savefig(p, fp * a[i] * "_12.html")
#         savefig(p, fp * a[i] * "_12.png")
#     end
# end

function sobol_bar_plots(res)
    fp = "plots/sobol/"
    mkpath(fp)
    v = ["Average wellbeing index", "GDP per person", "Inequality", "Observed warming", "Population", "Social tension"]
    a = ["awbi", "gdpp", "ineq", "ow", "pop", "ste"]
    for i in 1:lastindex(v)
        p = plot([bar(x=p_name[p_map], y=res.S1[2*i-1, p_map], name="Sobol first index"), bar(x=p_name[p_map], y=res.ST[2*i-1, p_map], name="Sobol total index")], Layout(title=v[i] * " (averaged over 6 years)"))
        savefig(p, fp * a[i] * "_06.html")
        savefig(p, fp * a[i] * "_06.png")
        p = plot([bar(x=p_name[p_map], y=res.S1[2*i, p_map], name="Sobol first index"), bar(x=p_name[p_map], y=res.ST[2*i, p_map], name="Sobol total index")], Layout(title=v[i] * " (averaged over 12 years)"))
        savefig(p, fp * a[i] * "_12.html")
        savefig(p, fp * a[i] * "_12.png")
    end
end
