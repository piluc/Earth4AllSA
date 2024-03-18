include("Earth4All.jl")

using ColorTypes
using ColorSchemes
using GlobalSensitivity
using ModelingToolkit
using DifferentialEquations
using PlotlyJS
using Serialization
using Statistics
using WorldDynamics

global p_name = ["ERDN2OKF2022", "ERDCH4KC2022", "DACCO22100", "FGDC2022", "EETF2022", "EGTRF2022", "EPTF2022", "ETGBW", "GEIC", "FETPO", "GFCO2SCCS", "EROCEPA2022", "GFNE", "GREF", "GCWR", "GFNRM", "GFRA", "EROCFSP", "USPIS2022", "USPUS2022", "GEFR", "MIROTA2022"]
global p_desc = ["Extra rate of decline in N2O per kg fertilizer from 2022", "Extra rate of decline in CH4 per kg crop after 2022 1/y", "Direct air capture of CO2 in 2100 GtCO2/y", "Fraction of government debt cancelled in 2022 1/y", "Extra empowerment tax from 2022 (share of NI)", "Extra general tax rate from 2022", "Extra pension tax from 2022 (share of NI)", "Extra transfer of government budget to workers", "Goal for extra income from commons (share of NI)", "Fraction of extra taxes paid by owners", "Goal for fraction of CO2-sources with CCS", "Extra ROC in energy productivity after 2022 1/y", "Goal for fraction new electrification", "Goal for renewable el fraction", "Goal for crop waste reduction", "Goal for fraction new red meat", "Goal for fraction regenerative agriculture", "Extra ROC in food sector productivity from 2022 1/y", "Unconventional stimulus in PIS from 2022 (share of GDP)", "Unconventional stimulus in PUS from 2022 (share of GDP)", "Goal for extra fertility reduction", "Max imported ROTA from 2022 1/y"]
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

function default_tltl_pars(x)
    y = copy(x)
    # Poverty turnaround
    y[4] = 0.0
    y[19] = 0.0
    y[20] = 0.0
    y[22] = 0.0
    # Inequality turnaround
    y[6] = 0.0
    y[10] = 0.5
    y[8] = 0.0
    y[9] = 0.0
    # Empowerment turnaround
    y[21] = 0.0
    y[5] = 0.0
    y[7] = 0.00
    # Food turnaround
    y[18] = 0.0
    y[15] = 0.05
    y[17] = 0.1
    y[16] = 0.1
    # Energy turnaround
    y[12] = 0.002
    y[13] = 0.5
    y[14] = 0.5
    y[11] = 0.2
    y[3] = 0.0
    # Other
    y[2] = 0.0
    y[1] = 0.0
    return y
end

function default_gl_pars(x)
    y = copy(x)
    # Poverty turnaround
    y[4] = 0.1
    y[19] = 0.01
    y[20] = 0.01
    y[22] = 0.005
    # Inequality turnaround
    y[6] = 0.01
    y[10] = 0.8
    y[8] = 0.2
    y[9] = 0.02
    # Empowerment turnaround
    y[21] = 0.2
    y[5] = 0.02
    y[7] = 0.02
    # Food turnaround
    y[18] = 0.0
    y[15] = 0.2
    y[17] = 0.5
    y[16] = 0.5
    # Energy turnaround
    y[12] = 0.004
    y[13] = 1.0
    y[14] = 1.0
    y[11] = 0.9
    y[3] = 8.0
    # Other
    y[2] = 0.01
    y[1] = 0.01
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

function plot_sol_var(e4a, _sol, var_ind, _title, sa, sl; kwargs...)
    @variables t
    return plot_variables(_sol, (t, 1980, 2100), getindex(_variables(e4a), var_ind); title=_title, showaxis=sa, showlegend=sl, kwargs...)
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

function plot_two_sol_var(e4a, _scen1, _sol1, _scen2, _sol2, var_ind, _title, sa, sl; kwargs...)
    @variables t
    return plot_two_sol_variables(_scen1, _sol1, _scen2, _sol2, (t, 1980, 2100), getindex(_variables(e4a), var_ind); title=_title, showaxis=sa, showlegend=sl, kwargs...)
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

function compute_gl18_sol(prob)
    a = fill(1.0, 22)
    a[1] = 0.0
    a[12] = 0.0
    a[15] = 0.0
    return compute_alpha_sol(prob, [], a, collect(1:22))
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

function retrieve_sobol(fn)
    return deserialize(fn)
end

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

function tornado_diagram(sys, prob, base_sol, var, delta_t)
    gl_p = default_gl_pars(prob.p)
    tltl_p = default_tltl_pars(prob.p)
    pars = parameters(sys)
    np = length(pars)
    bb = mean(base_sol[var][end-delta_t:end])
    lb = fill(0.0, np)
    ub = fill(0.0, np)
    for p in 1:np
        δ = gl_p[p] - tltl_p[p]
        if (δ > 0)
            a = fill(1.0, np)
            a[p] = 0.0
            lb_sol = compute_alpha_sol(prob, [], a, collect(1:np))
            lb[p] = (mean(lb_sol[var][end-delta_t:end]) - bb) / bb
            a = fill(1.0, np)
            if (tltl_p[p] + 2.0 * δ <= 1.0 || !(p in [4, 10, 11, 13, 14, 16, 17, 21]))
                a[p] = 2.0
            else
                a[p] = (0.99999999 - tltl_p[p]) / δ
            end
            ub_sol = compute_alpha_sol(prob, [], a, collect(1:np))
            ub[p] = (mean(ub_sol[var][end-delta_t:end]) - bb) / bb
        else
            lb[p] = 0.0
            ub[p] = 0.0
        end
    end
    td = []
    db = abs.(ub - lb)
    ldp = sortperm(db, rev=true)
    for p in 1:np
        push!(td, (ldp[p], pars[ldp[p]], round(lb[ldp[p]]; digits=3) * 100, round(ub[ldp[p]]; digits=3) * 100))
    end
    return td
end

function local_sensitivity_gl(e4a, prob, base_sol)
    α_values = [0.0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.0]
    gl_p = default_gl_pars(prob.p)
    tltl_p = default_tltl_pars(prob.p)
    prob = remake(prob, p=tltl_p)
    # tltl_sol = DifferentialEquations.solve(prob, Euler(); dt=0.015625, dtmax=0.015625)
    perc = []
    for p in 1:lastindex(gl_p)
        println("Computing percentages for parameter ", p)
        δ = gl_p[p] - tltl_p[p]
        sol = []
        for α in 1:lastindex(α_values)
            gl_p = default_gl_pars(prob.p)
            gl_p[p] = tltl_p[p] + δ * α_values[α]
            if (gl_p[p] <= 1.0 || !(p in [4, 10, 11, 13, 14, 16, 17, 21]))
                prob = remake(prob, p=gl_p)
                push!(sol, DifferentialEquations.solve(prob, Euler(); dt=0.015625, dtmax=0.015625))
            else
                gl_p = default_gl_pars(prob.p)
                if (α > 1 && tltl_p[p] + δ * α_values[α-1] < 1.0)
                    gl_p[p] = 0.99999999
                    prob = remake(prob, p=gl_p)
                    push!(sol, DifferentialEquations.solve(prob, Euler(); dt=0.015625, dtmax=0.015625))
                end
            end
        end
        perc_p = []
        for v in [e4a.AWBI, e4a.GDPP, e4a.INEQ, e4a.OW, e4a.POP, e4a.STE]
            perc_v = []
            for i in 1:lastindex(sol)
                li = lastindex(sol[i][v])
                push!(perc_v, (sol[i][v][li] - base_sol[v][li]) / base_sol[v][li])
            end
            push!(perc_p, perc_v)
        end
        push!(perc, perc_p)
    end
    return perc
end

function spider_plot(perc, sys, v, pars, tt, sl, acron)
    α_values = [0.0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.0]
    data = GenericTrace[]
    for p in pars
        scatter_name = p_desc[p]
        if (acron)
            scatter_name = p_name[p]
        end
        push!(data, scatter(; x=α_values, y=perc[p][v], name=scatter_name, mode="lines+markers"))
    end
    plot_title = ""
    if (tt)
        plot_title = v_name[v]
    end
    layout = Layout(; title=plot_title, xaxis=attr(tickmode="array", tickvals=[0.0, 1.0, 2.0], ticktext=["Low=TLTL", "Base=GL", "High"]), yaxis_range=[-0.25, 0.25], yaxis=attr(tickmode="array", tickvals=[-0.2, -0.1, 0.0, 0.1, 0.2], ticktext=["-20%", "-10%", "0%", "10%", "20%"]), showlegend=sl)
    PlotlyJS.plot(data, layout)
end

function mre_two_sol(sys, sol1, sol2, pepsi, verbose)
    max_re = 0.0
    max_v = nothing
    for v in states(sys)
        d = getdescription(v)
        if (d != "" && d != "Time instants" && !startswith(d, "LV functions") && !startswith(d, "RT functions"))
            re, _, _, _ = Earth4All.compare(sol1[v], sol2[v], pepsi)
            if (re > max_re)
                max_re = max(max_re, re)
                max_v = v
                if (verbose)
                    println(d, "\t", re)
                end
            end
        end
    end
    return max_v, max_re
end

function plot_two_sols(scen1, sol1, scen2, sol2, vars, plot_title)
    x = range(1, 7681, length=7681)
    traces = GenericTrace[]
    for v in 1:lastindex(vars)
        desc = getdescription(vars[v])
        trace1 = PlotlyJS.scatter(x=x, y=sol1[vars[v]], name=desc * "-" * scen1, line=attr(color="royalblue", dash="dash"))
        trace2 = PlotlyJS.scatter(x=x, y=sol2[vars[v]], name=desc * "-" * scen2, line=attr(color="firebrick", dash="dot"))
        push!(traces, trace1)
        push!(traces, trace2)
    end
    return PlotlyJS.plot(traces, Layout(legend=attr(x=0, y=1,), title=plot_title, xaxis=attr(tickmode="array", tickvals=collect(1:320:7681), ticktext=string.(collect(1980:5:2100)))))
end

