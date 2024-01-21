include("Earth4All.jl")

using ModelingToolkit
using DifferentialEquations
using Plots
using Serialization
using Statistics

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

function local_sensitivity_gl()
    @named e4a = Earth4All.earth4all()
    e4a_sys = structural_simplify(e4a)
    prob = ODEProblem(e4a_sys, [], (1980, 2100))
    gl_p = default_gl_pars(prob.p)
    tltl_p = default_tltl_pars(prob.p)
    prob = remake(prob, p=gl_p)
    gl_sol = solve(prob, Euler(); dt=0.015625, dtmax=0.015625)
    perc = []
    for p in 1:lastindex(gl_p)
        println("Computing percentages for parameter ", p)
        δ = gl_p[p] - tltl_p[p]
        sol = []
        for α in [0.0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.0]
            gl_p = default_gl_pars(prob.p)
            gl_p[p] = gl_p[p] + δ * α
            prob = remake(prob, p=gl_p)
            push!(sol, solve(prob, Euler(); dt=0.015625, dtmax=0.015625))
        end
        perc_p = []
        for v in [e4a.AWBI, e4a.GDPP, e4a.INEQ, e4a.OW, e4a.POP, e4a.STE]
            perc_v = []
            for i in 1:lastindex(sol)
                push!(perc_v, (sol[i][v][lastindex(sol[1][v])] - gl_sol[v][lastindex(sol[1][v])]) / gl_sol[v][lastindex(sol[1][v])])
            end
            push!(perc_p, perc_v)
        end
        push!(perc, perc_p)
    end
    return gl_sol, sol, perc
end

function local_sensitivity_tltl()
    @named e4a = Earth4All.earth4all()
    e4a_sys = structural_simplify(e4a)
    prob = ODEProblem(e4a_sys, [], (1980, 2100))
    gl_p = default_gl_pars(prob.p)
    tltl_p = default_tltl_pars(prob.p)
    prob = remake(prob, p=tltl_p)
    tltl_sol = solve(prob, Euler(); dt=0.015625, dtmax=0.015625)
    perc = []
    for p in 1:lastindex(gl_p)
        println("Computing percentages for parameter ", p)
        δ = gl_p[p] - tltl_p[p]
        sol = []
        for α in [0.0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.0]
            tltl_p = default_tltl_pars(prob.p)
            tltl_p[p] = tltl_p[p] + δ * α
            prob = remake(prob, p=tltl_p)
            push!(sol, solve(prob, Euler(); dt=0.015625, dtmax=0.015625))
        end
        perc_p = []
        for v in [e4a.AWBI, e4a.GDPP, e4a.INEQ, e4a.OW, e4a.POP, e4a.STE]
            perc_v = []
            for i in 1:lastindex(sol)
                push!(perc_v, (sol[i][v][lastindex(sol[1][v])] - tltl_sol[v][lastindex(sol[1][v])]) / tltl_sol[v][lastindex(sol[1][v])])
            end
            push!(perc_p, perc_v)
        end
        push!(perc, perc_p)
    end
    return tltl_sol, sol, perc
end

function local_sensitivity_tltl_all()
    @named e4a = Earth4All.earth4all()
    e4a_sys = structural_simplify(e4a)
    prob = ODEProblem(e4a_sys, [], (1980, 2100))
    gl_p = default_gl_pars(prob.p)
    tltl_p = default_tltl_pars(prob.p)
    prob = remake(prob, p=tltl_p)
    tltl_sol = solve(prob, Euler(); dt=0.015625, dtmax=0.015625)
    δ = gl_p - tltl_p
    sol = []
    for α in [0.0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.0]
        println("Computing percentages for alpha=", α)
        tltl_p = default_tltl_pars(prob.p)
        tltl_p = tltl_p + δ * α
        prob = remake(prob, p=tltl_p)
        push!(sol, solve(prob, Euler(); dt=0.015625, dtmax=0.015625))
    end
    perc = []
    for v in [e4a.AWBI, e4a.GDPP, e4a.INEQ, e4a.OW, e4a.POP, e4a.STE]
        perc_v = []
        for i in 1:lastindex(sol)
            push!(perc_v, (sol[i][v][lastindex(sol[1][v])] - tltl_sol[v][lastindex(sol[1][v])]) / tltl_sol[v][lastindex(sol[1][v])])
        end
        push!(perc, perc_v)
    end
    return tltl_sol, sol, perc
end

function local_sensitivity(p1, p2, p3, p4, p5, p6, p7)
    α_values = [0.0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.0]
    α_values = [0.0, 0.33, 0.66, 1, 1.33, 1.66, 2.0]
    @named e4a = Earth4All.earth4all()
    e4a_sys = structural_simplify(e4a)
    prob = ODEProblem(e4a_sys, [], (1980, 2100))
    gl_p = default_gl_pars(prob.p)
    tltl_p = default_tltl_pars(prob.p)
    prob = remake(prob, p=tltl_p)
    @time tltl_sol = solve(prob, Euler(); dt=0.015625, dtmax=0.015625)
    sol = []
    δ = gl_p - tltl_p
    for α1 in 1:lastindex(α_values)
        for α2 in 1:lastindex(α_values)
            for α3 in 1:lastindex(α_values)
                for α4 in 1:lastindex(α_values)
                    for α5 in 1:lastindex(α_values)
                        for α6 in 1:lastindex(α_values)
                            for α7 in 1:lastindex(α_values)
                                # for α8 in 1:lastindex(α_values)
                                #     for α9 in 1:lastindex(α_values)
                                st = time()
                                print(α1, α2, α3, α4, α5, α6, α7)
                                tltl_p = default_tltl_pars(prob.p)
                                tltl_p[p1] = tltl_p[p1] + δ[p1] * α_values[α1]
                                tltl_p[p2] = tltl_p[p2] + δ[p2] * α_values[α2]
                                tltl_p[p3] = tltl_p[p3] + δ[p3] * α_values[α3]
                                tltl_p[p4] = tltl_p[p4] + δ[p4] * α_values[α4]
                                tltl_p[p5] = tltl_p[p5] + δ[p5] * α_values[α5]
                                tltl_p[p6] = tltl_p[p6] + δ[p6] * α_values[α6]
                                tltl_p[p7] = tltl_p[p7] + δ[p7] * α_values[α7]
                                # tltl_p[p8] = tltl_p[p8] + δ[p8] * α_values[α8]
                                # tltl_p[p9] = tltl_p[p9] + δ[p9] * α_values[α9]
                                prob = remake(prob, p=tltl_p)
                                sol_α = solve(prob, Euler(); dt=0.015625, dtmax=0.015625)
                                li = lastindex(sol_α.t)
                                sol_v = [sol_α[e4a.AWBI][li], sol_α[e4a.GDPP][li], sol_α[e4a.INEQ][li], sol_α[e4a.OW][li], sol_α[e4a.POP][li], sol_α[e4a.STE][li]]
                                push!(sol, sol_v)
                                println(": ", time() - st, " seconds")
                                #     end
                                # end
                            end
                        end
                    end
                end
            end
        end
    end
    serialize("sols/sol_7_7.dat", sol)
end

local_sensitivity(3, 9, 11, 19, 20, 21, 22)
