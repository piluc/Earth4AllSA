### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 88115292-3f5b-405b-b611-2508d845b83d
begin
	using Pkg;
	Pkg.activate("Project.toml");
	Pkg.instantiate();
	using PlutoUI;
	println("Julia environment has been set up")
end

# ╔═╡ 6e2f37d1-7870-40b3-9e1a-d23eabb7b780
begin
	include("src/E4ASA.jl");
	println("E4ASA code has been included")
end

# ╔═╡ cafed8ec-e60d-11ee-1088-af7fa6542847
md"""
# Making the World Better with Fewer Turnarounds
## A sensitivity analysis of the Earth for All model
"""

# ╔═╡ b4501aec-6fb1-483a-9edf-3a9c45111cc1
	md"""
	Let's start by setting the Julia environment...
	"""

# ╔═╡ 8d5ac7c7-44ce-4ed9-a615-fa45b27aeb9e
md"""
...and by including the Julia code.
"""

# ╔═╡ d8aee08b-16b0-406c-a6d9-5392a333c488
md"""
We can now compute the seven scenarios analyzed in the paper.
"""

# ╔═╡ eeae7635-9488-476e-95a7-c6960cee892c
begin
	e4a, _, prob = e4a_sys_prob()
	tltl = compute_tltl_sol(prob)
	gl = compute_gl_sol(prob)
	gl18 = compute_gl18_sol(prob);
	gl6 = compute_alpha_sol(prob, [], [1.0,1.0,1.0,1.0,1.0,1.0], [3,9,11,19,20,22])
	dgl6 = compute_alpha_sol(prob, [], [2.0,2.0,2.0,2.0,2.0,2.0], [3,9,11,19,20,22])
	gl4 = compute_alpha_sol(prob, [], [0.25,2.0,0.5,2.0,2.0,2.0], [3,9,11,19,20,22])
	a = fill(1.0, 22)
	a[3] = 0.25
	a[11] = 0.5
	mgl = compute_alpha_sol(prob, [], a, collect(1:22))
	println("All scenarios have been computed")
end

# ╔═╡ 8284833a-0e0f-4867-b2b3-098470774562
md"""
Below we can select two scenarios (among the computed seven ones) and one output indicator, in order to compare the evolution of the indicator in the two scenarios. 
"""

# ╔═╡ 8c9cba78-e868-44a2-8663-8aa51929b4f5
md"""
First scenario $(@bind scen1 Select(["TLTL","GL","GL18","GL6","DGL6","GL4","MGL"]))
Second scenario $(@bind scen2 Select(["TLTL","GL","GL18","GL6","DGL6","GL4","MGL"]))
Indicator $(@bind indicator Select([e4a.AWBI,e4a.GDPP,e4a.INEQ,e4a.OW,e4a.POP,e4a.STE]))
"""

# ╔═╡ 2ee8fc9b-775c-4cae-b63b-23ccad3ff114
begin
	using Plots
	plotly()
	scens = [tltl,gl,gl18,gl6,dgl6,gl4,mgl]
	scen_names = ["TLTL","GL","GL18","GL6","DGL6","GL4","MGL"]
	i1 = findfirst(s -> s==scen1, scen_names)
	i2 = findfirst(s -> s==scen2, scen_names)
	plot_two_sols(scen1, scens[i1], scen2, scens[i2], [indicator], scen1 * " versus " * scen2)
end

# ╔═╡ 151c4b95-f492-4753-8c29-664606a74ccf
md"""
In order to compare the evolution of the other variables of the model, we can specify the acronym of the variable in the code below prefized by `e4a`.
"""

# ╔═╡ cadaad8a-452d-4489-95b5-c9f6f6e4c044
plot_two_sols(scen1, scens[i1], scen2, scens[i2], [e4a.OF], scen1 * " versus " * scen2)

# ╔═╡ a16f3874-4be5-48b6-a70a-26e8426603ab
md"""
Finally, below we can set the value of the twenty-one parameters involved in the five turnarounds and (by specifying its acronym prefixed by `e4a`) see the evolution of any variable in the TLTL scenario, in the GL scenario and in the scenario (called Pluto) produced by the chosen set of values.
"""

# ╔═╡ 0cb5be58-1df8-4266-8c83-88c3f85b92be
begin
struct MySlider 
    range::AbstractRange
    default::Number
end
function Base.show(io::IO, ::MIME"text/html", slider::MySlider)
    print(io, """
		<input type="range" 
		min="$(first(slider.range))" 
		step="$(step(slider.range))"
		max="$(last(slider.range))" 
		value="$(slider.default)"
		oninput="this.nextElementSibling.value=this.value">
		<output>$(slider.default)</output>""")
end
end

# ╔═╡ a4d3e6c6-10d7-4fde-aa86-431d5a54e95b
begin
	md"Direct air capture of CO2 in 2100 GtCO2/y $(@bind DACCO22100 MySlider(0.0:0.5:16.0,8.0))\
	Extra empowerment tax from 2022 $(@bind EETF2022 MySlider(0.0:0.001:0.04,0.02))\
	Extra general tax rate from 2022 $(@bind EGTRF2022 MySlider(0.0:0.001:0.02,0.01))\
	Extra pension tax from 2022 (share of NI) $(@bind EPTF2022 MySlider(0.0:0.001:0.04,0.02))\
	Extra rate of decline in CH4 per kg crop after 2022 1/y $(@bind ERDCH4KC2022 MySlider(0.0:0.001:0.02,0.01))\
	Extra rate of decline in N2O per kg fertilizer from 2022 $(@bind ERDN2OKF2022 MySlider(0.0:0.001:0.02,0.01))\
	Extra ROC in energy productivity after 2022 1/y $(@bind EROCEPA2022 MySlider(0.002:0.0001:0.006,0.004))\
	Extra transfer of govmnt budget to workers $(@bind ETGBW MySlider(0.0:0.01:0.4,0.2))\
	Fraction of extra taxes paid by owners $(@bind FETPO MySlider(0.5:0.01:0.99,0.8))\
	Fraction of govmnt debt cancelled in 2022 1/y $(@bind FGDC2022 MySlider(0.0:0.01:0.2,0.1))\
	Goal for crop waste reduction $(@bind GCWR MySlider(0.05:0.01:0.35,0.2))\
	Goal for extra fertility reduction $(@bind GEFR MySlider(0.0:0.01:0.4,0.2))\
	Goal for extra income from commons (share of NI) $(@bind GEIC MySlider(0.0:0.001:0.04,0.02))\
	Goal for fraction of CO2-sources with CCS $(@bind GFCO2SCCS MySlider(0.2:0.01:0.99,0.9))\
	Goal for fraction new electrification $(@bind GFNE MySlider(0.5:0.01:1.0,1.0))\
	Goal for fraction new red meat $(@bind GFNRM MySlider(0.1:0.01:0.9,0.5))\
	Goal for fraction regenerative agriculture $(@bind GFRA MySlider(0.1:0.01:0.9,0.5))\
	Goal for renewable el fraction $(@bind GREF MySlider(0.5:0.01:1.0,1.0))\
	Max imported ROTA from 2022 1/y $(@bind MIROTA2022 MySlider(0.0:0.0001:0.01,0.005))\
	Unconventional stimulus in PIS from 2022 (share of GDP) $(@bind USPIS2022 MySlider(0.0:0.001:0.02,0.01))\
	Unconventional stimulus in PUS from 2022 (share of GDP) $(@bind USPUS2022 MySlider(0.0:0.001:0.02,0.01))\
	"
end

# ╔═╡ 21abb389-f124-493b-8eda-251978ec2420
sol = compute_sol(prob, [DACCO22100,EETF2022,EGTRF2022,EPTF2022,ERDCH4KC2022,ERDN2OKF2022,EROCEPA2022,ETGBW,FETPO,FGDC2022,GCWR,GEFR,GEIC,GFCO2SCCS,GFNE,GFNRM,GFRA,GREF,MIROTA2022,USPIS2022,USPUS2022]);

# ╔═╡ db64d3d7-3f0d-4596-b159-11f2257e79c4
plot_three_sols("TLTL", tltl, "GL", gl, "Pluto", sol, [e4a.STE], "Pluto versus TLTL and GL")

# ╔═╡ Cell order:
# ╟─cafed8ec-e60d-11ee-1088-af7fa6542847
# ╟─b4501aec-6fb1-483a-9edf-3a9c45111cc1
# ╟─88115292-3f5b-405b-b611-2508d845b83d
# ╟─8d5ac7c7-44ce-4ed9-a615-fa45b27aeb9e
# ╟─6e2f37d1-7870-40b3-9e1a-d23eabb7b780
# ╟─d8aee08b-16b0-406c-a6d9-5392a333c488
# ╟─eeae7635-9488-476e-95a7-c6960cee892c
# ╟─8284833a-0e0f-4867-b2b3-098470774562
# ╟─8c9cba78-e868-44a2-8663-8aa51929b4f5
# ╟─2ee8fc9b-775c-4cae-b63b-23ccad3ff114
# ╟─151c4b95-f492-4753-8c29-664606a74ccf
# ╠═cadaad8a-452d-4489-95b5-c9f6f6e4c044
# ╟─a16f3874-4be5-48b6-a70a-26e8426603ab
# ╟─0cb5be58-1df8-4266-8c83-88c3f85b92be
# ╟─a4d3e6c6-10d7-4fde-aa86-431d5a54e95b
# ╟─21abb389-f124-493b-8eda-251978ec2420
# ╠═db64d3d7-3f0d-4596-b159-11f2257e79c4
