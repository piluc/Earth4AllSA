using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("Earth4All.jl")

# using ModelingToolkit
# sol = Earth4All.tltl_solution()
# vs_ds = Earth4All.read_vensim_dataset("e4a/vensim/tltl/population.txt", " : E4A-220501 TLTL")
# @named e4a = Earth4All.earth4all()
# Earth4All.compare_and_plot("TLTL", sol, "Population Mp", e4a.POP, vs_ds, 1980, 2100, 7681, 1, true)
