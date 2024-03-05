# Earth4AllSA
Code associated to the paper *Making the Word Better with Fewer Turnarounds. A sensitivity analysis of the Earth for All model* and based on the Julia implementation of the [Earth4All model](https://earth4all.life/the-science-rp/) using the [WorldDynamics framework](https://github.com/worlddynamics/WorldDynamics.jl) and referring to the [April 2023 Vensim version](https://web.archive.org/web/20220830093115/https://earth4all.life/the-science/) available [here](https://stockholmuniversity.app.box.com/s/uh7fjh52pvh7yx1mqfwqcyxdcvegrodf/folder/170558692760).

## How to run the code

### Prerequisites

[Install Julia](https://julialang.org/) and clone the repository: 
```sh
git clone https://github.com/piluc/Earth4AllSA
```

### Setting up the environment

After starting the Julia REPL in the repository folder, we can instantiate the environment by running
```jl
julia> using Pkg;
julia> Pkg.activate(".");
julia> Pkg.instantiate();
```

We can then include the sensitivity analysis code, which also loads the `Earth4All` module, by running
```jl
julia> include("src/E4ASA.jl");
```

### Reproducing the TLTL scenario and the GL scenario

The two scenarios can be reproduced by running
```jl
julia> e4a, _, prob = e4a_sys_prob();
julia> tltl = compute_tltl_sol(prob);
julia> gl = compute_gl_sol(prob);
```
The final values of the six indicators in the two scenarios can be retrieved by running
```jl
julia> tltl[e4a.AWBI][end]
julia> gl[e4a.AWBI][end]
julia> tltl[e4a.GDPP][end]
julia> gl[e4a.GDPP][end]
julia> tltl[e4a.INEQ][end]
julia> gl[e4a.INEQ][end]
julia> tltl[e4a.OW][end]
julia> gl[e4a.OW][end]
julia> tltl[e4a.POP][end]
julia> gl[e4a.POP][end]
julia> tltl[e4a.STE][end]
julia> gl[e4a.STE][end]
```
These values are the ones reported in the table of Fig. 1 of the paper. The figure itself can be reproduced by running

```jl
julia> plot_two_sol_var(e4a, "TLTL", tltl, "GL", gl, collect(1:6), "", false, true)
```
(actually the result is an interactive version of Fig. 1). One can also plot one scenario only by running, for example,
```jl
julia> plot_sol_var(e4a, tltl, collect(1:6), "The TLTL scenario", false, true)
```

### Reproducing the GL18 scenario

This scenario can be reproduced by running
```jl
julia> gl18 = compute_gl18_sol(prob);
```
The MRE shown in Table 1 can be computed by running

```jl
julia> compare(gl[e4a.AWBI], gl18[e4a.AWBI],1)
julia> compare(gl[e4a.GDPP], gl18[e4a.GDPP],1)
julia> compare(gl[e4a.INEQ], gl18[e4a.INEQ],1)
julia> compare(gl[e4a.OW], gl18[e4a.OW],1)
julia> compare(gl[e4a.POP], gl18[e4a.POP],1)
julia> compare(gl[e4a.STE], gl18[e4a.STE],1)
```

