using IfElse
using Formatting
using Graphs
using Latexify
using PlotlyJS
using Printf

using WorldDynamics

function add_equation!(eqs, equation)
   append!(eqs, [equation])
end

function delay_n!(eqs, x, rt, lv, delay, order)
   append!(eqs, [rt[1] ~ lv[1] / (delay / order)])
   append!(eqs, [D(lv[1]) ~ x - rt[1]])
   for d in 2:order
      append!(eqs, [rt[d] ~ lv[d] / (delay / order)])
      append!(eqs, [D(lv[d]) ~ rt[d-1] - rt[d]])
   end
end

ramp(x, slope, startx, endx) = IfElse.ifelse(x > startx, IfElse.ifelse(x < endx, slope * (x - startx), slope * (endx - startx)), 0)

function smooth!(eqs, x, input, delay_time)
   append!(eqs, [D(x) ~ (input - x) / delay_time])
end

clip(returnifgte, returniflt, inputvalue, threshold) = IfElse.ifelse(inputvalue ≥ threshold, returnifgte, returniflt)

step(inputvalue, returnifgte, threshold) = clip(returnifgte, zero(returnifgte), inputvalue, threshold)


pulse(inputvalue, start, width) = IfElse.ifelse(inputvalue >= start, 1, 0) * IfElse.ifelse(inputvalue < (start + width), 1, 0)

interpolate(x, x₁, xₙ, y₁, yₙ) = y₁ + (x - x₁) * ((yₙ - y₁) / (xₙ - x₁))

function interpolate(x, yvalues::Tuple{Vararg{Float64}}, xrange::Tuple{Float64,Float64})
   interpolate(x, collect(yvalues), collect(LinRange(xrange[1], xrange[2], length(yvalues))))
end

function interpolate(x, yvalues::Vector{Float64}, xvalues::Vector{Float64})
   # y gets the min y value if less than the min x range
   #   or the max y value if greater than the max x range
   y = (x < xvalues[1]) * yvalues[1] + (x ≥ xvalues[end]) * yvalues[end]

   # in case x is inside the range, y gets the interpolated value
   for i ∈ 1:length(yvalues)-1
      y += (x ≥ xvalues[i]) * (x < xvalues[i+1]) * interpolate(x, xvalues[i], xvalues[i+1], yvalues[i], yvalues[i+1])
   end

   return y
end

function withlookup(x, pairs::Vector{Tuple{Float64,Float64}})
   interpolate(x, map(t -> t[end], pairs), map(t -> t[1], pairs))
end