
using InteractNext
using AutoInteract
using Blink
using CSSUtil
import Plots.px
using Plots
using UnitfulPlots
using AxisArrays

plot(sol)
fields = fieldnames(organism.nodes[1].state)
labels = reshape(vcat(string.(fields, "1"), string.(fields, "2")), 1, 12)
p = plot(sol.t, sol', size=(1000, 700))

gr()
r1 = organism.records[1]; nothing
r2 = organism.records[2]; nothing
p = plot([r2.vars[i].assimilation.X_H for i in 1:8000], size=(1000, 700))
p = plot([r1.vars[i].rate for i in 1:8000], size=(1000, 700))
p = plot!(p,[r2.vars[i].rate for i in 1:8000], size=(1000, 700))
p = plot([r2.vars[i].scale for i in 1:8000], size=(1000, 700))
p = plot([r1.vars[i].assimilation.gs for i in 1:8000], size=(1000, 700))
p = plot([r1.vars[i].assimilation.ac for i in 1:8000], size=(1000, 700))
p = plot([r1.vars[i].assimilation.aj for i in 1:8000], size=(1000, 700))
p = plot([r1.vars[i].assimilation.par for i in 1:8000], size=(1000, 700))
p = plot([r1.vars[i].assimilation.km for i in 1:8000], size=(1000, 700))
p = plot([r1.vars[i].assimilation.rh for i in 1:8000], size=(1000, 700))
p = plot([r1.vars[i].assimilation.rd for i in 1:8000], size=(1000, 700))
p = plot([r1.vars[i].assimilation.vj for i in 1:8000], size=(1000, 700))
p = plot([r1.vars[i].assimilation.vcmax for i in 1:8000], size=(1000, 700))
p = plot([r1.vars[i].assimilation.jmax for i in 1:8000], size=(1000, 700))
p = plot([r1.vars[i].assimilation.aleaf for i in 1:8000], size=(1000, 700))
p = plot([r1.vars[i].assimilation.tleaf for i in 1:8000], size=(1000, 700))
p = plot!(p,[r1.vars[i].assimilation.tair for i in 1:8000], size=(1000, 700))
p = plot([r1.vars[i].assimilation.ca for i in 1:8000], size=(1000, 700))
p = plot([r1.vars[i].assimilation.cs for i in 1:8000], size=(1000, 700))
p = plot([r1.vars[i].assimilation.ci for i in 1:8000], size=(1000, 700))
p = plot([r1.vars[i].assimilation.gh for i in 1:8000], size=(1000, 700))
p = plot([r1.vars[i].assimilation.et for i in 1:8000], size=(1000, 700))
p = plot([ustrip.(r1.J[i][:,4]) for i in 1:8000], size=(1000, 700))
p = plot([r1.vars[i].height for i in 1:8000], size=(1000, 700))
p = plot!(p,[r2.vars[i].height for i in 1:8000], size=(1000, 700))


import AutoInteract: plotit!, widgetable
function plotit!(p, parent::AxisArray, data, label)
    range = inbounds_range(data)
    plot!(p, data[checked_range], label=label, layout=AutoInteract.deflayout)
end

widgetable(T, fieldname) = false#flattenable(T, fieldname)

wi = make_widgets(organism)
pl = make_plottables(organism)
# delete!(w, :vars)
pl = make_plottables(organism)
sp = get_signals(pl)
sw = get_signals(wi)
pl = make_plottables(organism)
sp = get_signals(pl)
ui = make_interface((pl, wi), box=vbox)
# settings = apply_all(settings, sw.value)
# plot_all(settings, sp.value, 1:10)
w = Window()
body!(w, ui)

map((sw, sp) -> debplot(p, scenario, u0, sw, sp), get_signals(wi), get_signals(pl))
