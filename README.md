These scripts accompany the paper:

"Integrating dynamic plant growth models and microclimates for species
distribution modelling" by Rafael Schouten, Peter Vesk and Michael Kearney, in the
Journal of Ecological modelling ([preprint](https://ecoevorxiv.org/ja4m6)).

![Transect simulations](https://media.githubusercontent.com/media/rafaqz/DEBplant/master/plots/transect_multiplot.png)

Scripts are provided for building al the plots from the paper, and also for
using the interactive user-interface used to examine and simplify the model.


To set up, run julia in this folder, then run:

```julia-repl
]
activate .
instantiate
```

Where `]` gets you into Pkg mode in the REPL. 

After that you can hit escape to leave the Pkg mode, and use the scripts:

- ui.jl loads the interactive user interface
- paper_plots.jl builds figures for the paper
- paper_maps.jl builds mapped figures for the paper. By default this will run usig pre-save
  simulaitons, but requires very large download (380 GB unzipped) and long build time
  to run from scratch.

As loading all the packages will take quite a while, load these in a julia
session that you keep open (ie in Atom/Juno), instead of running them from the
command line each time.

If you have any problems, open an issue in this repository, including your
version of Julia, description of the problem and error outputs where necessary.


![User interface](https://media.githubusercontent.com/media/rafaqz/DEBplant/master/plots/ui.png)
_User interface generated by the `ui.jl` script_
