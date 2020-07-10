These scripts accompany the paper:

"Integrating dynamic plant growth models and microclimates for species
distribution modelling" by Rafael Schouten, Peter Vesk and Michael Kearney, in the
Journal of Ecological modelling.

![Transect simulaitons over a decade](https://media.githubusercontent.com/media/rafaqz/DEBScripts/master/plots/all.png)

Scripts are provided for building al the plots from the paper, and also for
using the interactive user-interface used to examine and simplify the model.


To set up, run julia in this folder, then run:

```julia-repl
]
activate .
instantiate
```

Where `]` gets you into Pkg mode in the REPL. 

After that you can use the scripts:

- ui.jl loads the interactive user interface
- paper_plots.jl builds plots for the paper

As loading all the packages will take quite a while, load these in a julia
session that you keep open (ie in Atom/Juno), instead of running them from the
command line each time.


If you have any problems, open an issue in this repository, including your
version of Julia, description of the problem and error outputs where necessary.
