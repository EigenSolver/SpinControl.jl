module SpinControl

using Reexport
# include("operators.jl")

include("sequences/PulseSequences.jl")
include("ensembles/SpinEnsembles.jl")
@reexport using .SpinEnsembles, .PulseSequences

include("operators.jl")
export evolution, operation, rotation, channel, fidelity, isunitary
export statefidelity, processfidelity

end
