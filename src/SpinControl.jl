module SpinControl

using Reexport
# include("operators.jl")

include("sequences/PulseSequences.jl")
include("ensembles/SpinEnsembles.jl")
@reexport using .SpinEnsembles, .PulseSequences

include("operations.jl")
export evolution, operation, rotation, channel, krausoperators, isunitary
export statefidelity, processfidelity

end
