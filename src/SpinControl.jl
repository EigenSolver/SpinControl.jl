module SpinControl

using Reexport

using LinearAlgebra

include("sequences/PulseSequences.jl")
include("ensembles/SpinEnsembles.jl")
@reexport using .SpinEnsembles, .PulseSequences

include("operations.jl")
include("fidelities.jl")
export evolution, operation, rotation, channel, krausoperators, isunitary
export statefidelity, processfidelity, paulifidelity

end
