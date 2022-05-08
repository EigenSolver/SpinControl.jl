module SpinControl

using Reexport

using LinearAlgebra

include("sequences/PulseSequences.jl")
include("ensembles/SpinEnsembles.jl")
@reexport using .SpinEnsembles, .PulseSequences

include("operations.jl")
include("fidelities.jl")
include("controldynamics.jl")

export operate, unitary, rotation, krausoperators, isunitary, measure
export statefidelity, processfidelity, entanglementfidelity, paulifidelity
export carrfidelity, xyfidelity
export deploy 

end
