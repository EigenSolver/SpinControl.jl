"""
A julia package that provides necessary functionalities to study the 
quantum decoherence in disordered spin ensembles.
"""
module SpinEnsembles

using LinearAlgebra
using Statistics
import SpecialFunctions: gamma
import ProgressMeter: @showprogress

# datetpye
export SpinEnsemble, SpinCluster

# mathods and functions
export fid, fid, rabi, rabi, betasampling, coherencetime, randlocs, randcoefs, dipolarlinewidth

include("randloctions.jl")
include("spindynamics.jl")
include("dipolarcouplings.jl")

end