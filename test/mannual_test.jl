include("../src/RandLoctions.jl")
include("../src/Visualization.jl")
using SpinEnsembles: rand_bath_dipolar_coefs, visual_ensemble


## Check datatype passing
rand_bath_dipolar_coefs(1000,2,(1,3),method=:cubic)
rand_bath_dipolar_coefs(1000,3,(1,3),method=:spherical)
rand_bath_dipolar_coefs(1000,2,(1,3),method=:spherical)
rand_bath_dipolar_coefs(1000,1,(1,3),method=:spherical)

## Check by Visualization
visual_ensemble(rand_locs_spherical(1,3,N=300))
visual_ensemble(rand_locs_polar(1,10,N=300))

visual_ensemble(rand_locs_cubic(1,2,N=300,dim=3))
visual_ensemble(rand_locs_cubic(1,2,N=300,dim=2))
visual_ensemble(rand_locs_cubic(1,2,N=300,dim=1))

visual_ensemble(rand_locs(300,3,2))