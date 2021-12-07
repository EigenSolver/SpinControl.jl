include("../src/RandLoctions.jl")
include("../src/Visualization.jl")
using SpinEnsembles: randcoefs


## Check datatype passing
randcoefs(1000,2,(1,3),method=:cubic)
randcoefs(1000,3,(1,3),method=:spherical)
randcoefs(1000,2,(1,3),method=:spherical)
randcoefs(1000,1,(1,3),method=:spherical)

## Check by Visualization
visualensemble(randlocsspherical(1,3,N=300))
visualensemble(randlocspolar(1,3,N=300))

visualensemble(randlocscubic(1,2,N=300,dim=3))
visualensemble(randlocscubic(1,2,N=300,dim=2))
visualensemble(randlocscubic(1,2,N=300,dim=1))

visualensemble(randlocs(300,3,2))

M=randlocsspherical(1,10,N=1000)
minimum(x->abs(norm(x)), eachrow(M))>1