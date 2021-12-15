ensemble=SpinEnsemble(0.01,3,[0,0,1],0.1,30.0,:spherical)

spins=SpinCluster(ensemble)


rabi(spins,100)
averagerabi(ensemble)