ensemble=SpinEnsemble(0.01,3,[0,0,1],0.1,30.0,:spherical)

spins=SpinCluster(ensemble)


rabi(spins,100)
t=0:0.001:2
averagerabi(t, ensemble, 10; M=500)