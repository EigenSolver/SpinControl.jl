# SpinEnsembles.jl

**SpinEnsembles.jl** is a numerical package written in Julia that makes it easy to simulate the decoherence process and analyze the dynamical decoupling protocols in interacting spin ensembles.

The package implement functions that is widely used in magnetic resonance and spin dynamics, 
such as free induction decay (FID) and Rabi driving. Users can easily use these functions to build their own control pulse sequence and test it on various spin ensembles.

## Spin Ensemble
- Spin 1/2, randomly distributed, positional disorder.
- Dipole-dipole interactions.
- Lattice structure. (to be added).
- Driving via external magnetic field.

## Pulse Sequence
- Periodic external field.
- Tunable pulse duration and periodic pattern.
- Time modulated on different directions. 
- Custom noises. 

The package is currently under developing. 

## Development Schedule 
1. Fast prototype and draft documentation v0.1.0
2. `SpinEnsemble` type and method, v0.5.0
3. `PulseSequence` type and methods v0.9.0
4. Noise Simulation... v1.0.0
5. 

*Todo list*
- [ ] Add control noise for `PulseSequence`
- [ ] Add methods interaction of `PulseSequence` and `SpinEnsemble`
- [ ] Implement `PulseSequence` datatype
- [ ] Add unitests for dilute ensemble and linewidth
- [ ] Implement methods based on `SpinEnsemble`
- [ ] Add lattice structure for `SpinEnsemble`.
- [ ] Implement `SpinEnsemble` datatype 
- [ ] Implement ensemble averaged Rabi oscillation
- [ ] Add continous integration (CI)
- [ ] Deploy the documention on github IO
- [ ] Use `Union` to implement multi-type functions
- [ ] Use `@Macro` for easier visualization
- [ ] Make unified modeling language (UML) diagram, for the architecture
- [X] Reorganize the location functions
- [X] Optimize APIs
- [X] Write better documentation
- [X] Implement Rabi oscillation
- [X] Rename functions
- [X] Write documentation
- [X] Precompile the package and use it in the notebooks
- [X] Update dependencies
- [X] Encapsulate the package
- [X] Implement visualization functions
- [X] Add unittests
- [X] Add adaptive dephasingtime of FID
- [X] Implement FID
- [X] Implement random points generation
- [X] Git version control
- [ ] Read more on design pattern and software engineering