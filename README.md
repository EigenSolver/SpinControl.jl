# SpinControl.jl

**SpinControl.jl** is a numerical package written in Julia that makes it easy to simulate the dynamics and decoherence of a spin qubit interacting with environmental spin ensembles.

The package provides necessary data types for environmental spin bath, driving pulses, and quantum control protocols, as well as a collection of functions that are widely used in magnetic resonance and quantum control, such as free induction decay (FID) and Rabi driving. 
Users can easily use these abstractions to compose their own control sequences such as dynamical decoupling protocols and test their performance on various spin ensembles. 

**The package is archived as the support material of my master thesis.**

**[Quantum Control of Interacting Spins, Y.N.Zhang, Technische Universiteit Delft](http://resolver.tudelft.nl/uuid:9673e6ba-ff3e-402d-938c-f8d10508cddf)**

https://doi.org/10.4121/19766887

## Spin Ensemble

- Spin 1/2, randomly distributed, positional disorder.
- Dipole-dipole interactions.
- Tunable parameters like filling rate, cluster shape, cluster size, Zeeman field direction.  
- Driving via external magnetic field.
- Lattice structure. (To be added).

## Pulse Sequence

- Arbitrary driving axis, strength and duration.
- Various shape of driving pulse via composition.
- Customized pulse sequence with free combination of pulse shapes and intervals.
- Common dynamical decoupling protocols integrated.

## Spin Dynamics

- Efficient simulation via Monte-Carlo sampling.
- Highly optimized for common experiments like free induction decay and Rabi driving.
- Ideal evolution unitary for pulse sequence, composed from ideal rotaion unitaries of pulses.
- Monte-Carlo sampled Kraus operators for quantum operation on open subsystem.
- Realistic dynamical decoupled system dynamics via serial quantum operations. (To be added)

## Quantum Information

- Effective characterization of gate and operation fidelity.
- Realistic gate design.

## Development Schedule

1. Fast prototype and draft documentation v0.1.0
2. `SpinEnsemble` type and method, v0.2.0
3. `PulseSequence` type and methods v0.3.0
4. Noise Simulation... v0.4.0
