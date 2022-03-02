# SpinControl.jl

**SpinControl.jl** is a numerical package written in Julia that makes it easy to simulate the decoherence process and analyze the dynamical decoupling protocols in interacting spin ensembles.

The package provides necessary data types for spin ensembles and control pulses, as well as a collection of functions that are widely used in magnetic resonance and spin dynamics, such as free induction decay (FID) and Rabi driving.

Users can easily use these abstractions to build their own control pulse sequences such as dynamical decoupling protocols and test their performance on various spin ensembles.

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

- Effective characterization of gate and operation fidelity, via sampling.
- Realistic gate design.

## Development Schedule

**The package is currently under developing.**

1. Fast prototype and draft documentation v0.1.0
2. `SpinEnsemble` type and method, v0.2.0
3. `PulseSequence` type and methods v0.3.0
4. Noise Simulation... v0.4.0
