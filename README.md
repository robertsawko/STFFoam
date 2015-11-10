# Single Translational Frames for OpenFOAM

An equivalent of MFR for translational frames for OpenFOAM solvers.

 * At the moment instrumentation of solvers is necessary

## Solver instrumentation


```cpp
translationalFrame frameOfReference;

// This is in order to force the BC change on fixed boundaries
frameOfReference.registerVelocity(U); 
```
