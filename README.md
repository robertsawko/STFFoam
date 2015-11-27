# Translational frames for OpenFOAM

An equivalent of MFR for translational frames for OpenFOAM solvers.

 * At the moment instrumentation of solvers is necessary

## Compilation

Using GNU g++ should work fine by 
```bash
# Load OpenFOAM environement
cd <source>
mkdir build #for out of source build
cd build
cmake ../
make #or make -j 4 for a four thread build
```

On Astral we need to use Intel compilers and old libraries. I created a special
CMake variable which will switch between old and new compilation rules and files
for different stuff.

```bash
cd <source>
mkdir build #for out of source build
cd build
cmake -DCMAKE_CXX_COMPILER=icpc ../
make
```

## Solver instrumentation

Create the translational frame e.g.

```cpp
translationalFrame frameOfReference;

// This is in order to force the BC change on fixed boundaries
frameOfReference.registerVelocity(U); 
```

Update the translational frame somewhere in the time or steady state convergence
loop. Pressure field and the deviator of the stress tensor have to be passed
into the function.

```cpp
frameOfReference.update(p, turbulence->devReff());
```

