file(GLOB files
  $ENV{FOAM_APP}/solvers/multiphase/driftFluxFoam/incompressibleTwoPhaseInteractingMixture/incompressibleTwoPhaseInteractingMixture.C
  $ENV{FOAM_APP}/solvers/multiphase/driftFluxFoam/compressibleTurbulenceModels.C
  driftFluxFoam.C
)
add_executable(STFDriftFluxFoam ${files})

#Includes specfic to the project
include_directories ("$ENV{FOAM_APP}/solvers/multiphase/driftFluxFoam")
include_directories ("$ENV{FOAM_APP}/solvers/multiphase/driftFluxFoam/mixtureViscosityModels/lnInclude")
include_directories ("$ENV{FOAM_APP}/solvers/multiphase/driftFluxFoam/incompressibleTwoPhaseInteractingMixture")
include_directories ("$ENV{FOAM_APP}/solvers/multiphase/driftFluxFoam/relativeVelocityModels/lnInclude")
include_directories ("$ENV{FOAM_SRC}/transportModels/compressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/transportModels/twoPhaseMixture/lnInclude")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/turbulenceModels/lnInclude")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/compressible/lnInclude")

#TODO Is there a better way to add multiple OF internal libraries?
target_link_openfoam_libraries(STFDriftFluxFoam
  driftFluxTransportModels driftFluxRelativeVelocityModels
  compressibleTransportModels compressibleTurbulenceModels finiteVolume
  meshTools fvOptions sampling OpenFOAM)

target_link_libraries(STFDriftFluxFoam STF)

set_target_properties(STFDriftFluxFoam PROPERTIES RUNTIME_OUTPUT_DIRECTORY "$ENV{FOAM_USER_APPBIN}")

