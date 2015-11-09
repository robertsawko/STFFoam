add_executable(STFPimpleFoam pimpleFoam.C)

#Includes specfic to the project
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/turbulenceModels/lnInclude")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/incompressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/transportModels/")
include_directories ("$ENV{FOAM_SRC}/transportModels/incompressible/singlePhaseTransportModel")
include_directories ("$ENV{FOAM_SRC}/sampling/lnInclude")
include_directories ("${PROJECT_SOURCE_DIR}/src/translationalFrames/")

#TODO Is there a better way to add multiple OF internal libraries?
target_link_openfoam_libraries(STFPimpleFoam
  turbulenceModels incompressibleTurbulenceModels
  incompressibleTransportModels finiteVolume meshTools
  fvOptions sampling OpenFOAM)

#target_link_openfoam_userlib(STFInterFoam translationalFrames)

#target_link_libraries(STFPimpleFoam translationalFrames)

set_target_properties(STFPimpleFoam PROPERTIES RUNTIME_OUTPUT_DIRECTORY "$ENV{FOAM_USER_APPBIN}")

