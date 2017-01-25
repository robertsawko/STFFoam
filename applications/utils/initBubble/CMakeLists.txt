add_executable(initBubble initBubble.C)

#Includes specfic to the project

# include_directories ("$ENV{FOAM_APP}/solvers/multiphase/interFoam")
include_directories ("$ENV{FOAM_SRC}/transportModels/twoPhaseMixture/lnInclude")
include_directories ("$ENV{FOAM_SRC}/transportModels/incompressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/transportModels/interfaceProperties/lnInclude")
include_directories ("$ENV{FOAM_SRC}/transportModels/turbulenceModels/lnInclude")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/incompressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/turbulenceModels/lnInclude")
include_directories ("$ENV{FOAM_SRC}/transportModels/immiscibleIncompressibleTwoPhaseMixture/lnInclude")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/turbulenceModels/lnInclude")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/incompressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/transportModels/incompressible/singlePhaseTransportModel")
include_directories ("$ENV{FOAM_SRC}/transportModels/")
include_directories ("$ENV{FOAM_SRC}/sampling/lnInclude")

#TODO Is there a better way to add multiple OF internal libraries?
target_link_openfoam_libraries(initBubble
  immiscibleIncompressibleTwoPhaseMixture turbulenceModels
  incompressibleTurbulenceModels finiteVolume meshTools
  fvOptions sampling OpenFOAM)

target_link_libraries(initBubble STF)

set_target_properties(initBubble PROPERTIES RUNTIME_OUTPUT_DIRECTORY "$ENV{FOAM_USER_APPBIN}")

