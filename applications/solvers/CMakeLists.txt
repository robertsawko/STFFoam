include_directories ("$ENV{FOAM_SRC}/transportModels/")
include_directories ("$ENV{FOAM_SRC}/sampling/lnInclude")
include_directories ("${PROJECT_SOURCE_DIR}/src/")

add_subdirectory(STFPimpleFoam)
add_subdirectory(STFInterFoam)
add_subdirectory(STFDriftFluxFoam)
