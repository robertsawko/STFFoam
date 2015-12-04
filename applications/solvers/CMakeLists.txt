include_directories ("${PROJECT_SOURCE_DIR}/src/")
include_directories ("$ENV{FOAM_SRC}/transportModels/")
include_directories ("$ENV{FOAM_SRC}/sampling/lnInclude")

add_subdirectory(STFPimpleFoam)
add_subdirectory(STFInterFoam)
