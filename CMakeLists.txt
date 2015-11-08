cmake_minimum_required (VERSION 2.6)
project (STFSolvers)

#TODO: switch between Opt and Debug mode

#The build is based on OpenFOAM environmental varialbles Need to check 
#if they are loaded and if not TODO search and load an appropriate file?
if(DEFINED ENV{WM_PROJECT})
  message("\nOpenFOAM $ENV{WM_PROJECT_VERSION} found in $ENV{WM_PROJECT_DIR}.\n")
else()
  message("WM_PROJECT is not defined. You must set OpenFOAM specific environment
  variables for CMake to configure the build.")
  # exit early 
  return()
endif()

#TODO Is there a better way to break lines?
set(CMAKE_CXX_FLAGS "-m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWM_LABEL_SIZE=32 -std=c++11")

if (CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wold-style-cast -Wnon-virtual-dtor")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-variable")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNoRepository")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-depth-100")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Xlinker --add-needed")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xlinker --no-as-needed")
else ()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGTEST_HAS_TR1_TUPLE=0")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd327,654,819,1125,1476,1505,1572")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xSSE3 -O2 -no-prec-div")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNoRepository")
endif()

#define a macro to simplify adding OpenFOAM libraries
#TODO  add for each
macro (TARGET_LINK_OPENFOAM_LIBRARIES target)
  foreach(arg ${ARGN})
    if(NOT DEFINED ${arg}_LIB)
      find_library(${arg}_LIB lib${arg}.so $ENV{FOAM_LIBBIN})
    endif()
    target_link_libraries(${target} ${${arg}_LIB})
  endforeach(arg ${ARGN})
endmacro (TARGET_LINK_OPENFOAM_LIBRARIES)

#TODO  add for each
macro (TARGET_LINK_OPENFOAM_USERLIB arg1 arg2)
  if(NOT DEFINED ${arg2}_LIB)
    find_library(${arg2}_LIB lib${arg2}.so $ENV{FOAM_USER_LIBBIN})
  endif()
  #message("Using ${${arg2}_LIB}")
  target_link_libraries(${arg1} ${${arg2}_LIB})
endmacro (TARGET_LINK_OPENFOAM_USERLIB)

macro (ADD_MY_TEST testname executable)
  add_test(
      NAME ${testname}
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${testname}"
      COMMAND ${executable}
  )

  add_custom_command(
      TARGET ${executable}
      PRE_BUILD
      COMMAND mkdir -p ${testname}
      COMMAND cp -r ${PROJECT_SOURCE_DIR}/testing/${testname}/* ${testname}
      COMMENT "Creating and copying directories for ${testname}"
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

  add_custom_command(
      TARGET ${executable}
      POST_BUILD
      COMMAND blockMesh &> log.blockMesh
      COMMENT "Running blockMesh in ${CMAKE_CURRENT_BINARY_DIR}/${testname}"
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${testname}
  )
endmacro (ADD_MY_TEST)

macro (ALLRUN_TEST testname executable)
  add_test(
    NAME ${testname}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${testname}"
    COMMAND Allrun
  )
  set_tests_properties (${testname}
    PROPERTIES PASS_REGULAR_EXPRESSION "Pass"
  )

  add_custom_command(
    TARGET ${executable}
    PRE_BUILD
    COMMAND mkdir -p ${testname}
    COMMAND cp -r ${PROJECT_SOURCE_DIR}/testing/${testname}/* ${testname}
    COMMENT "Creating and copying directories for ${testname}"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endmacro (ALLRUN_TEST)

enable_testing()

option (UseOldLibs "Use old libraries" OFF)
#General includes required by all code
include_directories ("$ENV{FOAM_SRC}/fvOptions/lnInclude")
include_directories ("$ENV{FOAM_SRC}/finiteVolume/lnInclude")
include_directories ("$ENV{FOAM_SRC}/meshTools/lnInclude")
include_directories ("$ENV{FOAM_SRC}/OpenFOAM/lnInclude")
include_directories ("$ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude")

#include_directories(./gtest/include)

#We need to tell CMake to use use C files as CXX (OpenFOAM standard)
set_source_files_properties(filename.c PROPERTIES LANGUAGE CXX )

# gtest for testing
# add_subdirectory(ext/gtest/)
# include_directories(${GTEST_INCLUDE_DIRS})


# Library
add_subdirectory(src)
