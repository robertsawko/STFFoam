file(GLOB lib_source_files
  nonInertialFrame.C
)

add_library(STF SHARED ${lib_source_files})
