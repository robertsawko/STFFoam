file(GLOB lib_source_files
  translationalFrame.C
  constantTranslation.C
  fallingObject.C
  centreOfVolume.C
  bypassPig.C
  frameAwareBoundary.C
  translatingFixedVelocityFvPatchField.C
  STFScouringInflowVelocityFvPatchVectorField.C
)

add_library(STF SHARED ${lib_source_files})
set_target_properties (STF PROPERTIES LIBRARY_OUTPUT_DIRECTORY 
    "$ENV{FOAM_USER_LIBBIN}")
