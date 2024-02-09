###############################################################################
# cmake file for building the Higgs-BR-classes analysis package.
# @author Remi Ete, DESY (original)
# @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP (adaptation)
# @author Fabricio Jimenez Morales, LLR, CNRS, École Polytechnique, IPP (adaptation)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
###############################################################################

# Project name.
PROJECT( Digitization )

# Project version.
SET( ${PROJECT_NAME}_VERSION_MAJOR 1 )
SET( ${PROJECT_NAME}_VERSION_MINOR 0 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

### DEPENDENCIES ##############################################################
FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

# Load default settings from ILCSOFT_CMAKE_MODULES.
INCLUDE( ilcsoft_default_settings )

FIND_PACKAGE( Marlin 1.0 REQUIRED ) # Minimum required Marlin version.
INCLUDE_DIRECTORIES( SYSTEM ${Marlin_INCLUDE_DIRS} )
LINK_LIBRARIES( ${Marlin_LIBRARIES} )
ADD_DEFINITIONS( ${Marlin_DEFINITIONS} )

FIND_PACKAGE( MarlinUtil 1.0 REQUIRED ) # Minimum required MarlinUtil version.
INCLUDE_DIRECTORIES( SYSTEM ${MarlinUtil_INCLUDE_DIRS} )
LINK_LIBRARIES( ${MarlinUtil_LIBRARIES} )
ADD_DEFINITIONS( ${MarlinUtil_DEFINITIONS} )

FIND_PACKAGE( ROOT REQUIRED COMPONENTS GenVector)
INCLUDE_DIRECTORIES( SYSTEM ${ROOT_INCLUDE_DIRS} )
LINK_LIBRARIES( ${ROOT_LIBRARIES} )
ADD_DEFINITIONS( ${ROOT_DEFINITIONS} )

### LIBRARY ###################################################################
# Register all processors/ packages with cmake.
MACRO( ADD_PROCESSOR _input_dir )
    MESSAGE( STATUS "Including your processor: ${_input_dir}" )
    INCLUDE_DIRECTORIES( BEFORE ${_input_dir}/include )
    AUX_SOURCE_DIRECTORY( ${_input_dir}/src project_cxx_srcs )
ENDMACRO()

# -----------------------------------------------------------------------------
# Add the packages you want to compile here.
ADD_PROCESSOR( ./LCIO2BuildProcessor )
ADD_PROCESSOR( ./DigiLCIO2BuildProcessor )
ADD_PROCESSOR( ./ConversionProcessor )
ADD_PROCESSOR( ./ShapingProcessor    )
#ADD_PROCESSOR( ./ShapedConversionProcessor )

# Build the project library.
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${project_cxx_srcs} )
INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION lib )

# Display some variables and write them to cache.
DISPLAY_STD_VARIABLES()
