#[=======================================================================]
#FindPrecice
#-------

#Finds the Precice library.

#Imported Targets
#^^^^^^^^^^^^^^^^

#This module provides the following imported targets, if found:

#``Precice::Precice``
#The Precice library

#Result Variables
#^^^^^^^^^^^^^^^^

#This will define the following variables:

#``Precice_FOUND``
#True if the system has the Precice library.
#``Precice_VERSION``
#The version of the Precice library which was found.
#``Precice_INCLUDE_DIRS``
#Include directories needed to use Precice.
#``Precice_LIBRARIES``
#Libraries needed to link to Precice.

#Cache Variables
#^^^^^^^^^^^^^^^

#The following cache variables may also be set:

#``Precice_INCLUDE_DIR``
#The directory containing ``Precice.h``.
#``Precice_LIBRARY``
#The path to the Precice library.

#]=======================================================================]

# Try to grab information from PkgConfig
find_package(PkgConfig)
pkg_check_modules(PC_Precice QUIET Precice)

if(NOT Precice_DIR)
  find_path(Precice_INCLUDE_DIR
    NAMES precice/  PATHS ${PC_Precice_INCLUDE_DIRS}
    PATH_SUFFIXES precice
    )
  find_library(Precice_LIBRARY
    NAMES precice
    PATHS ${PC_Precice_LIBRARY_DIRS}
    )
else()
  find_path(Precice_INCLUDE_DIR
    NAMES precice/SolverInterface.hpp
    PATHS ${Precice_DIR}/include
    PATH_SUFFIXES precice
    NO_DEFAULT_PATH
    )
  find_library(Precice_LIBRARY
    NAMES precice
    PATHS  ${Precice_DIR}/lib
    NO_DEFAULT_PATH
    )
endif()

set(Precice_VERSION ${PC_Precice_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Precice
  FOUND_VAR Precice_FOUND
  REQUIRED_VARS
  Precice_LIBRARY
  Precice_INCLUDE_DIR
  VERSION_VAR Precice_VERSION
  )


if(NOT Precice_FOUND)
  message("Precice library not found. Define Precice_DIR with path to Precice.")
  return()
endif()

if(Precice_FOUND)
  set(Precice_LIBRARIES ${Precice_LIBRARY})
  set(Precice_INCLUDE_DIRS ${Precice_INCLUDE_DIR})
  set(Precice_DEFINITIONS ${PC_Precice_CFLAGS_OTHER})
endif()

if(Precice_FOUND AND NOT TARGET Precice::Precice)
  add_library(Precice::Precice UNKNOWN IMPORTED)
  set_target_properties(Precice::Precice PROPERTIES
    IMPORTED_LOCATION "${Precice_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_Precice_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${Precice_INCLUDE_DIR}"
    )
endif()

mark_as_advanced(
  Precice_INCLUDE_DIR
  Precice_LIBRARY
  )
