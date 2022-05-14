#FindMetis
#-------
#
#Finds the Metis library.
#
#Imported Targets
#^^^^^^^^^^^^^^^^
#
#This module provides the following imported targets, if found:
#
#``Metis::Metis``
#The Metis library
#
#Result Variables
#^^^^^^^^^^^^^^^^
#
#This will define the following variables:
#
#``Metis_FOUND``
#True if the system has the Metis library.
#``Metis_VERSION``
#The version of the Metis library which was found.
#``Metis_INCLUDE_DIRS``
#Include directories needed to use Metis.
#``Metis_LIBRARIES``
#Libraries needed to link to Metis.
#
#Cache Variables
#^^^^^^^^^^^^^^^
#
#The following cache variables may also be set:
#
#``Metis_INCLUDE_DIR``
#The directory containing ``metis.h``.
#``Metis_LIBRARY``
#The path to the Metis library.
#
#]=======================================================================]

# Try to grab information from PkgConfig
find_package(PkgConfig)
pkg_check_modules(PC_Metis QUIET Metis)

if(NOT Metis_DIR)
  find_path(Metis_INCLUDE_DIR
    NAMES metis.h
    PATHS ${PC_Metis_INCLUDE_DIRS}
    PATH_SUFFIXES metis
    )
  find_library(Metis_LIBRARY
    NAMES metis
    PATHS ${PC_Metis_LIBRARY_DIRS}
    )
else()
  find_path(Metis_INCLUDE_DIR
    NAMES metis.h
    PATHS ${Metis_DIR}/include
    PATH_SUFFIXES metis
    NO_DEFAULT_PATH
    )
  find_library(Metis_LIBRARY
    NAMES metis
    PATHS  ${Metis_DIR}/lib
    NO_DEFAULT_PATH
    )
endif()

set(Metis_VERSION ${PC_Metis_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Metis
  FOUND_VAR Metis_FOUND
  REQUIRED_VARS
  Metis_LIBRARY
  Metis_INCLUDE_DIR
  VERSION_VAR Metis_VERSION
  )

if(NOT Metis_FOUND)
  message("Metis library not found. Define Metis_DIR with path to Metis.")
  return()
endif()

if(Metis_FOUND)
  set(Metis_LIBRARIES ${Metis_LIBRARY})
  set(Metis_INCLUDE_DIRS ${Metis_INCLUDE_DIR})
  set(Metis_DEFINITIONS ${PC_Metis_CFLAGS_OTHER})
endif()

if(Metis_FOUND AND NOT TARGET Metis::Metis)
  add_library(Metis::Metis UNKNOWN IMPORTED)
  set_target_properties(Metis::Metis PROPERTIES
    IMPORTED_LOCATION "${Metis_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_Metis_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${Metis_INCLUDE_DIR}"
    )
endif()

mark_as_advanced(
  Metis_INCLUDE_DIR
  Metis_LIBRARY
  )
