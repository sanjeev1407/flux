#FindHypre
#-------
#
#Finds the Hypre library.
#
#Imported Targets
#^^^^^^^^^^^^^^^^
#
#This module provides the following imported targets, if found:
#
#``Hypre::Hypre``
#The Hypre library
#
#Result Variables
#^^^^^^^^^^^^^^^^
#
#This will define the following variables:
#
#``Hypre_FOUND``
#True if the system has the Hypre library.
#``Hypre_VERSION``
#The version of the Hypre library which was found.
#``Hypre_INCLUDE_DIRS``
#Include directories needed to use Hypre.
#``Hypre_LIBRARIES``
#Libraries needed to link to Hypre.
#
#Cache Variables
#^^^^^^^^^^^^^^^
#
#The following cache variables may also be set:
#
#``Hypre_INCLUDE_DIR``
#The directory containing ``HYPRE.h``.
#``Hypre_LIBRARY``
#The path to the Hypre library.
#
#]=======================================================================]

# Try to grab information from PkgConfig
find_package(PkgConfig)
pkg_check_modules(PC_Hypre QUIET Hypre)

if(NOT Hypre_DIR)
  find_path(Hypre_INCLUDE_DIR
    NAMES HYPRE.h
    PATHS ${PC_Hypre_INCLUDE_DIRS}
    PATH_SUFFIXES hypre
    )
  find_library(Hypre_LIBRARY
    NAMES HYPRE
    PATHS ${PC_Hypre_LIBRARY_DIRS}
    )
else()
  find_path(Hypre_INCLUDE_DIR
    NAMES HYPRE.h
    PATHS ${Hypre_DIR}/include
    PATH_SUFFIXES hypre
    NO_DEFAULT_PATH
    )
  find_library(Hypre_LIBRARY
    NAMES HYPRE
    PATHS ${Hypre_DIR}/lib
    NO_DEFAULT_PATH
    )
endif()

set(Hypre_VERSION ${PC_Hypre_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Hypre
  FOUND_VAR Hypre_FOUND
  REQUIRED_VARS
  Hypre_LIBRARY
  Hypre_INCLUDE_DIR
  VERSION_VAR Hypre_VERSION
  )

if(NOT Hypre_FOUND)
  message("Hypre library not found. Define Hypre_DIR with path to Hypre.")
  return()
endif()

if(Hypre_FOUND)
  set(Hypre_LIBRARIES ${Hypre_LIBRARY})
  set(Hypre_INCLUDE_DIRS ${Hypre_INCLUDE_DIR})
  set(Hypre_DEFINITIONS ${PC_Hypre_CFLAGS_OTHER})
endif()

if(Hypre_FOUND AND NOT TARGET Hypre::Hypre)
  add_library(Hypre::Hypre UNKNOWN IMPORTED)
  set_target_properties(Hypre::Hypre PROPERTIES
    IMPORTED_LOCATION "${Hypre_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_Hypre_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${Hypre_INCLUDE_DIR}"
    )
endif()

mark_as_advanced(
  Hypre_INCLUDE_DIR
  Hypre_LIBRARY
  )
