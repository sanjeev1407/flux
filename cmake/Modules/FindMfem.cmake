#FindMfem
#-------
#
#Finds the MFEM library.
#
#Imported Targets
#^^^^^^^^^^^^^^^^
#
#This module provides the following imported targets, if found:
#
#``Mfem::Mfem``
#The MFEM library
#
#Result Variables
#^^^^^^^^^^^^^^^^
#
#This will define the following variables:
#
#``Mfem_FOUND``
#True if the system has the Mfem library.
#``Mfem_VERSION``
#The version of the Mfem library which was found.
#``Mfem_INCLUDE_DIRS``
#Include directories needed to use Mfem.
#``Mfem_LIBRARIES``
#Libraries needed to link to Mfem.
#
#Cache Variables
#^^^^^^^^^^^^^^^
#
#The following cache variables may also be set:
#
#``Mfem_INCLUDE_DIR``
#The directory containing ``mfem.h``.
#``Mfem_LIBRARY``
#The path to the Mfem library.
#
#]=======================================================================]

# Try to grab information from PkgConfig
find_package(PkgConfig)
pkg_check_modules(PC_Mfem QUIET Mfem)

if(NOT Mfem_DIR)
  find_path(Mfem_INCLUDE_DIR
    NAMES mfem.hpp
    PATHS ${PC_Mfem_INCLUDE_DIRS}
    PATH_SUFFIXES mfem
    )
  find_library(Mfem_LIBRARY
    NAMES mfem
    PATHS ${PC_Mfem_LIBRARY_DIRS}
    )
else()
  find_path(Mfem_INCLUDE_DIR
    NAMES mfem.hpp
    PATHS ${Mfem_DIR}/include
    NO_DEFAULT_PATH
    )
  find_library(Mfem_LIBRARY
    NAMES mfem
    PATHS ${Mfem_DIR}/lib
    NO_DEFAULT_PATH
    )
endif()

set(Mfem_VERSION ${PC_Mfem_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Mfem
  FOUND_VAR Mfem_FOUND
  REQUIRED_VARS
  Mfem_LIBRARY
  Mfem_INCLUDE_DIR
  VERSION_VAR Mfem_VERSION
  )

if(NOT Mfem_FOUND)
  message("MFEM library not found. Define Mfem_DIR with path to MFEM.")
  return()
endif()

if(Mfem_FOUND)
  set(Mfem_LIBRARIES ${Mfem_LIBRARY})
  set(Mfem_INCLUDE_DIRS ${Mfem_INCLUDE_DIR})
  set(Mfem_DEFINITIONS ${PC_Mfem_CFLAGS_OTHER})
endif()

if(Mfem_FOUND AND NOT TARGET Mfem::Mfem)
  add_library(Mfem::Mfem UNKNOWN IMPORTED)
  set_target_properties(Mfem::Mfem PROPERTIES
    IMPORTED_LOCATION "${Mfem_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_Mfem_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${Mfem_INCLUDE_DIR}"
    )
endif()

mark_as_advanced(
  Mfem_INCLUDE_DIR
  Mfem_LIBRARY
  )
