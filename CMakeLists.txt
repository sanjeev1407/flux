cmake_minimum_required(VERSION 3.9)

project(flux2D VERSION 1.0
  DESCRIPTION "Finite-element soLver for Unsteady electromagnetiX"
  LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules")

### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
  message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory.")
endif()


find_package(MPI REQUIRED)
find_package(Mfem REQUIRED)
find_package(Lapack REQUIRED)
find_package(Metis REQUIRED)
find_package(Hypre REQUIRED)
find_package(Precice REQUIRED)

if(USE_ADIOS2)
find_package(ADIOS2 REQUIRED)
endif()

# Packages held in this repo in ./external
set(DEBUG_ASSERT_SINGLE_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/external/debug_assert")

add_executable(flux2D flux/flux2D.hpp flux/flux2D.cpp flux/flux2D_solver.hpp flux/flux2D_solver.cpp flux/precice_adapter.hpp flux/precice_adapter.cpp flux/mfem_coefficients.hpp)
add_executable(dummy2D flux/precice_adapter.hpp flux/precice_adapter.cpp flux/dummy2D.cpp)

target_include_directories(flux2D PUBLIC "./")
target_include_directories(flux2D PUBLIC "${DEBUG_ASSERT_SINGLE_HEADER}")
target_link_libraries(flux2D PUBLIC MPI::MPI_CXX)
target_link_libraries(flux2D PUBLIC Mfem::Mfem)
target_link_libraries(flux2D PUBLIC Metis::Metis)
target_link_libraries(flux2D PUBLIC Hypre::Hypre)
target_link_libraries(flux2D PUBLIC Precice::Precice)
target_link_libraries(flux2D PUBLIC ${LAPACK_LIBRARIES})


target_include_directories(dummy2D PUBLIC "./")
target_include_directories(dummy2D PUBLIC "${DEBUG_ASSERT_SINGLE_HEADER}")
target_link_libraries(dummy2D PUBLIC MPI::MPI_CXX)
target_link_libraries(dummy2D PUBLIC Mfem::Mfem)
target_link_libraries(dummy2D PUBLIC Metis::Metis)
target_link_libraries(dummy2D PUBLIC Hypre::Hypre)
target_link_libraries(dummy2D PUBLIC Precice::Precice)
target_link_libraries(dummy2D PUBLIC ${LAPACK_LIBRARIES})

if(USE_ADIOS2)
  target_link_libraries(flux2D PUBLIC adios2::cxx11_mpi)
  target_link_libraries(dummy2D PUBLIC adios2::cxx11_mpi)
endif()

