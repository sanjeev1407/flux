# FLUX
  Finite-element soLver for Unsteady elctromagnetiX

Building FLUX requires MPI and the packages MFEM, LAPACK, METIS, HYPRE, and PRECICE. Additionally, ADIOS2 can be included in the compilation via the CMake flag `-D USE_ADIOS2`. Standard CMake flags can be used during compilation, such as `-D Mfem_DIR=/path/to/mfem`, if the libraries are not in standard locations for your operating system. Example flags for compiling FLUX : cmake ../ -D CMAKE_BUILD_TYPE=Release -D Hypre_DIR=/path/to/hypre -D Metis_DIR=/path/to/metis -D Mfem_DIR=/path/to/mfem -D Precice_DIR=/path/to/precice

Currently, only the executable `flux2D` is built, which is the 2D axi-symmetric solver. Precice can be used with it via the command line options. An example of available command line options : mpirun -np 4 flux2D -m ../mesh/icp_torch_chamber_farfield.msh -p precice-config_efield.xml -cr 0.08 -cl 0.5 -f 0.37e6 -pdt 0.00000001 -pr 100 -power 150000.0 -nc 6 -co_r 0.109 -co_loc '0.127 0.177 0.227 0.277 0.327 0.377'
Explanation:
            
             -p : flag for activating precice coupling 
             -cr : cut-off radius upto which volume coupling will be done
             -cl : cut-off length upto which volume coupling will be done
             -f : frequency of inductor coils in Hz
             -pdt : precice window time size in seconds
             -pr : frequency of printing output (100 mean output will be printed after every 100 precice coupling time steps)
             -power : power dissipated by inductor coils in Watts
             -o : order of finite elements for Electric field (By default it is 1)
             -nc : number of coils
             -co_r : coil radius in metres
             -co_loc : locations of coils in metres (e.g. -co_loc '116.2e-3 124.2e-3 132.2e-3')

The subset of the domain that will be communicated is determined by the C++ lambda passed to the function BuildPreciceToMfemMap. See the documentation of this function in `flux/flux2D.hpp` for more information.

Note: When using preCICE, the window time size must be given via the command line option `-pdt`. This is the size of the fake timestep FLUX will take, essentially leading to one call to flux per every preCICE communication.
