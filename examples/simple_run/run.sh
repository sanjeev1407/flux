#bin/bash

mpiexec -np $2 $1 -m icp_2d_farfield.msh -o 1 -freq 0.37e6 -power 150000.0
