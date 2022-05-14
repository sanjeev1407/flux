#bin/bash

$1/dummy2D -m icp_2d_farfield.msh -p precice_config.xml -cr 0.08 -cl 0.5 -t 5 &
$1/flux2D -m icp_2d_farfield.msh -p precice_config.xml -cr 0.08 -cl 0.5 -pdt 1.0 -freq 0.45e6 -power 150000.0
