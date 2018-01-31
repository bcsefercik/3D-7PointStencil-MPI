# Parallel Computation of 7-Point Stencil via 3D Decomposition - MPI

## Description

Implementation of a parallel algorithm using 7 point stencil computation is performed using MPI.

## Compilation

- 3d_slice.cpp:		$ mpic++ -Wall 3d_slice.cpp -o slice

- 3d_serial.cpp:	$ gcc -o 3d_serial serial.cpp

## Usage

$ mpirun -np "p" ./slice "width" "height" "depth" "iteration"
where: 
- "p": number of processees
- "width": width of the 3d structure
- "height": height of the 3d structure
- "depth": depth of the 3d structure
- ("width"+2) must be divisible by "p" and result must be bigger than 2

e.g. $ mpirun -np 4 ./slice 126 62 62 100
	 $ ./serial 126 62 62 100
