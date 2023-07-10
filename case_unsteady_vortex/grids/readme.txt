#
# Compile and run the program to generate a grid.
#

gfortran hw1_generate_grids.f90

./a.out

# or

./a.exe

# and you'll use grid_129x129_tria.grid
# copy it to the above directory as 
# "vortex.grid",

cp grid_129x129_tria.grid ../vortex.grid

# You will put -> project_name = "vortex"
# in the input.nml file and the CFD code will
# read vortex.grid.

