#
# Compile and run the program to generate a grid.
#

gfortran edu2d_grid_cylinder.f90

./a.out

# or

./a.exe

# and you'll use cylinder_080x040_tria.grid
# copy it to the above directory as 
# "cylinder.grid",

cp cylinder_080x040_tria.grid ../cylinder.grid

# You will put -> project_name = "cylinder"
# in the input.nml file and the CFD code will
# read cylinder.grid.

