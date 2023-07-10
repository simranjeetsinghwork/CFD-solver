###################################################
#
#
#
# To compile, just type make. This will use Makefile.

%make

# To compile with Makefile_debug, specify the name by -f

%make -f Makefile_debug

# To run, type edu2d

%edu2d

# or

%./edu2d

# The code will read example.grid and example.bc and
# write out a Tecplot file 'example_tec.dat' for viewing.


# To clean up .o .mod and the executable, type "make clean".

%make clean

