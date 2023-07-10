!********************************************************************************
! This program generates 2D quad and triangular grids in a rectangular domain.
!
! Example:
!
!     o--------o--------o--------o--------o
!     |        |        |        |        |
!     |        |        |        |        |
!   16|      17|      18|      19|      20|           nx = 5
!     o--------o--------o--------o--------o           ny = 5
!     |        |        |        |        |
!     |        |        |        |        |   nnodes = 5x5 = 25
!   11|      12|      13|      14|      15|
!     o--------o--------o--------o--------o
!     |        |        |        |        |
!     |        |        |        |        |
!    6|       7|       8|       9|      10|
!     o--------o--------o--------o--------o
!     |        |        |        |        |
!     |        |        |        |        |
!    1|       2|       3|       4|       5|
!     o--------o--------o--------o--------o
!
!
!
! 3-Step Generation:
!
! 1. Generate a temporary structured grid data for nodes: xs(i,j) and ys(i,j)
! 2. Generate a 1D node array: x(1:nnodes), y(1:nnodes)
! 3. Generate element connectivity data: tria(1:ntria,3), quad(1:nquad,4)
!
!
! Inuput: 
!        xmin, xmax = x-coordinates of the left and right ends
!        ymin, ymax = y-coordinates of the bottom and top ends
!                nx = number of nodes in x-direction
!                ny = number of nodes in y-direction
!
!        (NOTE: All input parameters are defined inside the program.)
!
! Output: Aseries of grids with
!
!        tria_grid_tecplot.dat = tecplot file of the triangular grid
!        quad_grid_tecplot.dat = tecplot file of the quadrilateral grid
!                tria_grid.vtk = .vtk file of the triangular grid
!                quad_grid.vtk = .vtk file of the quadrilateral grid
!                    tria.grid = triangular-grid file for a solver
!                    quad.grid = quadrilateral-grid file for a solver
!                   project.bc = file that contains boundary condition info
!
!        (NOTE: The grid file format is specific to the EDU2D solvers.)
!
!
! To compile: 'gfortran hw1_generate_grids.f90', and run by typing 'a.out'.
!
!
!        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!
! the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!
! This is Version 0 (September 2018).
!
! This F90 code is written and made available for an educational purpose as well
! as for generating grids for the OSSAN-Euler2D code.
! This file may be updated in future.
!
! Katate Masatsuka,  http://www.cfdbooks.com
!********************************************************************************
 program twod_rectangular_grid

 implicit none

!Parameters
  integer, parameter :: sp = kind(1.0) !Single precision
  integer, parameter :: p2 = selected_real_kind(2*precision(1.0_sp)) !Double precision

 !real(p2)=double precision.

 real(p2) :: zero=0.0_p2, half=0.5_p2, one=1.0_p2

!Input  - domain size and grid dimensions
 real(p2) :: xmin, xmax !Minimum x and Max x.
 real(p2) :: ymin, ymax !Minimum y and Max y
 integer  :: nx         !Number of nodes in x-direction
 integer  :: ny         !Number of nodes in y-direction

 integer :: iperturb, irandom, ismooth

!Output - grid files
 character(80) :: datafile_tria_tec = "tria_grid_tecplot.dat"
 character(80) :: datafile_quad_tec = "quad_grid_tecplot.dat"
 character(80) :: datafile_tria_vtk = "tria_grid.vtk"
 character(80) :: datafile_quad_vtk = "quad_grid.vtk"
 character(80) :: datafile_tria     = "tria.grid" !Triangular grid file for a solver
 character(80) :: datafile_quad     = "quad.grid" !      Quad grid file for a solver
 character(80) :: datafile_bc       = "project.bc"!BC file for a solver


!Local variables
 real(p2), dimension(:,:), allocatable :: xs, ys !Structured grid data

 integer :: nnodes !Total number of nodes
 integer ::  ntria !Total number of triangles
 integer ::  nquad !Total number of quadrilaterals
 integer ::  inode !Local variables used in the 1D nodal array

 integer , dimension(:,:), allocatable :: tria !Triangle connectivity data
 integer , dimension(:,:), allocatable :: quad !Quad connectivity data
 real(p2), dimension(:)  , allocatable :: x, y !Nodal coordinates, 1D array

 real(p2) :: dx !Uniform grid spacing in x-direction = (xmax-xmin)/nx
 real(p2) :: dy !Uniform grid spacing in y-direction = (ymax-ymin)/ny
 integer  :: i, j, os

 real(p2) :: rn, ar

 integer :: igrid

!  Define the grid size: the number of nodes in each direction.
!  NOTE: "ny" is better to be an odd number to place a node at the midpoint in y.

   iperturb = 0  !=0 no perturbation  , >0 perturb nodal coordinates
   irandom  = 0  !=0 regular splitting, >0 random splitting
   ismooth  = 0  !=0 no smoothing     , >0 smooth node distribution by averaging

!Generate grids.

 grids: do igrid = 5, 5

  if     (igrid==1) then
   nx = 9
   datafile_tria     = "grid_009x009_tria.grid"
   datafile_quad     = "grid_009x009_quad.grid"
   datafile_tria_tec = "grid_009x009_tria_tec.dat"
   datafile_quad_tec = "grid_009x009_quad_tec.dat"
   datafile_tria_vtk = "grid_009x009_tria.vtk"
   datafile_quad_vtk = "grid_009x009_quad.vtk"
  elseif (igrid==2) then
   nx = 17
   datafile_tria     = "grid_017x017_tria.grid"
   datafile_quad     = "grid_017x017_quad.grid"
   datafile_tria_tec = "grid_017x017_tria_tec.dat"
   datafile_quad_tec = "grid_017x017_quad_tec.dat"
   datafile_tria_vtk = "grid_017x017_tria.vtk"
   datafile_quad_vtk = "grid_017x017_quad.vtk"
  elseif (igrid==3) then
   nx = 33
   datafile_tria     = "grid_033x033_tria.grid"
   datafile_quad     = "grid_033x033_quad.grid"
   datafile_tria_tec = "grid_033x033_tria_tec.dat"
   datafile_quad_tec = "grid_033x033_quad_tec.dat"
   datafile_tria_vtk = "grid_033x033_tria.vtk"
   datafile_quad_vtk = "grid_033x033_quad.vtk"
  elseif (igrid==4) then
   nx = 65
   datafile_tria     = "grid_065x065_tria.grid"
   datafile_quad     = "grid_065x065_quad.grid"
   datafile_tria_tec = "grid_065x065_tria_tec.dat"
   datafile_quad_tec = "grid_065x065_quad_tec.dat"
   datafile_tria_vtk = "grid_065x065_tria.vtk"
   datafile_quad_vtk = "grid_065x065_quad.vtk"
  elseif (igrid==5) then
   nx = 129
   datafile_tria     = "grid_129x129_tria.grid"
   datafile_quad     = "grid_129x129_quad.grid"
   datafile_tria_tec = "grid_129x129_tria_tec.dat"
   datafile_quad_tec = "grid_129x129_quad_tec.dat"
   datafile_tria_vtk = "grid_129x129_tria.vtk"
   datafile_quad_vtk = "grid_129x129_quad.vtk"
  endif

! In this program, we set ny=nx.

   ny = nx

!--------------------------------------------------------------------------------
! 0. Define the grid size and allocate the structured grid data array.
!
! ymax --------------------------------------
!      .                                    .
!      .                                    .
!      .                                    .
!      .                                    .
!      .                                    .
!      .                                    .
!      .                                    .
!      .                                    .
!      .                                    .
!      .                                    .
!      .                                    .
!      .                                    .
!      .                                    .
!      .                                    .
!      .                                    .
! ymin --------------------------------------
!    xmin                                  xmax

!  Define the domain: here we define a unit square.

      xmin = -20.0_p2
      xmax =  10.0_p2

      ymin = -20.0_p2
      ymax =  10.0_p2

!  Allocate arrays.
   allocate(xs(nx,ny),ys(nx,ny))

!--------------------------------------------------------------------------------
! 1. Generate a structured 2D grid data, (i,j) data: go up in y-direction!
!
! j=5 o--------o--------o--------o--------o
!     |        |        |        |        |
!     |        |        |        |        |   On the left is an example:
!     |        |        |        |        |         nx = 5
! j=4 o--------o--------o--------o--------o         ny = 5
!     |        |        |        |        |
!     |        |        |        |        |
!     |        |        |        |        |
! j=3 o--------o--------o--------o--------o
!     |        |        |        |        |
!     |        |        |        |        |
!     |        |        |        |        |
! j=2 o--------o--------o--------o--------o
!     |        |        |        |        |
!     |        |        |        |        |
!     |        |        |        |        |
! j=1 o--------o--------o--------o--------o
!     i=1      i=2      i=3      i=4      i=5

     write(*,*) "Generating structured data..."

!  Compute the grid spacing in x-direction
    dx = (xmax-xmin)/real(nx-1)

!  Compute the grid spacing in y-direction
    dy = (ymax-ymin)/real(ny-1)

!  Generate nodes in the domain.

   do j = 1, ny  ! Go up in y-direction.
    do i = 1, nx ! Go to the right in x-direction.

     xs(i,j) = xmin + dx*real(i-1)
     ys(i,j) = ymin + dy*real(j-1)

    end do
   end do

! Perturbation

 if (iperturb==1) then

  do j = 2, ny-1   !Go up in y-direction.
   do i = 2, nx-1  !Go to the right in x-direction.

    call random_number(rn)
    xs(i,j) = xs(i,j) + (rn-0.5_p2)*( xs(i+1,j)-xs(i-1,j) )*0.37_p2
    ys(i,j) = ys(i,j) - (rn-0.5_p2)*( ys(i,j+1)-ys(i,j-1) )*0.37_p2

   end do
  end do

 endif

!--------------------------------------------------------------------------------
! 2. Generate unstructured data: 1D array to store the node information.
!
!    - The so-called lexcographic ordering -
!
!   21       22       23       24       25
!     o--------o--------o--------o--------o
!     |        |        |        |        |
!     |        |        |        |        |   On the left is an example:
!   16|      17|      18|      19|      20|           nx = 5
!     o--------o--------o--------o--------o           ny = 5
!     |        |        |        |        |
!     |        |        |        |        |   nnodes = 5x5 = 25
!   11|      12|      13|      14|      15|
!     o--------o--------o--------o--------o
!     |        |        |        |        |
!     |        |        |        |        |
!    6|       7|       8|       9|      10|
!     o--------o--------o--------o--------o
!     |        |        |        |        |
!     |        |        |        |        |
!    1|       2|       3|       4|       5|
!     o--------o--------o--------o--------o
!
   write(*,*) "Generating 1D node array for unstructured grid data..."

!  Total number of nodes
   nnodes = nx*ny

!  Allocate the arrays
   allocate(x(nnodes),y(nnodes))

! Node data: the nodes are ordered in 1D array.

  do j = 1, ny   !Go up in y-direction.
   do i = 1, nx  !Go to the right in x-direction.

    inode = i + (j-1)*nx   !<- Node number in the lexcographic ordering
      x(inode) =   xs(i,j)
      y(inode) =   ys(i,j)

   end do
  end do

! Rescaling: ar = aspect ratio. E.g., ar=100 very thin domain and grid.

    ar = 1.0_p2 !Let's use aspect ratio = 1, i.e., a square.

   do i = 1, nnodes
     x(i) = x(i)      !No change in x.
     y(i) = y(i) / ar !Rescale y.
   end do

! Deallocate the structured data: xs and ys, which are not needed any more.
! - You guys helped me create the 1D array. Thanks!
  deallocate(xs, ys)

 write(*,*)
 write(*,*) " Nodes have been generated:"
 write(*,*) "       nx  =", nx
 write(*,*) "       ny  =", ny
 write(*,*) "    nx*ny  =", nx*ny
 write(*,*) "    nnodes =", nnodes
 write(*,*)
 write(*,*) " Now, generate elements..."
 write(*,*)

!--------------------------------------------------------------------------------
! 3. Generate unstructured element data:
!
!    We generate both quadrilateral and triangular grids.
!    Both grids are constructed in the unstructured (finite-element) data.
!

! Allocate arrays of triangular and quad connectivity data.

!   Number of quadrilaterals = (nx-1)(ny-1)

    allocate( quad((nx-1)*(ny-1)  ,4) )  !(nx-1) cells in x-direction, (ny-1) cells in y-direction

!   Number of triangles = 2*(nx-1)*(ny-1)

    allocate(  tria(2*(nx-1)*(ny-1),3) ) !Twice the # of quads since we split a quad into 2 triangles.

!--------------------------------------------------------------------------------
! (1)Genearte a triangular grid

  write(*,*) "Generating triangular grid..."
  call generate_tria_grid

   if (ismooth > 0) call smoothing !To make it more unstructured and irregular.

  write(*,*)
  write(*,*) " Number of triangles =", ntria
  write(*,*)

  write(*,*) "Writing a tecplot file for the triangular grid..."
  call write_tecplot_file(datafile_tria_tec)
  write(*,*) " --> File generated: ", datafile_tria_tec

  write(*,*) "Writing a .vtk file for the triangular grid..."
  call write_vtk_file(datafile_tria_vtk)
  write(*,*) " --> File generated: ", datafile_tria_vtk

  write(*,*) "Writing a grid file for the triangular grid..."
  call write_grid_file(datafile_tria)
  write(*,*) " --> File generated: ", datafile_tria

  write(*,*)

!--------------------------------------------------------------------------------
! (2)Generate a quadrilateral grid

  write(*,*) "Generating quad grid..."
  call generate_quad_grid
  write(*,*)
  write(*,*) " Number of quads =", nquad
  write(*,*)
  write(*,*) "Writing a tecplot file for the quadrilateral grid..."
  call write_tecplot_file(datafile_quad_tec)
  write(*,*) " --> File generated: ", datafile_quad_tec

  write(*,*) "Writing a .vtk file for the quadrilatera grid..."
  call write_vtk_file(datafile_quad_vtk)
  write(*,*) " --> File generated: ", datafile_quad_vtk

  write(*,*) "Writing a grid file for the quadrilateral grid..."
  call write_grid_file(datafile_quad)
  write(*,*) " --> File generated: ", datafile_quad

!--------------------------------------------------------------------------------
! (3)Write a boundary condition file: to be read by a solver.
!    For now, let's say 'dirichlet' conditions. You can edit this file for
!    a specific problem you're trying to solve, and probably need to change
!    the file name, also.

  write(*,*) "Writing a BC file... ", trim(datafile_bc)

  open(unit=3, file=datafile_bc, status="unknown", iostat=os)
   write(*,*) "Generating bcmap file..."
   write(3,*) "Boundary Segment  Boundary Condition"
   write(3,*) "               1          dirichlet  !Left"
   write(3,*) "               2          dirichlet  !Bottom"
   write(3,*) "               3          dirichlet  !Right"
   write(3,*) "               4          dirichlet  !Top"
  close(3)

!--------------------------------------------------------------------------------

 deallocate(x,y,quad,tria)


 end do grids

 write(*,*)
 write(*,*) "Successfully completed. Stop."

 stop

 contains



!********************************************************************************
! Smooth the nodal coordinates of a triangular grid.
!********************************************************************************
 subroutine smoothing

 implicit none

!Local variables
 integer  :: i
 real(p2) :: zero = 0.0_p2, omega
 integer  :: v(4),k

 real(p2), dimension(:), allocatable :: dx, dy, total

 allocate(dx(nnodes), dy(nnodes), total(nnodes))

 !Smoothing factor

  omega = 0.2_p2

 !We average nodal coordinates locally, 'ismooth' times.
  do k = 1, ismooth

      dx = zero
      dy = zero
   total = zero

   do i = 1, ntria

    v(1:3) = tria(i,1:3)

    dx(v(1)) = dx(v(1)) + x(v(2))-x(v(1))+ x(v(3))-x(v(1))
    dy(v(1)) = dy(v(1)) + y(v(2))-y(v(1))+ y(v(3))-y(v(1))

    dx(v(2)) = dx(v(2)) + x(v(3))-x(v(2))+ x(v(1))-x(v(2))
    dy(v(2)) = dy(v(2)) + y(v(3))-y(v(2))+ y(v(1))-y(v(2))

    dx(v(3)) = dx(v(3)) + x(v(1))-x(v(3))+ x(v(2))-x(v(3))
    dy(v(3)) = dy(v(3)) + y(v(1))-y(v(3))+ y(v(2))-y(v(3))

    total(v) = total(v) + 2

   end do

  do i = 1, nnodes

  !Apply smoothign only in the interior nodes:
  if (x(i)>zero .and. x(i) < 1.0_p2) then
   if (y(i)>zero .and. y(i) < 1.0_p2) then

    x(i) = x(i) + omega * dx(i)/total(i)  !Averaging in x.
    y(i) = y(i) + omega * dy(i)/total(i)  !Averaging in y.

   endif
  endif

  end do

 end do


 deallocate(dx, dy, total)

 end subroutine smoothing

!********************************************************************************
! This subroutine generates triangles by constructing the connectivity data.
!********************************************************************************
 subroutine generate_tria_grid
 implicit none
!Local variables
 integer  :: i, j, inode, i1, i2, i3, i4
 real(p2) :: rn

! No quads
 nquad = 0
  quad = 0

! Trianguler grid with right-up diagonals (i.e., / ).
!
!  inode+nx   inode+nx+1     i4      i3
!       o--------o           o--------o
!       |     .  |           |     .  |
!       |   .    |     or    |   .    |
!       | .      |           | .      |
!       o--------o           o--------o
!    inode    inode+1        i1      i2
!
! Triangle is defined by the counterclockwise ordering of nodes.

 ntria = 0

 do j = 1, ny-1
  do i = 1, nx-1

   inode = i + (j-1)*nx

!     Define the local numbers (see figure above)
      i1 = inode
      i2 = inode + 1
      i3 = inode + nx + 1
      i4 = inode + nx

    call random_number(rn)

 !If randome diagonal is not requested, just use this splitting.
  if (irandom==0) then

           ntria = ntria + 1
    tria(ntria,1) = i1
    tria(ntria,2) = i2
    tria(ntria,3) = i3

           ntria = ntria + 1
    tria(ntria,1) = i1
    tria(ntria,2) = i3
    tria(ntria,3) = i4

 !Here is randome diagonal splitting if irandom/=0.
  else

   if (rn < half) then

           ntria = ntria + 1
    tria(ntria,1) = i1
    tria(ntria,2) = i2
    tria(ntria,3) = i3

           ntria = ntria + 1
    tria(ntria,1) = i1
    tria(ntria,2) = i3
    tria(ntria,3) = i4

   else

           ntria = ntria + 1
    tria(ntria,1) = i1
    tria(ntria,2) = i2
    tria(ntria,3) = i4

           ntria = ntria + 1
    tria(ntria,1) = i2
    tria(ntria,2) = i3
    tria(ntria,3) = i4

   endif

  endif

  end do
 end do

 end subroutine generate_tria_grid
!********************************************************************************


!********************************************************************************
! This subroutine generates quads by constructing the connectivity data.
!********************************************************************************
 subroutine generate_quad_grid
 implicit none
!Local variables
 integer :: i, j, inode, i1, i2, i3, i4

! No triangles
  ntria = 0
   tria = 0

!
!  inode+nx   inode+nx+1     i4      i3
!       o--------o           o--------o
!       |        |           |        |
!       |        |     or    |        |
!       |        |           |        |
!       o--------o           o--------o
!     inode   inode+1        i1      i2
!
! Quad is defined by the counterclockwise ordering of nodes.

! Quadrilateral grid

 nquad = 0

 do j = 1, ny-1
  do i = 1, nx-1

   inode = i + (j-1)*nx
!     Define the local numbers (see figure above)
      i1 = inode
      i2 = inode + 1
      i3 = inode + nx + 1
      i4 = inode + nx

!  Order the quad counterclockwise:
          nquad = nquad + 1
   quad(nquad,1) = i1
   quad(nquad,2) = i2
   quad(nquad,3) = i3
   quad(nquad,4) = i4

  end do
 end do

 end subroutine generate_quad_grid
!********************************************************************************



!********************************************************************************
! This subroutine writes a tecplot file.
!********************************************************************************
 subroutine write_tecplot_file(datafile)
 implicit none
 character(80),            intent(in) :: datafile
 integer :: os
!--------------------------------------------------------------------------------
 open(unit=1, file=datafile, status="unknown", iostat=os)
 write(1,*) 'title = "grid"'
 write(1,*) 'variables = "x","y",'
 write(1,*) 'zone N=',nnodes,',E=',ntria+nquad,',ET=quadrilateral,F=FEPOINT'
!--------------------------------------------------------------------------------
 do i = 1, nnodes
  write(1,*) x(i),y(i)
 end do
!--------------------------------------------------------------------------------
!Triangles
 if (ntria > 0) then
  do i = 1, ntria
   write(1,*)  tria(i,1),  tria(i,2), tria (i,3),  tria(i,3) !The last one is a dummy.
  end do
 endif

!Quadrilaterals
 if (nquad > 0) then
  do i = 1, nquad
   write(1,*) quad(i,1), quad(i,2), quad(i,3), quad(i,4)
  end do
 endif
!--------------------------------------------------------------------------------
 close(1)
 end subroutine write_tecplot_file
!********************************************************************************


!*******************************************************************************
! This subroutine writes a .vtk file for the grid whose name is defined by
! filename_vtk.
!
! Use Paraview to read .vtk and visualize it.  https://www.paraview.org
!
! Search in Google for 'vkt format' to learn .vtk file format.
!*******************************************************************************
 subroutine write_vtk_file(datafile)

  implicit none
  character(80),            intent(in) :: datafile
  integer :: i, j, os

  write(*,*)
  write(*,*) "-------------------------------------------------------"
  write(*,*) ' Writing .vtk file = ', trim(datafile)
  write(*,*)

 !Open the output file.
  open(unit=8, file=datafile, status="unknown", iostat=os)

!---------------------------------------------------------------------------
! Header information

  write(8,'(a)') '# vtk DataFile Version 3.0'
  write(8,'(a)')  datafile
  write(8,'(a)') 'ASCII'
  write(8,'(a)') 'DATASET UNSTRUCTURED_GRID'

!---------------------------------------------------------------------------
! Nodal information
!
! Note: These nodes i=1,nnodes are interpreted as i=0,nnodes-1 in .vtk file.
!       So, later below, the connectivity list for tria and quad will be
!       shifted by -1.

   write(8,*) 'POINTS ', nnodes, ' double'

   do j = 1, nnodes
    write(8,'(3es25.15)') x(j), y(j), zero
   end do

!---------------------------------------------------------------------------
! Cell information.

  !CELLS: # of total cells (tria+quad), total size of the cell list.

  write(8,'(a,i12,i12)') 'CELLS ',ntria+nquad, (3+1)*ntria + (4+1)*nquad

  ! Note: The latter is the number of integer values written below as data.
  !           4 for triangles (# of vertices + 3 vertices), and
  !           5 for quads     (# of vertices + 4 vertices).

  !---------------------------------
  ! 2.1 List of triangles (counterclockwise vertex ordering)

   if (ntria > 0) then
                         ! (# of vertices = 3), 3 vertices in counterclockwise
    do i = 1, ntria
     write(8,'(a,4i12)') '3', tria(i,1)-1, tria(i,2)-1, tria(i,3)-1
                              ! -1 since VTK reads the nodes as 0,1,2,3,..., not 1,2,3,..
    end do
   endif

  !---------------------------------
  ! 2.2 List of quads (counterclockwise vertex ordering)

   if (nquad > 0) then
                         ! (# of vertices = 4), 4 vertices in counterclockwise
    do i = 1, nquad
     write(8,'(a,4i12)') '4', quad(i,1)-1, quad(i,2)-1, quad(i,3)-1, quad(i,4)-1
                             ! -1 since VTK reads the nodes as 0,1,2,3,..., not 1,2,3,..
    end do
   endif

!---------------------------------------------------------------------------
! Cell type information.

                                   !# of all cells
  write(8,'(a,i11)') 'CELL_TYPES ', ntria+nquad

  !Triangle is classified as the cell type 5 in the .vtk format.

  if (ntria > 0) then
   do i = 1, ntria
    write(8,'(i3)') 5
   end do
  endif

  !Triangle is classified as the cell type 9 in the .vtk format.

  if (nquad > 0) then
   do i = 1, nquad
    write(8,'(i3)') 9
   end do
  endif

  close(8)

 end subroutine write_vtk_file
!********************************************************************************


!********************************************************************************
! This subroutine writes a grid file to be read by a solver.
! NOTE: Unlike the tecplot file, this files contains boundary info.
!********************************************************************************
 subroutine write_grid_file(datafile)
 implicit none
 character(80),            intent(in) :: datafile
 integer :: os
!--------------------------------------------------------------------------------
 open(unit=1, file=datafile, status="unknown", iostat=os)

!--------------------------------------------------------------------------------
! Grid size: # of nodes, # of triangles, # of quadrilaterals
  write(1,*) nnodes, ntria, nquad

!--------------------------------------------------------------------------------
! Node data
  do i = 1, nnodes
   write(1,*) x(i), y(i)
  end do

!--------------------------------------------------------------------------------
! Triangle connectivity
  if (ntria > 0) then
   do i = 1, ntria
    write(1,*) tria(i,1), tria(i,2), tria(i,3) 
   end do
  endif

! Quad connectivity
  if (nquad > 0) then
   do i = 1, nquad
    write(1,*) quad(i,1), quad(i,2), quad(i,3), quad(i,4) 
   end do
  endif

! Boundary data:
!
!  Example: nx=ny=7
!
!  T=Top, B=Bottom, L=Left, R=Right
!
!    LT----T----T----T----T----T----TR
!     |    |    |    |    |    |    |
!     L----o----o----o----o----o----R
!     |    |    |    |    |    |    |
!     L----o----o----o----o----o----R
!     |    |    |    |    |    |    |
!     L----o----o----o----o----o----R
!     |    |    |    |    |    |    |
!     L----o----o----o----o----o----R
!     |    |    |    |    |    |    |
!     L----o----o----o----o----o----R
!     |    |    |    |    |    |    |
!    LB----B----B----B----B----B----BR
!

! Number of boundary segments
  write(1,*) 4

  write(1,*)  ny         !Inflow
  write(1,*)  nx         !Bottom Outflow
  write(1,*)  ny         !Right  Outflow
  write(1,*)  nx         !Top Wall

  write(1,*)

! Left
  do j = ny, 1, -1
   i = 1
    write(1,*) i + (j-1)*nx
  end do

! Bottom
  do i = 1, nx
   j = 1
    write(1,*) i + (j-1)*nx
  end do

! Right
  do j = 1, ny
   i = nx
    write(1,*) i + (j-1)*nx
  end do

! Top
  do i = nx, 1, -1
   j = ny
    write(1,*) i + (j-1)*nx
  end do

!--------------------------------------------------------------------------------
 close(1)
 end subroutine write_grid_file
!********************************************************************************


 end program twod_rectangular_grid
