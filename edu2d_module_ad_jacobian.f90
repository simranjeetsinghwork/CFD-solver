!********************************************************************************
!* Educationally-Designed Unstructured 3D (EDU3D) Code
!*
!*  --- EDU3D Euler
!*
!* This module containes subroutines that comptutes the residual Jacobian
!* for the node-centered edge-based discretization. The Jacobian is exact
!* for the first-order node-centered edge-based discretization. So, Newton's
!* method can be constructed for the first-order scheme. This Jacobian can serve
!* as a preconditioning matrix for Jacobian-Free Newton-Krylov solver for
!* higher-order edge-based schemes.
!*
!*
!*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!*
!* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!*
!* 
!* This F90 program is written and made available for an educational purpose.
!*
!* This file may be updated in future.
!*
!* Katate Masatsuka, October 2017. http://www.cfdbooks.com
!********************************************************************************

 module module_jacobian

  implicit none

 !This module contains the following subroutine:

  public :: compute_jacobian ! compute the jacobian

  private :: kth_nghbr        ! find a nghbr index
  private :: gewp_solve       ! Gauss elimination to invert the diagonal blocks

 contains

!********************************************************************************
! Compute the residual Jacobian exactly for the first-order residual.
!
! - Jac = [V/dtau+dRes1/dU] for the linear system: Jac*dU = -Residual.
!
!    where Res1     = 1st-order accurate residual (compact stencil),
!          Residual = 1st/2nd-order residual.
!          V/dtau   = global diagonal matrix of dual_vol/dtau.
!
! - This is the linear system that needs to be solved per nonlinear iteration.
! - This is very similar to the residual subroutine.
! - Exact flux Jacobian will be computed by Automatic Differentiation.
! - Diagonal blocks are inverted and stored at the end.
!
! - Off-diagonal blocks are stored as follows:
!
!     node1 with 5 nghbrs : jac_off( 1,:,:), jac_off( 2,:,:), ..., jac_off( 5,:,:),
!     node2 with 6 nghbrs : jac_off( 6,:,:), jac_off( 7,:,:), ..., jac_off(11,:,:),
!     node3 with 4 nghbrs : jac_off(12,:,:), jac_off(13,:,:), ..., jac_off(15,:,:),
!     etc.
!
!    Loop over nghbrs can be performed by using nnghrs array, which has the
!    number of nghbr nodes for each node: e.g.,
!
!      nnghbrs(node1)=5, nnghbrs(node2)=6, nnghbrs(node3)=4, etc,
!
!    and k0 array, which has the starting index for each node:
!
!      k0(node1)=0, k0(node2)=5, k0(node3)=11, etc.
! 
!    The off-diagonal blocks can be accessed at a node j as
!
!     do k = 1, nnghbrs(j)
!      jac_off(i,k0(j)+k,:,:)
!     end do
!
! - The Jacobian matirx is stored in two parts: diagonal and off-diagonal blocks.
!
!       For example, consider the residual at node j, Rj, which is a vector,
!       and suppose j has 5 neighbor nodes k=(1, 2, 3, 4, 5).
!       For implicit defect-correction solver, the residual is taken as 1st-order scheme,
!       so that Res_j is a function of Uj, U1, U2, U3, U4, and U5 (compact stencil).
!
!       Then, the j-th component of [V/dt+dR/dU]*dU = -Residual is
!
!       [        dRj/dU1]*dU1 + [dRj/dU2]*dU2 + [dRj/dU3]*dU3 +
!       [V/dtauj+dRj/dUj]*dUj + [dRj/dU4]*dU4 + [dRj/dU5]*dU5 = - Rj
!
!       where dtauj is the local time step at j, dRj/dUj is the derivative of Rj
!       with respect to the solution at j, i.e., Uj, and similarly for others.
!       The Gauss-Seidel relaxation can be written as
!
!       dUj = [V/dtj+dR/dUj]^{-1}*(- [dRj/dU1]*dU1 - [dRj/dU2]*dU2
!                                  - [dRj/dU3]*dU3 - [dRj/dU4]*dU4
!                                  - [dRj/dU5]*dU5 - Res_j )
!
!       To perform this operation, we store two types of data:
!
!        1.  5x5     Diagonal block:  jac_diag(j, :,:) = [V/dtauj+dRj/dUj]
!        2.  5x5 Off-Diagonal Block:  jac_off(j,k,:,:) = [        dRj/dUk], k = 1,2,3,4,5
!
!       so, we can perform the relaxation at every node by looping over the neighbors:
!
!         b = -Res_j
!        do k = 1, nnghbrs(j)
!         b = b - jac_off(i,k0(i)+k,:,:)*dU( nghbr(k0(i)+k), :)
!        end do
!         dUj = jac_diag(j, :,:)^{-1}*b,
!
!       where nghbr(k0(i)+k) stores the actual node number for kth nghbr.
!
!       Actual relaxation will be performed in edu3d_euler_linear_solver.f90.
!
!       To compute dR/dUj and dR/dUk, just like the residual is computed as a sum of
!       numerical fluxes over edges, we compute the Jacobian blocks as a sum of
!       derivatives of the numerical flux over edges. So, we loop over edges,
!       compute the derivatives of the numerical flux with respect to each end node
!       of the edge, and accumulate the contribution at the two nodes.
!       At the end of the edge loop, at every node j, we have
!
!               dRj/dUj = sum_over_edges ( dflux/dUj )
!               dRj/dUk = dflux/dUk
!
!       which forms the j-th row of the global Jacobian. Here, 'flux' is the
!       numerical flux between nodes j and k.
!
! - The flux Jacobian, dflux/dUj, is computed by automatic differentiation.
!   Numerical flux subroutines in edu3d_flux_functions_ddt.f90 are all written
!   by automatic differentiation. They return the flux vector as well as
!   the derivative. This is very convenient, and it makes it so easy to try out
!   a new numerical flux without coding the Jacobian subroutine.
!   See edu3d_flux_functions_ddt.f90 and edu3d_euler_fluxjac.f90 to learn
!   how it is implemented.
!
!********************************************************************************
 subroutine compute_jacobian(jac)

 use module_common_data     , only : p2, zero, half, one

 use module_common_data     , only : x, y, bc_type

 use module_ccfv_data_grid  , only : nfaces, face, ncells, cell, bound, nbound, face_nrml, face_nrml_mag

 use module_ccfv_data_soln  , only : res, u, gradw, wsn

 use module_flux            , only : interface_flux
 use module_bc_states       , only : get_right_state
 use module_ccfv_gradient   , only : compute_gradients
 use module_input_parameter , only : second_order, use_limiter
 use module_ccfv_limiter    , only : compute_limiter, phi
 use module_ccfv_data_soln , only : dtau

 implicit none

 !Custom data type for Jacobian
  type jac_data_type
   integer, dimension(:,:)  , pointer :: diag  !Diagonal block
   integer, dimension(:,:,:), pointer :: off   !Off-diagonal block
   integer, dimension(:,:)  , pointer :: inverse_diag !Inverse of diagonal block
  end type jac_data_type

 type(jac_data_type), dimension(ncells), intent(out) :: jac

!Local variables

 real(p2)                 :: xm, ym
 integer                  :: i
 integer                  :: c1, c2

 real(p2), dimension(4)   :: u1, u2
 real(p2), dimension(2)   :: unit_face_normal

 real(p2), dimension(4)   :: num_flux
 integer                  :: j, ib, v1, v2
 real(p2), dimension(4)   :: ub
 real(p2)                 :: wave_speed
 real(p2)                 :: phi1, phi2

 real(p2), dimension(4,4) :: duRduL
 real(p2), dimension(4,4) :: dFnduL, dFnduR

 integer                  :: idestat
 real(p2), dimension(4,4) :: inverse

! Initialization

  do i = 1, ncells
   jac(i)%diag = zero
   jac(i)%off  = zero
  end do

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
! Residual computation: interior faces

!--------------------------------------------------------------------------------
! Flux computation across internal faces (to be accumulated in res(:))
!
!        o---o---------o       face(j,:) = [i,k,v2,v1]
!       .    .          .
!      .     .           .
!     .      .normal      .
!    .  Left .--->  Right  .
!   .   c1   .       c2     .
!  .         .               .
! o----------o----------------o
!
!
! 1. Extrapolate the solutions to the face-midpoint from centroids 1 and 2.
! 2. Compute the numerical flux.
! 3. Add it to the residual for 1, and subtract it from the residual for 2.
!
!--------------------------------------------------------------------------------
  loop_faces : do i = 1, nfaces

 ! Left and right cells of the i-th face

     c1 = face(i,1)  ! Left cell of the face
     c2 = face(i,2)  ! Right cell of theface

     v1 = face(i,3)  ! Left  node of the face
     v2 = face(i,4)  ! Right node of the face

     u1 = u(c1,1:4) !Conservative variables at c1
     u2 = u(c2,1:4) !Conservative variables at c2

   unit_face_normal = face_nrml(i,1:2) !Unit face normal vector: c1 -> c2.

 ! Compute the left and right flux Jacobians.

   call interface_jac_ad(          u1,       u2   , & !<- Left/right states
                                  unit_face_normal, & !<- unit face normal
                                    dFnduL, dFnduR  ) !<- Output

 !  Add the Jacobian multiplied by the magnitude of the directed area vector to node1.

     jac(c1)%diag(:,:)  = jac(c1)%diag(c1,:,:) + (dFnduL) * face_nrml_mag(i)
     k = kth_nghbr(c1,c2)
     jac(c1)%off(k,:,:) = jac(c1)%off(k,:,:)   + (dFnduR) * face_nrml_mag(i)

 !  Subtract the Jacobian multiplied by the magnitude of the directed area vector from node2.
 !  NOTE: Subtract because the outward face normal is -da for the node2.

     jac(c2)%diag(:,:)  = jac(c2)%diag(:,:)  - (dFnduL) * face_nrml_mag(i)
     k = kth_nghbr(c2,c1)
     jac(c1)%off(k,:,:) = jac(c1)%off(k,:,:) - (dFnduR) * face_nrml_mag(i)

  end do loop_edges
!--------------------------------------------------------------------------------


!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
! Residual computation: boundary faces:
!
! Close the residual by looping over boundary faces and distribute a contribution
! to the corresponding cell.

! Boundary face j consists of nodes j and j+1.
!
!  Interior domain      /
!                      /
!              /\     o
!             /  \   /
!            / c1 \ /   Outside the domain
! --o-------o------o
!           j   |  j+1
!               |   
!               v Face normal for the face j.
!
! c = bcell, the cell having the boundary face j.
!

  boundary_part : do ib = 1, nbound

   bface : do j = 1, bound(ib)%nbfaces

     v1 = bound(ib)%bnode(j)
     v2 = bound(ib)%bnode(j+1)

    !Cell having a boundary face defined by the set of nodes j and j+1.
     c1 = bound(ib)%bcell(j)

         u1 = u(c1,1:4)
     unit_face_normal = bound(ib)%bface_nrml(j,1:2)

    call get_right_state_ad(xm,ym, u1, unit_face_normal, bc_type(ib), ub, duRduL)

    call interface_jac_ad(          u1,       ub   , & !<- Left/right states
                                   unit_face_normal, & !<- unit face normal
                                   dFnduL, dFnduR  ) !<- Output


     jac(c1)%diag(:,:)  = jac(c1)%diag(:,:)  + ( dFnduL + matmul(dFnduR,duRduL) ) * face_nrml_mag(i)

   end do bface

  end do boundary_part 

!--------------------------------------------------------------------------------
! Add pseudo time term vol/dtau to the diagonals of the diagonal blocks
! to form: Jacobian = P*{-1}*V/dtau + dRes1/dU, where V/dtau is a global diagonal
! matrix having vol(i)/dtau(i) for each node, and Res1 is the first-order
! residual.

  do i = 1, nnodes
   do k = 1, 4
    jac(i)%diag(k,k) = jac(i)%diag(k,k) + vol(i)/dtau(i)
   end do
  end do

!--------------------------------------------------------------------------------
! Invert the diagonal blocks.
!  - I will invert the diagonal blocks here, and store them.
!  - Is it inefficient? Think about it.

  invert_at_cell : do i = 1, nnodes

     inverse = zero
     idestat = 0

 !  Invert the diagonal block at node i by Gauss elimination with pivoting.
 !  Note: gewp_solve() actually solves a linear system, Ax=b, but here
 !        we use it only to obtain the inverse of A. So, b is a dummy.

 !                    A                dim  A^{-1}
    call gewp_solve( jac(i)%diag(:,:),  4 , inverse, idestat )

 !  Report errors

    if (idestat/=0) then
     write(*,*) " Error in inverting the diagonal block... Stop"
     write(*,*) "  Node number = ", i
     write(*,*) "  Location    = ", x(i), y(i)
     do k = 1, 4
      write(*,'(12(es8.1))') ( jac(i)%diag(k,j), j=1,4 )
     end do
     stop
    endif

 !  Save the inverse in jac_dinv.

    jac(i)%inverse_diag(i,:,:) = inverse

  end do invert_at_cell

    if (debug_mode) stop

 end subroutine compute_jacobian
!********************************************************************************


!********************************************************************************
! This function is useful to find 'k', such that 
! cell 'n2' is the k-th nghbr of the cell 'n1'.
!********************************************************************************
 function kth_nghbr(n1,n2)

 use module_ccfv_data_grid, only : cell

 implicit none

 integer, intent(in) :: n1, n2

 integer :: k
 integer :: kth_nghbr
 logical :: found

  kth_nghbr = 0
      found = .false.

 !Find k, such that n2 is the k-th neghbor of n1.

  n1_nghbr_loop : do k = 1, cell(n1)%nnghbrs

   if ( cell(n1)%nghbr(k) == n2 ) then

    kth_nghbr = k
    found = .true.
    exit n1_nghbr_loop

   endif

  end do n1_nghbr_loop

  if (.not. found) then
   write(*,*) " Neighbor not found.. Error. Stop"
   stop
  endif

 end function kth_nghbr






!****************************************************************************
!* ------------------ GAUSS ELIMINATION WITH PIVOTING ---------------------
!*
!*  This computes the inverse of an (nm)x(nm) matrix "ai" and also
!*  computes the solution to a given lienar system.
!*
!*  IN :       ai = An (nm)x(nm) matrix whoise inverse is sought.
!*             bi = A vector of (nm): Right hand side of the linear sytem
!*             nm = The size of the matrix "ai"
!*
!* OUT :
!*            sol = Solution to the linear system: ai*sol=bi
!*        inverse = the inverse of "ai".
!*       idetstat = 0 -> inverse successfully computed
!*                  1 -> THE INVERSE DOES NOT EXIST (det=0).
!*                  2 -> No unique solutions exist.
!*****************************************************************************
  subroutine gewp_solve(ai,nm, inverse,idetstat)

  implicit none

 integer , parameter ::    p2 = selected_real_kind(15) ! Double precision
 real(p2), parameter ::  zero = 0.0_p2
 real(p2), parameter ::   one = 1.0_p2

  integer ,                   intent( in) :: nm
  real(p2), dimension(nm,nm), intent( in) :: ai

  real(p2), dimension(nm,nm), intent(out) :: inverse
  integer ,                   intent(out) :: idetstat

  real(p2), dimension(nm,nm+1) :: a
  real(p2), dimension(nm)      :: x
  integer , dimension(nm)      :: nrow
  integer                      :: I,J,K,pp,m

  do m = 1, nm
!*****************************************************************************
!* Set up the matrix a
!*****************************************************************************

       do J=1,nm
        do I=1,nm
          a(I,J) = ai(I,J)
        end do
       end do

       do k=1,nm
          a(k,nm+1)=zero; nrow(k)=k
       end do
          a(m,nm+1)=one

!*****************************************************************************
!* Let's go.
!*****************************************************************************
    do j=1,nm-1
!*****************************************************************************
!* FIND SMALLEST pp FOR a(pp,j) IS MAXIMUM IN JTH COLUMN.
!***************************************************************************** 
      call findmax(nm,j,pp,a,nrow)
!*****************************************************************************
!* IF a(nrow(p),j) IS zero, THERE'S NO UNIQUE SOLUTIONS      
!*****************************************************************************
      if (abs(a(nrow(pp),j)) < epsilon(one)) then
       write(6,*) 'THE INVERSE DOES NOT EXIST.'
        idetstat = 1
        return
      endif
!*****************************************************************************
!* IF THE MAX IS NOT A DIAGONAL ELEMENT, SWITCH THOSE ROWS       
!*****************************************************************************
      if (nrow(pp) .ne. nrow(j)) then
      call switch(nm,j,pp,nrow)
      else
      endif  
!*****************************************************************************
!* ELIMINATE ALL THE ENTRIES BELOW THE DIAGONAL ONE
!***************************************************************************** 
      call eliminate_below(nm,j,a,nrow)

    end do
!*****************************************************************************
!* CHECK IF a(nrow(N),N)=0.0 .
!*****************************************************************************
      if (abs(a(nrow(nm),nm)) < epsilon(one)) then
        write(6,*) 'NO UNIQUE SOLUTION EXISTS!'
        idetstat = 2
        return
      else
      endif
!*****************************************************************************
!* BACKSUBSTITUTION!
!*****************************************************************************
      call backsub(nm,x,a,nrow)
!*****************************************************************************
!* STORE THE SOLUTIONS, YOU KNOW THEY ARE INVERSE(i,m) i=1...
!*****************************************************************************
      do i=1,nm
         inverse(i,m)=x(i)
      end do
!*****************************************************************************
  end do

      idetstat = 0

    return

!*****************************************************************************
 end subroutine gewp_solve

!*****************************************************************************
!* Four subroutines below are used in gewp_solve() above.
!*****************************************************************************
!* FIND MAXIMUM ELEMENT IN jth COLUMN 
!***************************************************************************** 
      subroutine findmax(nm,j,pp,a,nrow)

      implicit none

      integer , parameter   :: p2 = selected_real_kind(15) ! Double precision
      integer , intent( in) :: nm
      real(p2), intent( in) :: a(nm,nm+1)
      integer , intent( in) :: j,nrow(nm)
      integer , intent(out) :: pp
      real(p2)              :: max
      integer               :: i

            max=abs(a(nrow(j),j)); pp=j

           do i=j+1,nm

             if (max < abs(a(nrow(i),j))) then

                  pp=i; max=abs(a(nrow(i),j))

             endif

           end do

      return

      end subroutine findmax
!*****************************************************************************
!* SWITCH THOSE ROWS       
!*****************************************************************************
      subroutine switch(nm,j,pp,nrow)

      implicit none

      integer, intent(   in) :: nm,j,pp
      integer, intent(inout) :: nrow(nm)
      integer                :: ncopy

      if (nrow(pp).ne.nrow(j)) then

         ncopy=nrow(j)
         nrow(j)=nrow(pp)
         nrow(pp)=ncopy

      endif

      return

      end subroutine switch
!*****************************************************************************
!* ELIMINATE ALL THE ENTRIES BELOW THE DIAGONAL ONE
!*(Give me j, the column you are working on now)
!***************************************************************************** 
      subroutine eliminate_below(nm,j,a,nrow)

      implicit none

      integer , parameter     :: p2 = selected_real_kind(15) ! Double precision
      real(p2), parameter     :: zero = 0.0_p2
      integer , intent(   in) :: nm
      real(p2), intent(inout) :: a(nm,nm+1)
      integer , intent(   in) :: j,nrow(nm)
      real(p2)                :: m
      integer                 :: k,i

      do i=j+1,nm

        m=a(nrow(i),j)/a(nrow(j),j)
        a(nrow(i),j)=zero

          do k=j+1,nm+1
            a(nrow(i),k)=a(nrow(i),k)-m*a(nrow(j),k)
          end do

      end do

      return

      end subroutine eliminate_below
!*****************************************************************************
!* BACKSUBSTITUTION!
!*****************************************************************************
      subroutine backsub(nm,x,a,nrow)

      implicit none

      integer , parameter   :: p2 = selected_real_kind(15) ! Double precision
      real(p2), parameter   :: zero = 0.0_p2

      integer , intent( in) :: nm
      real(p2), intent( in) :: a(nm,nm+1)
      integer , intent( in) :: nrow(nm)
      real(p2), intent(out) :: x(nm)
      real(p2)              :: sum
      integer               :: i,k

      x(nm)=a(nrow(nm),nm+1)/a(nrow(nm),nm)

      do i=nm-1,1,-1

         sum=zero

           do k=i+1,nm

              sum=sum+a(nrow(i),k)*x(k)

           end do

      x(i)=(a(nrow(i),nm+1)-sum)/a(nrow(i),i)

      end do

      return

      end subroutine backsub
!*********************************************************************


 end module module_jacobian
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
