!********************************************************************************
!* Educationally-Designed Unstructured 3D (EDU3D) Code
!*
!*  --- EDU3D Euler
!*
!*
!* This module containes a linear relaxtion solver. Currently, the sequential
!* Gauss-Seidel scheme is implemented. Other schemes may be implemented and
!* explored: e.g., multi-color GS, symmetric GS, Jacobi, etc.
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

 module module_linear_solver

  implicit none

 !This module contains the following subroutine:

  public :: linear_relaxation ! perform linear relaxations.

 contains

!********************************************************************************
!
! This subroutine relaxes the linear system by Sequential Gauss-Seidel method.
!
! - Other methods are possible: multi-color GS, symmetric GS, etc.
! - This subroutine is used in the preconditioner for the JFNK solver.
! - Convergence results will be printed in the file: fort.1000.
!
!********************************************************************************
 subroutine linear_relaxation(jac,res, du)

  use module_common_data    , only : p2, zero, one
  use module_ccfv_data_grid , only : ncells, cell

 use data_module    , only : du, res, jac_off, jac_dinv
 use data_module    , only : lrelax_sweeps_actual, lrelax_roc

 use input_parameter, only : lrelax_sweeps, lrelax_tolerance

 implicit none

 !Custom data type for Jacobian
  type jac_data_type
   integer, dimension(:,:)  , pointer :: diag  !Diagonal block
   integer, dimension(:,:,:), pointer :: off   !Off-diagonal block
   integer, dimension(:,:)  , pointer :: inverse_diag !Inverse of diagonal block
  end type jac_data_type

!Input
 real(p2)          , dimension(ncells,4), intent(in) :: res
 type(cc_data_type), dimension(ncells)  , intent(in) :: jac

!Input
 real(p2), dimension(ncells,4),  intent(out) :: du     !correction
 integer                      ,  intent(out) :: irelax !# of relaxations performed.
 real(p2)                     ,  intent(out) :: roc    !Residual reduction from the initial.


!Local variables
 real(p2), dimension(4)   :: b
 integer                  :: i, k, inghbr

!Linear system residual: linear_res = A^{-1}*(b-Ax) at each cell.
 real(p2), dimension(4  ) :: linear_res

!Residual norms(L1,L2,Linf) for the linear system
 real(p2), dimension(4,3) :: linear_res_norm
 real(p2), dimension(4,3) :: linear_res_norm_initial

  lrelax_sweeps_actual = 0      ! to count actual number of relaxations

!Note: write(1000,*) -> This writes in a file named 'fort.1000'.
  write(1000,*)
  write(1000,*) " GS linear relaxaton scheme "
  write(1000,*)

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
! 0. Initialize the correction
!--------------------------------------------------------------------------------

     du = zero

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
! 1. Initial linear res.
!--------------------------------------------------------------------------------

    do i = 1, nnodes

    ! Form the right hand side of GS: b = [ sum( off_diagonal_block*du ) - residual ]
    ! which is just b = -residual since du=0 initially.

        b = -res(i,:)

       ! To compute the linear residual norm: A^{-1}*(b-Ax) at each cell.

                    linear_res = matmul( jac(i)%inverse_diag(:,:), b ) - du(i,:)

          linear_res_norm(:,1) = linear_res_norm(:,1) + abs( linear_res )

    end do

   !Initial linear residual norm

    linear_res_norm(:,1) = linear_res_norm(:,1) / real(nnodes, p2)

     write(1000,'(a,i10,a,es12.5)') " before relax ", 0, &
                " max(L1 norm) = ", maxval(linear_res_norm(:,1))

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
! 2. Linear Relaxation (Sweep)
!--------------------------------------------------------------------------------

  irelax = 0

  relax : do

   irelax = irelax + 1
 
   linear_res_norm(:,1) = zero

 !---------------------------------------------------------
 !---------------------------------------------------------
 ! A sequential Gauss-Seidel Relaxation(sweep)

    gs_node_loop : do i = 1, nnodes

    ! Form the right hand side of GS: [ sum( off_diagonal_block*du ) - residual ]!

        b = -res(i,:)
       gs_node_nghbrs : do k = 1, nnghbrs(i)
        inghbr = cell(i)%nghbr(k)
        b = b - matmul( jac(i)%off(k,:,:), du(inghbr,:) )
       end do gs_node_nghbrs

    ! Update du by the GS relaxation:
    !
    ! e.g., for 3 nghbrs, perform the relaxation in the form:
    !
    !                     diagonal block        sum of off-diagonal block contributions
    !       duj = omega*{ [V/dtj+dR/duj]^{-1}*(-[dRj/du1]*du1 -[dRj/du2]*du2 -[dRj/du3]*du3 -Res_j) - duj }

        linear_res = matmul( jac_dinv(i,:,:), b ) - du(i,:)

           du(i,:) = du(i,:) + linear_res

    ! To compute and store the linear residual norm
       linear_res_norm(:,1) = linear_res_norm(:,1) + abs( linear_res )

    end do gs_node_loop

 ! End of A sequential Gauss-Seidel Relaxation(sweep)
 !---------------------------------------------------------
 !---------------------------------------------------------

 ! This is the current linear residual norm (L1): |res| = |A^{-1}*(b-A*x)|.

    linear_res_norm(:,1) = linear_res_norm(:,1) / real(nnodes, p2)

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
! 3. Check the linear residual.
!--------------------------------------------------------------------------------

 !--------------------------------------------------------------------
 !  After the first relaxation

    if (irelax==1) then

     !-----------------------------------------------------------------
     ! Print the initial linear residual norm in the file 'fort.1000'.

     linear_res_norm_initial = linear_res_norm
     write(1000,'(a,i10,a,es12.5)') "  after relax ", irelax, &
                " max(L1 norm) = ", maxval(linear_res_norm(:,1))

     if ( maxval(linear_res_norm(:,1)) < my_eps ) then
      write(1000,*) " Machine zero res reached. Exit GS relaxation. Total sweeps = ", irelax
      write(1000,*) " tolerance_linear = ", lrelax_tolerance
      lrelax_sweeps_actual = irelax
      exit relax
     endif

 !--------------------------------------------------------------------
 !  After the second relaxation

    else

     !Convergence rate: (current linear res)/(initial linear res)
      roc = maxval(linear_res_norm(1:nq,1)/linear_res_norm_initial(1:nq,1))

     !-----------------------------------------------------------------
     ! Print the current linear residual norm in the file 'fort.1000'.

       if (roc < one) then

        write(1000,'(a,i10,a,es12.5,2x,2x,f8.3)')   "  after relax ", irelax,          &
                  " max(L1 norm) = ", maxval(linear_res_norm(1:nq,1)), roc
       else

        write(1000,'(a,i10,a,es12.5,2x,f6.3,2x,a)') "  after relax ", irelax,          &
                  " max(L1 norm) = ", maxval(linear_res_norm(1:nq,1)), roc," <- diverge"
       endif


     !-----------------------------------------------------------------
     ! Tolerance met: Exit

      if (roc < lrelax_tolerance) then

       write(1000,*)
       write(1000,*) " Tolerance met. Exit GS relaxation. Total sweeps = ", irelax
       write(1000,*) " tolerance_linear = ", lrelax_tolerance

       exit relax

     !-----------------------------------------------------------------
     ! If tolerance is NOT met.

      else

       !---------------------------------------------
       ! Stop if the linear residual is too small.

        if ( maxval(linear_res_norm(1:nq,1)) < my_eps ) then

         write(1000,*)
         write(1000,*) " Residuals too small. Exit GS relaxation. Total sweeps = ", irelax
         write(1000,*) " maxval(linear_res_norm(1:nq,1)) = ", maxval(linear_res_norm(1:nq,1))

         exit relax

        endif

       !-------------------------------------------------
       ! Stop if the maximum number of sweeps is reached.

        if (irelax == lrelax_sweeps) then

         write(1000,*)
         write(1000,*) " Tolerance not met... sweeps = ", lrelax_sweeps

        endif

      endif

     ! End of Tolerance met or NOT met.
     !-----------------------------------------------------------------

    endif

!--------------------------------------------------------------------
! End of 3. Check the linear residual.
!--------------------------------------------------------------------

  end do relax

 return

 end subroutine linear_relaxation
!********************************************************************************



 end module module_linear_solver
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
