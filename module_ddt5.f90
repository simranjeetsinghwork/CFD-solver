!********************************************************************************
!* Automatic Differentiation (AD) module, Version 2 (2018),
!*
!*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!*
!* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!*
!* This module contains the derivative data type (ddt) and associated operators,
!* specifically written for storing 1 derivatives (thus, ddt1).
!*
!* To use the data type and operators, declare "use derivative_data".
!* Then, declare a scalar derivative_data_type variable, say 'a', by
!*
!*  type(derivative_data_type) :: a
!*
!* Then, a has components, a%f to store a funciton value, and the array of
!* dimension 1, a%df(1:1) to store derivatives.
!*
!* Once declared, various operators can be applied to these variables, e.g.,
!*
!*  type(derivative_data_type) :: a1, a2, a3
!*   a3 = a1 + a2
!*   a1 = a3**a2 - a2*a1
!*
!* Operators are overloaded as defined in this module.
!* Currently available operators are
!*
!* [ assignment, +, -, *, /, **, ddt_sqrt, ddt_abs, ddt_sign, ddt_max, ddt_min,
!*   <, <=, >, >=, == ]
!*
!* More operators may be added as needed.
!*
!* This module is used in the program ("ad_driver.f90") for computing dg/dx of
!* the vector function g(x)=(g1(x), g2(x), g3(x)), where x=(x1,x2,x3).
!*
!* Note: This module can be easily extended to accomodate more derivatives.
!*       Replace "df1" by "df5", for example, to define the derivative data type
!*       having %df of dimension 5 ( also modify dimension(5) -> dimension(5) ).
!*
!* Katate Masatsuka http://www.cfdbooks.com
!********************************************************************************
 module derivative_data_df5

  implicit none

! Select the precision:
!  selected_real_kind( 5) = single precision
!  selected_real_kind(15) = double precision

  integer           , parameter :: my_precision = selected_real_kind(15)
  real(my_precision), parameter :: zero = 0.0_my_precision
  real(my_precision), parameter :: half = 0.5_my_precision
  real(my_precision), parameter ::  one = 1.0_my_precision

  private

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
! Derivative data type definition
! Note: This data type has 1 function value and derivatives w.r.t. five
!       variables. Intended for use with a vector of 5 variables.
  public :: derivative_data_type_df5

  type derivative_data_type_df5
   real(my_precision)               ::  f ! function value
   real(my_precision), dimension(5) :: df ! df/du
  end type derivative_data_type_df5
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------


! Seeding must be performed always before the derivatives are computed.

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
! Set derivatives w.r.t. which derivatives of a functional are computed.
! Basically, this sets variable(i)%df(i)=1.0
! E.g.,  call ddt_seed(ddt variable)

  public :: ddt_seed
  interface ddt_seed
   module procedure ddt_seed_scalar
   module procedure ddt_seed_vector
  end interface
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------


! Below, various operators are defined for the new data type.

!--------------------------------------------------------------------------------
! Assignment operator
! E.g.,  ddt_variable1 = ddt_variable2
! E.g.,  ddt_variable  = real value
! E.g.,    real value  = ddt_variable

  public :: assignment(=)
  interface assignment(=)
   module procedure d_assign_d
   module procedure d_assign_r
   module procedure r_assign_d
  end interface

!--------------------------------------------------------------------------------
! Addition operator
! E.g.,  ddt_variable1 + ddt_variable2
! E.g.,  ddt_variable  + real value
! E.g.,    real value  + ddt_variable

  public :: operator(+)
  interface operator(+)
   module procedure d_plus_d
   module procedure d_plus_r
   module procedure r_plus_d
  end interface

!--------------------------------------------------------------------------------
! Subtraction operator
! E.g.,    ddt_variable1 - ddt_variable2
! E.g.,    ddt_variable  - real value
! E.g.,      real value  - ddt_variable
! E.g.,  - ddt_variable

  public :: operator(-)
  interface operator(-)
   module procedure d_minus_d
   module procedure d_minus_r
   module procedure r_minus_d
   module procedure   minus_d
  end interface

!--------------------------------------------------------------------------------
! Multiplication operator
! E.g.,    ddt_variable1 * ddt_variable2
! E.g.,    ddt_variable  * real value
! E.g.,      real value  * ddt_variable

  public :: operator(*)
  interface operator(*)
   module procedure d_times_d
   module procedure d_times_r
   module procedure r_times_d
  end interface

!--------------------------------------------------------------------------------
! Division operator
! E.g.,    ddt_variable1 / ddt_variable2
! E.g.,    ddt_variable  / real value
! E.g.,      real value  / ddt_variable

  public :: operator(/)
  interface operator(/)
   module procedure d_divided_by_d
   module procedure d_divided_by_r
   module procedure r_divided_by_d
  end interface

!--------------------------------------------------------------------------------
! Exponent operator
! E.g.,    ddt_variable1 ** ddt_variable2
! E.g.,    ddt_variable  ** (real value)
! E.g.,     (real value) ** ddt_variable
! E.g.,    ddt_variable  ** (integer)

  public :: operator(**)
  interface operator(**)
   module procedure d_power_d
   module procedure d_power_r
   module procedure r_power_d
   module procedure d_power_i
  end interface

!--------------------------------------------------------------------------------
! Square root operator
! E.g.,   ddt_sqrt(ddt_variable)

  public :: ddt_sqrt
  interface ddt_sqrt
   module procedure sqrt_d
  end interface

!--------------------------------------------------------------------------------
! Absolute value operator
! E.g.,  ddt_abs(ddt_variable)

  public :: ddt_abs
  interface ddt_abs
   module procedure abs_d
  end interface

!--------------------------------------------------------------------------------
! Sign operator
! E.g.,  ddt_sign(    one, ddt_variable)
! E.g.,  ddt_sign(rael value, ddt_variable)

  public :: ddt_sign
  interface ddt_sign
   module procedure sign_rd
  end interface

!--------------------------------------------------------------------------------
! Max operator
! E.g.,  ddt_max(ddt_variable1, ddt_variable2)
! E.g.,  ddt_max(ddt_variable1,    real value)
! E.g.,  ddt_max(   real value, ddt_variable2)

  public :: ddt_max
  interface ddt_max
   module procedure max_of_d_and_d
   module procedure max_of_d_and_r
   module procedure max_of_r_and_d
  end interface

!--------------------------------------------------------------------------------
! Min operator
! E.g.,  ddt_min(ddt_variable1, ddt_variable2)
! E.g.,  ddt_min(ddt_variable1,    real value)
! E.g.,  ddt_min(   real value, ddt_variable2)

  public :: ddt_min
  interface ddt_min
   module procedure min_of_d_and_d
   module procedure min_of_d_and_r
   module procedure min_of_r_and_d
  end interface

!--------------------------------------------------------------------------------
! Less-than (logical). It compares the function values.
! E.g.,  ddt_variable1 < ddt_variable2
! E.g.,  ddt_variable  < real value
! E.g.,    real value  < ddt_variable

  public :: operator(<)
  interface operator(<)
   module procedure d_less_than_d
   module procedure d_less_than_r
   module procedure r_less_than_d
  end interface

!--------------------------------------------------------------------------------
! Less-than-equal (logical). It compares the function values.
! E.g.,  ddt_variable1 <= ddt_variable2
! E.g.,  ddt_variable  <= real value
! E.g.,    real value  <= ddt_variable

  public :: operator(<=)
  interface operator(<=)
   module procedure d_less_than_equal_d
   module procedure d_less_than_equal_r
   module procedure r_less_than_equal_d
  end interface

!--------------------------------------------------------------------------------
! Greater-than (logical). It compares the function values.
! E.g.,  ddt_variable1 > ddt_variable2
! E.g.,  ddt_variable  > real value
! E.g.,    real value  > ddt_variable

  public :: operator(>)
  interface operator(>)
   module procedure d_greater_than_d
   module procedure d_greater_than_r
   module procedure r_greater_than_d
  end interface

!--------------------------------------------------------------------------------
! Greater-than-equal (logical). It compares the function values.
! E.g.,  ddt_variable1 >= ddt_variable2
! E.g.,  ddt_variable  >= real value
! E.g.,    real value  >= ddt_variable

  public :: operator(>=)
  interface operator(>=)
   module procedure d_greater_than_equal_d
   module procedure d_greater_than_equal_r
   module procedure r_greater_than_equal_d
  end interface

!--------------------------------------------------------------------------------
! Equal (logical). It compares the function values.
! E.g.,  ddt_variable1 == ddt_variable2
! E.g.,  ddt_variable  == real value
! E.g.,    real value  == ddt_variable

  public :: operator(==)
  interface operator(==)
   module procedure d_equal_d
   module procedure d_equal_r
   module procedure r_equal_d
  end interface

 contains

!*******************************************************************************
! Seeding: Set the diagoanl components to be 1.0, so that derivatives are
!          computed w.r.t. these variables. In a way, set the correct derivative.
!
!  Input: d = the derivative_data_type variable to be seeded.
! Output: d with d%df = 1.0 (only the diagonals, d(f_k)/d(f_k) = 1.0)
!
!*******************************************************************************
  subroutine ddt_seed_scalar(d)

   type(derivative_data_type_df5), intent(inout) :: d

    d%df = one

  end subroutine ddt_seed_scalar
!-------------------------------------------------------------------------------
  subroutine ddt_seed_vector(d)

   type(derivative_data_type_df5), dimension(:), intent(inout) :: d
   integer                                   :: n, i

    n = size(d(1)%df)

    do i = 1, n
     d(i)%df(i) = one
    end do

  end subroutine ddt_seed_vector

!*******************************************************************************
! Assignment
!*******************************************************************************
  pure elemental subroutine d_assign_d(d1,d2)

   type(derivative_data_type_df5),  intent(out) :: d1
   type(derivative_data_type_df5),  intent( in) :: d2

    d1%f  = d2%f
    d1%df = d2%df

  end subroutine d_assign_d
!-------------------------------------------------------------------------------
  pure elemental subroutine d_assign_r(d,r)

   type(derivative_data_type_df5),  intent(out) :: d
   real(my_precision)        ,  intent( in) :: r

    d%f  = r
    d%df = zero

  end subroutine d_assign_r
!-------------------------------------------------------------------------------
  pure elemental subroutine r_assign_d(r,d)

   real(my_precision)        ,  intent(out) :: r
   type(derivative_data_type_df5),  intent( in) :: d

    r = d%f

  end subroutine r_assign_d

!*******************************************************************************
! Addition
!*******************************************************************************
  pure elemental function d_plus_d(d1, d2)

   type(derivative_data_type_df5), intent(in) :: d1, d2
   type(derivative_data_type_df5)             :: d_plus_d

    d_plus_d%f  = d1%f  + d2%f
    d_plus_d%df = d1%df + d2%df

  end function d_plus_d
!-------------------------------------------------------------------------------
  pure elemental function d_plus_r(d, r)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)            , intent(in) :: r
   type(derivative_data_type_df5)             :: d_plus_r

    d_plus_r%f  = d%f  + r
    d_plus_r%df = d%df

  end function d_plus_r
!-------------------------------------------------------------------------------
  pure elemental function r_plus_d(r, d)

   type(derivative_data_type_df5),  intent(in) :: d
   real(my_precision)            ,  intent(in) :: r
   type(derivative_data_type_df5)              :: r_plus_d

    r_plus_d%f  = d%f  + r
    r_plus_d%df = d%df

  end function r_plus_d

!*******************************************************************************
! Subtraction
!*******************************************************************************
  pure elemental function d_minus_d(d1, d2)

   type(derivative_data_type_df5), intent(in) :: d1, d2
   type(derivative_data_type_df5)             :: d_minus_d

    d_minus_d%f  = d1%f  - d2%f
    d_minus_d%df = d1%df - d2%df

  end function d_minus_d
!-------------------------------------------------------------------------------
  pure elemental function d_minus_r(d, r)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)            , intent(in) :: r
   type(derivative_data_type_df5)             :: d_minus_r

    d_minus_r%f  = d%f - r
    d_minus_r%df = d%df

  end function d_minus_r
!-------------------------------------------------------------------------------
  pure elemental function r_minus_d(r, d)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)            , intent(in) :: r
   type(derivative_data_type_df5)             :: r_minus_d

    r_minus_d%f  = r - d%f
    r_minus_d%df =   - d%df

  end function r_minus_d
!-------------------------------------------------------------------------------
  pure elemental function minus_d(d)

   type(derivative_data_type_df5), intent(in) :: d
   type(derivative_data_type_df5)             :: minus_d

    minus_d%f  = -d%f
    minus_d%df = -d%df

  end function minus_d

!*******************************************************************************
! Multiplication
!*******************************************************************************
  pure elemental function d_times_d(d1, d2)

   type(derivative_data_type_df5), intent(in) :: d1, d2
   type(derivative_data_type_df5)             :: d_times_d

    d_times_d%f  = d1%f*d2%f
    d_times_d%df = d1%f*d2%df + d2%f*d1%df

  end function d_times_d
!-------------------------------------------------------------------------------
  pure elemental function d_times_r(d, r)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)            , intent(in) :: r
   type(derivative_data_type_df5)             :: d_times_r

    d_times_r%f  = r*d%f
    d_times_r%df = r*d%df

  end function d_times_r
!-------------------------------------------------------------------------------
  pure elemental function r_times_d(r, d)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)            , intent(in) :: r
   type(derivative_data_type_df5)             :: r_times_d

    r_times_d%f  = r*d%f
    r_times_d%df = r*d%df

  end function r_times_d

!*******************************************************************************
! Division (NOTE: derivative of r is zero: r'=0)
!*******************************************************************************
  pure elemental function d_divided_by_d(d1, d2)

   type(derivative_data_type_df5), intent(in) :: d1, d2
   type(derivative_data_type_df5)             :: d_divided_by_d

    d_divided_by_d%f  = d1%f / d2%f
    d_divided_by_d%df = ( d1%df*d2%f - d2%df*d1%f ) / (d2%f*d2%f)

  end function d_divided_by_d
!-------------------------------------------------------------------------------
  pure elemental function d_divided_by_r(d, r)

   type(derivative_data_type_df5),  intent(in) :: d
   real(my_precision)            ,  intent(in) :: r
   type(derivative_data_type_df5)              :: d_divided_by_r

    d_divided_by_r%f  = d%f  / r
    d_divided_by_r%df = d%df / r

  end function d_divided_by_r
!-------------------------------------------------------------------------------
  pure elemental function r_divided_by_d(r, d)

   real(my_precision),    intent(in) :: r
   type(derivative_data_type_df5),  intent(in) :: d
   type(derivative_data_type_df5)              :: r_divided_by_d

    r_divided_by_d%f  = r / d%f
    r_divided_by_d%df = -d%df*r / (d%f*d%f) ! (r/f)'=(r'*f-f'*r)/f^2=-f'*r/f^2

  end function r_divided_by_d

!*******************************************************************************
! Square root
!*******************************************************************************
  pure elemental function sqrt_d(d)
 
   type(derivative_data_type_df5), intent(in) :: d
   type(derivative_data_type_df5)             :: sqrt_d

    sqrt_d%f = sqrt(d%f)

    avoid_zero_denominator : if ( sqrt_d%f > epsilon(one) ) then
      sqrt_d%df = half * d%df / sqrt_d%f
    else
      sqrt_d%df = half * d%df / epsilon(one)
    end if avoid_zero_denominator

  end function sqrt_d

!*******************************************************************************
! Power: Let F=f1^f2. We want F'.
!  1. Set Z=log(f1^f2)=f2*log(f1) -> Differentiate -> Z'=f2'*log(f1)+f2/f1*f1'
!  2. Since Z=log(F), we have Z'=F'/F, so, F'=Z'*F.
!  1. Hence, F' = f1^f2*( f2'*log(f1) + f2/f1*f1' )
!*******************************************************************************
  pure elemental function d_power_d(d1, d2)

  type(derivative_data_type_df5), intent(in) :: d1, d2
  type(derivative_data_type_df5)             :: d_power_d

   d_power_d%f  = d1%f**d2%f
   d_power_d%df = d1%f**d2%f * ( d2%df*log(d1%f) + d2%f*d1%df/d1%f  )
 
  end function d_power_d
!-------------------------------------------------------------------------------
  pure elemental function d_power_r(d, r)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)            , intent(in) :: r
   type(derivative_data_type_df5)             :: d_power_r

    d_power_r%f  = d%f**r
    d_power_r%df = r*d%f**(r-one)*d%df

  end function d_power_r
!-------------------------------------------------------------------------------
  pure elemental function r_power_d(r, d)

  real(my_precision)            , intent(in) :: r
  type(derivative_data_type_df5), intent(in) :: d
  type(derivative_data_type_df5)             :: r_power_d

   r_power_d%f  = r**d%f
   r_power_d%df = r**d%f * d%df*log(r)
 
  end function r_power_d
!-------------------------------------------------------------------------------
  pure elemental function d_power_i(d, i)

   type(derivative_data_type_df5), intent(in) :: d
   integer                       , intent(in) :: i
   type(derivative_data_type_df5)             :: d_power_i

    d_power_i%f  = d%f**i
    d_power_i%df = real(i,my_precision)*d%f**(i-1)*d%df

  end function d_power_i

!*******************************************************************************
! Absolute value fucntion
!*******************************************************************************
  pure elemental function abs_d(d)

   type(derivative_data_type_df5), intent(in) :: d
   type(derivative_data_type_df5)             :: abs_d

    if ( d%f < zero ) then
      abs_d%f  = -d%f
      abs_d%df = -d%df
    else
      abs_d%f  =  d%f
      abs_d%df =  d%df
    endif

  end function abs_d

!*******************************************************************************
! Sign function
!*******************************************************************************
  pure elemental function sign_rd(r,d)

   real(my_precision)            , intent(in) :: r
   type(derivative_data_type_df5), intent(in) :: d
   type(derivative_data_type_df5)             :: sign_rd

    if ( d%f >= zero) then
      sign_rd%f  =  abs(r)
      sign_rd%df =  zero
    else
      sign_rd%f  = -abs(r)
      sign_rd%df =  zero
    end if

  end function sign_rd

!*******************************************************************************
! Max
!*******************************************************************************
  pure elemental function max_of_d_and_d(d1,d2)
 
   type(derivative_data_type_df5), intent(in) :: d1, d2
   type(derivative_data_type_df5)             :: max_of_d_and_d

    if ( d1%f > d2%f ) then
      max_of_d_and_d = d1
    else
      max_of_d_and_d = d2
    endif

  end function max_of_d_and_d
!-------------------------------------------------------------------------------
  pure elemental function max_of_d_and_r(d,r)

   real(my_precision)            , intent(in) :: r
   type(derivative_data_type_df5), intent(in) :: d
   type(derivative_data_type_df5)             :: max_of_d_and_r

    if ( r > d%f ) then
     max_of_d_and_r%f  = r
     max_of_d_and_r%df = zero
    else
     max_of_d_and_r = d
    endif

  end function max_of_d_and_r
!-------------------------------------------------------------------------------
  pure elemental function max_of_r_and_d(r,d2)

   real(my_precision)            , intent(in) :: r
   type(derivative_data_type_df5), intent(in) :: d2
   type(derivative_data_type_df5)             :: max_of_r_and_d

    if ( r > d2%f ) then
     max_of_r_and_d%f  = r
     max_of_r_and_d%df = zero
    else
      max_of_r_and_d = d2
    endif

  end function max_of_r_and_d

!*******************************************************************************
! Min
!*******************************************************************************
  pure elemental function min_of_d_and_d(d1,d2)

   type(derivative_data_type_df5), intent(in) :: d1, d2
   type(derivative_data_type_df5)             :: min_of_d_and_d

    if ( d1%f < d2%f ) then
      min_of_d_and_d = d1
    else
      min_of_d_and_d = d2
    endif

  end function min_of_d_and_d
!-------------------------------------------------------------------------------
  pure elemental function min_of_d_and_r(d,r)

   real(my_precision)            , intent(in) :: r
   type(derivative_data_type_df5), intent(in) :: d
   type(derivative_data_type_df5)             :: min_of_d_and_r

    if ( r < d%f ) then
     min_of_d_and_r%f  = r
     min_of_d_and_r%df = zero
    else
     min_of_d_and_r = d
    endif

  end function min_of_d_and_r
!-------------------------------------------------------------------------------
  pure elemental function min_of_r_and_d(r,d2)

   real(my_precision)            , intent(in) :: r
   type(derivative_data_type_df5), intent(in) :: d2
   type(derivative_data_type_df5)             :: min_of_r_and_d

    if ( r < d2%f ) then
     min_of_r_and_d%f  = r
     min_of_r_and_d%df = zero
    else
      min_of_r_and_d = d2
    endif

  end function min_of_r_and_d

!*******************************************************************************
! Less than (Logical)
!*******************************************************************************
  pure elemental function d_less_than_d(d1, d2)

   type(derivative_data_type_df5), intent(in) :: d1, d2
   logical                                    :: d_less_than_d

    d_less_than_d = ( d1%f < d2%f )

  end function d_less_than_d
!-------------------------------------------------------------------------------
  pure elemental function d_less_than_r(d, r)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)            , intent(in) :: r
   logical                                    :: d_less_than_r

    d_less_than_r = ( d%f < r )

  end function d_less_than_r
!-------------------------------------------------------------------------------
  pure elemental function r_less_than_d(r, d)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)                  , intent(in) :: r
   logical                                :: r_less_than_d

    r_less_than_d = ( r < d%f )

  end function r_less_than_d

!*******************************************************************************
! Less than or equal (Logical)
!*******************************************************************************
  pure elemental function d_less_than_equal_d(d1, d2)

   type(derivative_data_type_df5), intent(in) :: d1, d2
   logical                                    :: d_less_than_equal_d

    d_less_than_equal_d = ( d1%f <= d2%f )

  end function d_less_than_equal_d
!-------------------------------------------------------------------------------
  pure elemental function d_less_than_equal_r(d, r)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)            , intent(in) :: r
   logical                                    :: d_less_than_equal_r

    d_less_than_equal_r = ( d%f <= r )

  end function d_less_than_equal_r
!-------------------------------------------------------------------------------
  pure elemental function r_less_than_equal_d(r, d)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)            , intent(in) :: r
   logical                                    :: r_less_than_equal_d

    r_less_than_equal_d = ( r <= d%f )

  end function r_less_than_equal_d

!*******************************************************************************
! Greater than (Logical)
!*******************************************************************************
  pure elemental function d_greater_than_d(d1, d2)

   type(derivative_data_type_df5), intent(in) :: d1, d2
   logical                                    :: d_greater_than_d

    d_greater_than_d = ( d1%f > d2%f )

  end function d_greater_than_d
!-------------------------------------------------------------------------------
  pure elemental function d_greater_than_r(d, r)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)        , intent(in) :: r
   logical                                :: d_greater_than_r

    d_greater_than_r = ( d%f > r )

  end function d_greater_than_r
!-------------------------------------------------------------------------------
  pure elemental function r_greater_than_d(r, d)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)            , intent(in) :: r
   logical                                    :: r_greater_than_d

    r_greater_than_d = ( r > d%f )

  end function r_greater_than_d

!*******************************************************************************
! Greater than or equal (Logical)
!*******************************************************************************
  pure elemental function d_greater_than_equal_d(d1, d2)

   type(derivative_data_type_df5), intent(in) :: d1, d2
   logical                                    :: d_greater_than_equal_d

    d_greater_than_equal_d = ( d1%f >= d2%f )

  end function d_greater_than_equal_d
!-------------------------------------------------------------------------------
  pure elemental function d_greater_than_equal_r(d, r)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)        , intent(in) :: r
   logical                                :: d_greater_than_equal_r

    d_greater_than_equal_r = ( d%f >= r )

  end function d_greater_than_equal_r
!-------------------------------------------------------------------------------
  pure elemental function r_greater_than_equal_d(r, d)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)            , intent(in) :: r
   logical                                    :: r_greater_than_equal_d

    r_greater_than_equal_d = ( r >= d%f )

  end function r_greater_than_equal_d

!*******************************************************************************
! Equal (Logical). This is tricky. Good to avoid using it.
!*******************************************************************************
  pure elemental function d_equal_d(d1, d2)

   type(derivative_data_type_df5), intent(in) :: d1, d2
   logical                                    :: d_equal_d

  d_equal_d = ( abs(d1%f-d2%f)/(max(abs(d1%f),abs(d2%f))+epsilon(1.0)) < epsilon(1.0) )

  end function d_equal_d
!-------------------------------------------------------------------------------
  pure elemental function d_equal_r(d, r)

   type(derivative_data_type_df5), intent(in) :: d
   real(my_precision)            , intent(in) :: r
   logical                                    :: d_equal_r

  d_equal_r = ( abs(d%f-r)/(max(abs(d%f),abs(r))+epsilon(1.0)) < epsilon(1.0) )

  end function d_equal_r
!-------------------------------------------------------------------------------
  pure elemental function r_equal_d(r, d2)

   type(derivative_data_type_df5), intent(in) :: d2
   real(my_precision)            , intent(in) :: r
   logical                                    :: r_equal_d

  r_equal_d = ( abs(d2%f-r)/(max(abs(d2%f),abs(r))+epsilon(1.0)) < epsilon(1.0) )

  end function r_equal_d

 end module derivative_data_df5

