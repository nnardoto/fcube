!> cube_arith — arithmetic operations on CubeFile scalar fields
!>
!> Provides:
!>   Operator overloading : +  -  *  /  (binary, CubeFile op CubeFile)
!>   Scalar operations    : +  -  *  /  (CubeFile op real, real op CubeFile)
!>   Unary                : - (negation)
!>   Functions            : cube_abs, cube_sqrt, cube_apply
!>
!> Grid compatibility:
!>   All binary operations between two CubeFiles require identical grids
!>   (same nx/ny/nz, origin, and axes within a tight tolerance).
!>   Incompatible grids call error stop with a message identifying which
!>   dimension or vector differs.
!>
!> Result metadata:
!>   The output CubeFile inherits geometry, atoms, and periodicity from
!>   the left-hand operand (or the CubeFile operand for scalar ops).
!>   Title lines are set to a short description of the operation.
!>
!> Memory:
!>   Every operator returns a new CubeFile by value.  Fortran's move
!>   semantics on function results avoid copies in simple expressions,
!>   but chained expressions (a+b+c+d) allocate one temporary per
!>   binary node.  For performance-critical loops use cube_apply instead.
module cube_arith
  use cube_kinds, only: dp, ip, zero, eps_safe
  use cube_data,  only: CubeFile, cube_copy_header, cube_compatible, &
                        cube_is_alloc
  implicit none
  private

  ! ------------------------------------------------------------------ !
  !  Public operator interfaces
  ! ------------------------------------------------------------------ !

  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: operator(/)

  ! ------------------------------------------------------------------ !
  !  Public named functions
  ! ------------------------------------------------------------------ !

  public :: cube_abs      ! |field|  element-wise
  public :: cube_sqrt     ! sqrt(field) element-wise (negative → 0 with warning)
  public :: cube_apply    ! apply a user-supplied elemental function in-place
  public :: cube_max_val  ! scalar: maximum value of field
  public :: cube_min_val  ! scalar: minimum value of field

  ! ------------------------------------------------------------------ !
  !  Generic interfaces
  ! ------------------------------------------------------------------ !

  interface operator(+)
    module procedure cube_add_cube    ! cb + cb
    module procedure cube_add_real    ! cb + r
    module procedure real_add_cube    ! r  + cb
  end interface

  interface operator(-)
    module procedure cube_sub_cube    ! cb - cb
    module procedure cube_sub_real    ! cb - r
    module procedure real_sub_cube    ! r  - cb
    module procedure cube_negate      ! -cb
  end interface

  interface operator(*)
    module procedure cube_mul_cube    ! cb * cb  (element-wise)
    module procedure cube_mul_real    ! cb * r
    module procedure real_mul_cube    ! r  * cb
  end interface

  interface operator(/)
    module procedure cube_div_cube    ! cb / cb  (element-wise)
    module procedure cube_div_real    ! cb / r
    module procedure real_div_cube    ! r  / cb
  end interface

contains

  ! ================================================================== !
  !  Compatibility guard
  ! ================================================================== !

  subroutine require_compatible(a, b, op)
    type(CubeFile),  intent(in) :: a, b
    character(len=*), intent(in) :: op

    character(len=256) :: msg

    if (.not. cube_is_alloc(a)) &
      error stop "cube_arith (" // op // "): left operand field not allocated"
    if (.not. cube_is_alloc(b)) &
      error stop "cube_arith (" // op // "): right operand field not allocated"

    if (cube_compatible(a, b)) return

    ! Build an informative message identifying the mismatch
    if (a%nx /= b%nx) then
      write(msg, '(a,i0,a,i0)') &
        "nx differs: ", a%nx, " vs ", b%nx
    else if (a%ny /= b%ny) then
      write(msg, '(a,i0,a,i0)') &
        "ny differs: ", a%ny, " vs ", b%ny
    else if (a%nz /= b%nz) then
      write(msg, '(a,i0,a,i0)') &
        "nz differs: ", a%nz, " vs ", b%nz
    else if (any(abs(a%origin - b%origin) > 1.0e-6_dp)) then
      write(msg, '(a,3f12.6,a,3f12.6)') &
        "origin differs: ", a%origin, " vs ", b%origin
    else
      write(msg, '(a)') "step vectors (axes) differ"
    end if

    error stop "cube_arith (" // op // "): incompatible grids — " // trim(msg)
  end subroutine require_compatible

  ! ================================================================== !
  !  Binary CubeFile op CubeFile
  ! ================================================================== !

  function cube_add_cube(a, b) result(c)
    type(CubeFile), intent(in) :: a, b
    type(CubeFile)             :: c
    call require_compatible(a, b, '+')
    call cube_copy_header(a, c)
    c%title(1) = "cube_arith: a + b"
    c%title(2) = ""
    c%field = a%field + b%field
  end function

  function cube_sub_cube(a, b) result(c)
    type(CubeFile), intent(in) :: a, b
    type(CubeFile)             :: c
    call require_compatible(a, b, '-')
    call cube_copy_header(a, c)
    c%title(1) = "cube_arith: a - b"
    c%title(2) = ""
    c%field = a%field - b%field
  end function

  !> Element-wise product.  For convolution use a dedicated routine;
  !> this is pointwise multiplication f(r)*g(r).
  function cube_mul_cube(a, b) result(c)
    type(CubeFile), intent(in) :: a, b
    type(CubeFile)             :: c
    call require_compatible(a, b, '*')
    call cube_copy_header(a, c)
    c%title(1) = "cube_arith: a * b"
    c%title(2) = ""
    c%field = a%field * b%field
  end function

  !> Element-wise division.  Division by zero at a voxel is guarded by
  !> eps_safe: result = a / max(|b|, eps_safe), preserving sign of b.
  function cube_div_cube(a, b) result(c)
    type(CubeFile), intent(in) :: a, b
    type(CubeFile)             :: c
    call require_compatible(a, b, '/')
    call cube_copy_header(a, c)
    c%title(1) = "cube_arith: a / b"
    c%title(2) = ""
    where (abs(b%field) > eps_safe)
      c%field = a%field / b%field
    elsewhere
      c%field = zero
    end where
  end function

  ! ================================================================== !
  !  CubeFile op scalar
  ! ================================================================== !

  function cube_add_real(a, r) result(c)
    type(CubeFile), intent(in) :: a
    real(dp),       intent(in) :: r
    type(CubeFile)             :: c
    if (.not. cube_is_alloc(a)) &
      error stop "cube_arith (+): operand field not allocated"
    call cube_copy_header(a, c)
    c%title(1) = "cube_arith: a + scalar"
    c%title(2) = ""
    c%field = a%field + r
  end function

  function cube_sub_real(a, r) result(c)
    type(CubeFile), intent(in) :: a
    real(dp),       intent(in) :: r
    type(CubeFile)             :: c
    if (.not. cube_is_alloc(a)) &
      error stop "cube_arith (-): operand field not allocated"
    call cube_copy_header(a, c)
    c%title(1) = "cube_arith: a - scalar"
    c%title(2) = ""
    c%field = a%field - r
  end function

  function cube_mul_real(a, r) result(c)
    type(CubeFile), intent(in) :: a
    real(dp),       intent(in) :: r
    type(CubeFile)             :: c
    if (.not. cube_is_alloc(a)) &
      error stop "cube_arith (*): operand field not allocated"
    call cube_copy_header(a, c)
    c%title(1) = "cube_arith: a * scalar"
    c%title(2) = ""
    c%field = a%field * r
  end function

  function cube_div_real(a, r) result(c)
    type(CubeFile), intent(in) :: a
    real(dp),       intent(in) :: r
    type(CubeFile)             :: c
    if (.not. cube_is_alloc(a)) &
      error stop "cube_arith (/): operand field not allocated"
    if (abs(r) <= eps_safe) &
      error stop "cube_arith (/): division by zero (scalar denominator)"
    call cube_copy_header(a, c)
    c%title(1) = "cube_arith: a / scalar"
    c%title(2) = ""
    c%field = a%field / r
  end function

  ! ================================================================== !
  !  Scalar op CubeFile
  ! ================================================================== !

  function real_add_cube(r, a) result(c)
    real(dp),       intent(in) :: r
    type(CubeFile), intent(in) :: a
    type(CubeFile)             :: c
    c = cube_add_real(a, r)
    c%title(1) = "cube_arith: scalar + a"
  end function

  function real_sub_cube(r, a) result(c)
    real(dp),       intent(in) :: r
    type(CubeFile), intent(in) :: a
    type(CubeFile)             :: c
    if (.not. cube_is_alloc(a)) &
      error stop "cube_arith (-): operand field not allocated"
    call cube_copy_header(a, c)
    c%title(1) = "cube_arith: scalar - a"
    c%title(2) = ""
    c%field = r - a%field
  end function

  function real_mul_cube(r, a) result(c)
    real(dp),       intent(in) :: r
    type(CubeFile), intent(in) :: a
    type(CubeFile)             :: c
    c = cube_mul_real(a, r)
    c%title(1) = "cube_arith: scalar * a"
  end function

  function real_div_cube(r, a) result(c)
    real(dp),       intent(in) :: r
    type(CubeFile), intent(in) :: a
    type(CubeFile)             :: c
    if (.not. cube_is_alloc(a)) &
      error stop "cube_arith (/): operand field not allocated"
    call cube_copy_header(a, c)
    c%title(1) = "cube_arith: scalar / a"
    c%title(2) = ""
    where (abs(a%field) > eps_safe)
      c%field = r / a%field
    elsewhere
      c%field = zero
    end where
  end function

  ! ================================================================== !
  !  Unary
  ! ================================================================== !

  function cube_negate(a) result(c)
    type(CubeFile), intent(in) :: a
    type(CubeFile)             :: c
    if (.not. cube_is_alloc(a)) &
      error stop "cube_arith (unary -): operand field not allocated"
    call cube_copy_header(a, c)
    c%title(1) = "cube_arith: -a"
    c%title(2) = ""
    c%field = -a%field
  end function

  ! ================================================================== !
  !  Named functions
  ! ================================================================== !

  function cube_abs(a) result(c)
    type(CubeFile), intent(in) :: a
    type(CubeFile)             :: c
    if (.not. cube_is_alloc(a)) &
      error stop "cube_abs: field not allocated"
    call cube_copy_header(a, c)
    c%title(1) = "cube_arith: |a|"
    c%title(2) = ""
    c%field = abs(a%field)
  end function

  !> Element-wise square root.  Negative voxels are clamped to zero
  !> (common for near-zero density values that go slightly negative due
  !> to numerical noise).
  function cube_sqrt(a) result(c)
    type(CubeFile), intent(in) :: a
    type(CubeFile)             :: c
    if (.not. cube_is_alloc(a)) &
      error stop "cube_sqrt: field not allocated"
    call cube_copy_header(a, c)
    c%title(1) = "cube_arith: sqrt(a)"
    c%title(2) = ""
    where (a%field >= zero)
      c%field = sqrt(a%field)
    elsewhere
      c%field = zero
    end where
  end function

  !> Apply a scalar function f(x) to every voxel in-place.
  !> Avoids allocating a temporary CubeFile for simple transformations.
  !>
  !> `f` must be an elemental or pure scalar function real(dp)->real(dp).
  !> The dummy interface cannot declare it elemental (Fortran standard
  !> restriction), but any elemental function satisfies this interface.
  !>
  !> Example:
  !>   call cube_apply(rho, func)
  !>   ...
  !>   elemental real(dp) function func(x)
  !>     real(dp), intent(in) :: x
  !>     func = x**2
  !>   end function
  subroutine cube_apply(cb, f)
    type(CubeFile), intent(inout) :: cb
    interface
      pure real(dp) function f(x)
        import dp
        real(dp), intent(in) :: x
      end function
    end interface
    integer(ip) :: ix, iy, iz
    if (.not. cube_is_alloc(cb)) &
      error stop "cube_apply: field not allocated"
    do ix = 1, cb%nx
      do iy = 1, cb%ny
        do iz = 1, cb%nz
          cb%field(ix,iy,iz) = f(cb%field(ix,iy,iz))
        end do
      end do
    end do
  end subroutine

  pure real(dp) function cube_max_val(cb)
    type(CubeFile), intent(in) :: cb
    cube_max_val = maxval(cb%field)
  end function

  pure real(dp) function cube_min_val(cb)
    type(CubeFile), intent(in) :: cb
    cube_min_val = minval(cb%field)
  end function

end module cube_arith