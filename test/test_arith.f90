!> Tests for cube_arith: operator overloading, scalar ops, unary,
!> named functions, compatibility guards, and metadata propagation.
program test_arith
  use cube_kinds
  use cube_data
  use cube_arith
  implicit none

  integer :: fail = 0

  call test_binary_ops()
  call test_scalar_ops()
  call test_unary()
  call test_named_funcs()
  call test_metadata()
  call test_incompatible()

  if (fail == 0) then
    print '(a)', "cube_arith: all tests passed."
  else
    print '(a,i0,a)', "cube_arith: ", fail, " test(s) FAILED."
    error stop 1
  end if

contains

  ! ------------------------------------------------------------------ !
  !  Build a simple NxNxN cube filled with a constant value
  ! ------------------------------------------------------------------ !
  subroutine make_cube(cb, n, val)
    type(CubeFile), intent(inout) :: cb
    integer(ip),    intent(in)    :: n
    real(dp),       intent(in)    :: val

    cb%nx = n;  cb%ny = n;  cb%nz = n;  cb%natoms = 0
    cb%origin = zero
    cb%axes   = zero
    cb%axes(1,1) = 0.2_dp;  cb%axes(2,2) = 0.2_dp;  cb%axes(3,3) = 0.2_dp
    call cube_alloc(cb)
    cb%field = val
  end subroutine

  ! ------------------------------------------------------------------ !
  subroutine test_binary_ops()
    type(CubeFile) :: a, b, c

    call make_cube(a, 4, 3.0_dp)
    call make_cube(b, 4, 2.0_dp)

    ! Addition
    c = a + b
    call check("add: value",    all_close(c%field, 5.0_dp))
    call check("add: shape x",  size(c%field,1) == 4)
    call check("add: shape z",  size(c%field,3) == 4)
    call cube_free(c)

    ! Subtraction
    c = a - b
    call check("sub: value",  all_close(c%field, 1.0_dp))
    call cube_free(c)

    ! Element-wise multiply
    c = a * b
    call check("mul: value",  all_close(c%field, 6.0_dp))
    call cube_free(c)

    ! Element-wise divide
    c = a / b
    call check("div: value",  all_close(c%field, 1.5_dp))
    call cube_free(c)

    ! Chained expression
    c = (a + b) * (a - b)   ! (3+2)*(3-2) = 5*1 = 5
    call check("chain: value", all_close(c%field, 5.0_dp))
    call cube_free(c)

    call cube_free(a);  call cube_free(b)
  end subroutine

  ! ------------------------------------------------------------------ !
  subroutine test_scalar_ops()
    type(CubeFile) :: a, c

    call make_cube(a, 3, 4.0_dp)

    ! cb + scalar
    c = a + 1.0_dp
    call check("cb+r: value",  all_close(c%field, 5.0_dp))
    call cube_free(c)

    ! cb - scalar
    c = a - 1.0_dp
    call check("cb-r: value",  all_close(c%field, 3.0_dp))
    call cube_free(c)

    ! cb * scalar
    c = a * 2.0_dp
    call check("cb*r: value",  all_close(c%field, 8.0_dp))
    call cube_free(c)

    ! cb / scalar
    c = a / 2.0_dp
    call check("cb/r: value",  all_close(c%field, 2.0_dp))
    call cube_free(c)

    ! scalar + cb
    c = 10.0_dp + a
    call check("r+cb: value",  all_close(c%field, 14.0_dp))
    call cube_free(c)

    ! scalar - cb
    c = 10.0_dp - a
    call check("r-cb: value",  all_close(c%field, 6.0_dp))
    call cube_free(c)

    ! scalar * cb
    c = 3.0_dp * a
    call check("r*cb: value",  all_close(c%field, 12.0_dp))
    call cube_free(c)

    ! scalar / cb (with non-zero denominator)
    c = 8.0_dp / a    ! 8/4 = 2
    call check("r/cb: value",  all_close(c%field, 2.0_dp))
    call cube_free(c)

    call cube_free(a)
  end subroutine

  ! ------------------------------------------------------------------ !
  subroutine test_unary()
    type(CubeFile) :: a, c

    call make_cube(a, 3, 5.0_dp)
    c = -a
    call check("neg: value",  all_close(c%field, -5.0_dp))
    call cube_free(c);  call cube_free(a)
  end subroutine

  ! ------------------------------------------------------------------ !
  subroutine test_named_funcs()
    type(CubeFile) :: a, c
    integer(ip) :: ix, iy, iz

    ! cube_abs
    call make_cube(a, 3, -2.5_dp)
    c = cube_abs(a)
    call check("abs: positive",  all_close(c%field, 2.5_dp))
    call cube_free(c);  call cube_free(a)

    ! cube_sqrt — normal values
    call make_cube(a, 3, 4.0_dp)
    c = cube_sqrt(a)
    call check("sqrt: value",   all_close(c%field, 2.0_dp))
    call cube_free(c)

    ! cube_sqrt — negative values clamped to zero
    a%field = -1.0_dp
    c = cube_sqrt(a)
    call check("sqrt: negative→0", all_close(c%field, 0.0_dp))
    call cube_free(c);  call cube_free(a)

    ! cube_max_val / cube_min_val
    call make_cube(a, 4, 0.0_dp)
    do ix = 1, 4; do iy = 1, 4; do iz = 1, 4
      a%field(ix,iy,iz) = real(ix + iy + iz, dp)
    end do; end do; end do
    call check("maxval", abs(cube_max_val(a) - 12.0_dp) < 1.0e-14_dp)
    call check("minval", abs(cube_min_val(a) -  3.0_dp) < 1.0e-14_dp)
    call cube_free(a)

    ! cube_apply — in-place, no allocation
    call make_cube(a, 3, 3.0_dp)
    call cube_apply(a, square)
    call check("apply: x^2",  all_close(a%field, 9.0_dp))
    call cube_free(a)

    ! Division by zero guard: scalar / cube where some voxels are zero
    call make_cube(a, 2, 0.0_dp)   ! all zeros
    c = 1.0_dp / a                  ! guarded: result should be zero everywhere
    call check("r/0: guarded zero", all_close(c%field, 0.0_dp))
    call cube_free(c);  call cube_free(a)

    ! Element-wise division guard
    call make_cube(a, 2, 0.0_dp)
    block
      type(CubeFile) :: b
      call make_cube(b, 2, 5.0_dp)
      c = b / a    ! denominator zero → result zero
      call check("cb/0: guarded zero", all_close(c%field, 0.0_dp))
      call cube_free(b)
    end block
    call cube_free(c);  call cube_free(a)
  end subroutine

  ! ------------------------------------------------------------------ !
  subroutine test_metadata()
    type(CubeFile) :: a, b, c

    call make_cube(a, 3, 1.0_dp)
    call make_cube(b, 3, 2.0_dp)
    a%periodicity = [1_ip, 1_ip, 0_ip]   ! slab

    c = a + b
    ! Geometry inherited from left operand
    call check("meta: compatible with a",  cube_compatible(c, a))
    ! Periodicity inherited from left operand
    call check("meta: periodicity x", c%periodicity(1) == 1_ip)
    call check("meta: periodicity z", c%periodicity(3) == 0_ip)
    ! natoms preserved
    call check("meta: natoms", c%natoms == a%natoms)

    call cube_free(a);  call cube_free(b);  call cube_free(c)
  end subroutine

  ! ------------------------------------------------------------------ !
  subroutine test_incompatible()
    type(CubeFile) :: a, b

    call make_cube(a, 4, 1.0_dp)

    ! Different nx
    b%nx = 5;  b%ny = 4;  b%nz = 4;  b%natoms = 0
    b%origin = zero
    b%axes   = zero
    b%axes(1,1) = 0.2_dp;  b%axes(2,2) = 0.2_dp;  b%axes(3,3) = 0.2_dp
    call cube_alloc(b)
    b%field = 1.0_dp

    call check("incompat: not compatible", .not. cube_compatible(a, b))

    ! We cannot test error stop directly (it terminates the program),
    ! but we confirm the guard function itself detects the mismatch.
    call cube_free(a);  call cube_free(b)
  end subroutine

  ! ================================================================== !
  !  Helpers
  ! ================================================================== !

  !> True if all elements of arr are within 1e-12 of val.
  logical function all_close(arr, val)
    real(dp), intent(in) :: arr(:,:,:), val
    all_close = all(abs(arr - val) < 1.0e-12_dp)
  end function

  pure real(dp) function square(x)
    real(dp), intent(in) :: x
    square = x * x
  end function

  subroutine check(label, condition)
    character(*), intent(in) :: label
    logical,      intent(in) :: condition
    if (condition) then
      print '("  OK  ",a)', label
    else
      print '("  FAIL ",a)', label
      fail = fail + 1
    end if
  end subroutine

end program test_arith