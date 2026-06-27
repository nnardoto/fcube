!> Tests for cube_diff: axis derivatives, Cartesian gradient, Laplacian.
!>
!> Test strategy
!> -------------
!> 1. Linear field f = C*x  ->  df/dx = C  (exact at ALL interior points
!>    because the 6th-order stencil reproduces polynomials of degree <= 6)
!>
!> 2. Quadratic field f = x^2 + y^2 + z^2  ->  nabla^2 f = 6  (exact)
!>
!> 3. Gradient direction: f = x  ->  gx=1, gy=gz=0 at interior points
!>
!> 4. Periodic grid with f = sin(2*pi*n/N)  ->  exact derivative via FD
!>    (periodic BC removes all boundary effects)
!>
!> 5. Gaussian density on non-periodic grid:
!>    - gradient at grid centre = 0  (by symmetry)
!>    - integral of nabla^2 rho ~ 0  (Gauss theorem, density decays at walls)
program test_diff
  use cube_kinds
  use cube_data
  use cube_diff
  implicit none

  integer :: fail = 0

  call test_linear_deriv()
  call test_quadratic_laplacian()
  call test_gradient_direction()
  call test_periodic_deriv()
  call test_gaussian_properties()
  call test_second_axis_deriv()
  call test_metadata()

  if (fail == 0) then
    print '(a)', "cube_diff: all tests passed."
  else
    print '(a,i0,a)', "cube_diff: ", fail, " test(s) FAILED."
    error stop 1
  end if

contains

  ! ================================================================== !
  !  Test helpers
  ! ================================================================== !

  !> Build an N^3 orthogonal non-periodic cube with step h, origin at 0.
  subroutine make_cube(cb, n, h)
    type(CubeFile), intent(inout) :: cb
    integer(ip),    intent(in)    :: n
    real(dp),       intent(in)    :: h

    cb%nx = n;  cb%ny = n;  cb%nz = n;  cb%natoms = 0
    cb%origin = zero
    cb%axes   = zero
    cb%axes(1,1) = h;  cb%axes(2,2) = h;  cb%axes(3,3) = h
    cb%periodicity = [PERIODIC_NONE, PERIODIC_NONE, PERIODIC_NONE]
    call cube_alloc(cb)
    cb%field = zero
  end subroutine make_cube

  !> Build an N^3 orthogonal PERIODIC cube with step h, origin at 0.
  subroutine make_periodic_cube(cb, n, h)
    type(CubeFile), intent(inout) :: cb
    integer(ip),    intent(in)    :: n
    real(dp),       intent(in)    :: h

    cb%nx = n;  cb%ny = n;  cb%nz = n;  cb%natoms = 0
    cb%origin = zero
    cb%axes   = zero
    cb%axes(1,1) = h;  cb%axes(2,2) = h;  cb%axes(3,3) = h
    cb%periodicity = [PERIODIC_ALL, PERIODIC_ALL, PERIODIC_ALL]
    call cube_alloc(cb)
    cb%field = zero
  end subroutine make_periodic_cube

  subroutine check(label, condition)
    character(*), intent(in) :: label
    logical,      intent(in) :: condition
    if (condition) then
      print '("  OK  ",a)', label
    else
      print '("  FAIL ",a)', label
      fail = fail + 1
    end if
  end subroutine check

  ! ================================================================== !
  !  Test 1: d/dx of C*x = C  (exact for 6th-order stencil)
  ! ================================================================== !
  subroutine test_linear_deriv()
    integer(ip), parameter :: N = 16
    real(dp),    parameter :: H = 0.3_dp, C = 2.5_dp

    type(CubeFile) :: cb, dc
    integer(ip)    :: ix, iy, iz
    real(dp)       :: x, err

    call make_cube(cb, N, H)

    ! f(x,y,z) = C * x
    do ix = 1, N
      x = cb%origin(1) + (ix - 1) * H
      cb%field(ix, :, :) = C * x
    end do

    dc = cube_diff_axis(cb, 1)    ! d/dx

    ! Check INTERIOR points only (stencil reaches 3 points from boundary)
    err = 0.0_dp
    do ix = 4, N-3
      do iy = 1, N
        do iz = 1, N
          err = max(err, abs(dc%field(ix,iy,iz) - C))
        end do
      end do
    end do
    call check("linear df/dx: interior error < 1e-12", err < 1.0e-12_dp)

    ! Derivative along y and z must be zero (f does not depend on y,z)
    err = maxval(abs(cube_diff_axis(cb, 2)%field(4:N-3, :, :)))
    call check("linear df/dy = 0 (interior)", err < 1.0e-12_dp)

    err = maxval(abs(cube_diff_axis(cb, 3)%field(:, :, 4:N-3)))
    call check("linear df/dz = 0 (interior)", err < 1.0e-12_dp)

    call cube_free(cb);  call cube_free(dc)
  end subroutine test_linear_deriv

  ! ================================================================== !
  !  Test 2: Laplacian of x^2 + y^2 + z^2 = 6  (exact)
  ! ================================================================== !
  subroutine test_quadratic_laplacian()
    integer(ip), parameter :: N = 14
    real(dp),    parameter :: H = 0.4_dp

    type(CubeFile) :: cb, lap
    integer(ip)    :: ix, iy, iz
    real(dp)       :: x, y, z, err

    call make_cube(cb, N, H)

    do ix = 1, N
      x = cb%origin(1) + (ix - 1) * H
      do iy = 1, N
        y = cb%origin(2) + (iy - 1) * H
        do iz = 1, N
          z = cb%origin(3) + (iz - 1) * H
          cb%field(ix,iy,iz) = x**2 + y**2 + z**2
        end do
      end do
    end do

    lap = cube_laplacian(cb)

    ! nabla^2(x^2+y^2+z^2) = 2+2+2 = 6  exactly at interior points
    err = 0.0_dp
    do ix = 4, N-3
      do iy = 4, N-3
        do iz = 4, N-3
          err = max(err, abs(lap%field(ix,iy,iz) - 6.0_dp))
        end do
      end do
    end do
    call check("laplacian of r^2: interior error < 1e-11", err < 1.0e-11_dp)

    ! Check title was set
    call check("laplacian: title set", len_trim(lap%title(1)) > 0)

    call cube_free(cb);  call cube_free(lap)
  end subroutine test_quadratic_laplacian

  ! ================================================================== !
  !  Test 3: gradient direction -- f = x  ->  gx=1, gy=gz=0
  ! ================================================================== !
  subroutine test_gradient_direction()
    integer(ip), parameter :: N = 12
    real(dp),    parameter :: H = 0.25_dp

    type(CubeFile) :: cb, gx, gy, gz
    integer(ip)    :: ix
    real(dp)       :: x, err

    call make_cube(cb, N, H)
    do ix = 1, N
      x = cb%origin(1) + (ix - 1) * H
      cb%field(ix, :, :) = x
    end do

    call cube_gradient(cb, gx, gy, gz)

    err = maxval(abs(gx%field(4:N-3, :, :) - 1.0_dp))
    call check("gradient: gx = 1 for f=x (interior)", err < 1.0e-12_dp)

    err = maxval(abs(gy%field(4:N-3, :, :)))
    call check("gradient: gy = 0 for f=x (interior)", err < 1.0e-12_dp)

    err = maxval(abs(gz%field(4:N-3, :, :)))
    call check("gradient: gz = 0 for f=x (interior)", err < 1.0e-12_dp)

    call cube_free(cb);  call cube_free(gx)
    call cube_free(gy);  call cube_free(gz)
  end subroutine test_gradient_direction

  ! ================================================================== !
  !  Test 4: periodic grid -- f = sin(2*pi*n/N)
  !  df/dx = (2*pi/(N*h)) * cos(2*pi*n/N)  -- no boundary effects
  ! ================================================================== !
  subroutine test_periodic_deriv()
    integer(ip), parameter :: N = 20
    real(dp),    parameter :: H = 0.5_dp

    type(CubeFile) :: cb, dc
    integer(ip)    :: ix, iy, iz
    real(dp)       :: theta, err, k_phys, analytic

    ! k_phys = 2*pi / (N*H) is the physical wavenumber
    k_phys = 2.0_dp * 3.14159265358979323846_dp / (real(N, dp) * H)

    call make_periodic_cube(cb, N, H)

    do ix = 1, N
      theta = 2.0_dp * 3.14159265358979323846_dp * (ix - 1) / real(N, dp)
      cb%field(ix, :, :) = sin(theta)
    end do

    dc = cube_diff_axis(cb, 1)

    ! Analytical: d sin(2*pi*(ix-1)/N) / dx = k_phys * cos(2*pi*(ix-1)/N)
    err = 0.0_dp
    do ix = 1, N
      do iy = 1, N
        do iz = 1, N
          analytic = k_phys * cos(2.0_dp * 3.14159265358979323846_dp &
                                  * (ix - 1) / real(N, dp))
          err = max(err, abs(dc%field(ix,iy,iz) - analytic))
        end do
      end do
    end do

    ! Spectral FFT derivative of sin: exact to floating-point precision.
    ! Error ~ N * epsilon_machine ~ 20 * 2.2e-16 ~ 4e-15
    call check("periodic deriv (FFT spectral): max err < 1e-12", err < 1.0e-12_dp)

    call cube_free(cb);  call cube_free(dc)
  end subroutine test_periodic_deriv

  ! ================================================================== !
  !  Test 5: Gaussian density -- symmetry and integral checks
  ! ================================================================== !
  subroutine test_gaussian_properties()
    integer(ip), parameter :: N = 30
    real(dp),    parameter :: H = 0.4_dp, ALPHA_G = 0.5_dp

    type(CubeFile) :: cb, gnorm, lap
    integer(ip)    :: ix, iy, iz
    real(dp)       :: x, y, z, r2, dvol, int_lap
    real(dp)       :: cx, cy, cz      ! grid centre coords
    integer(ip)    :: icx, icy, icz   ! grid centre indices

    call make_cube(cb, N, H)
    ! Centre of the grid
    icx = N/2 + 1;  icy = N/2 + 1;  icz = N/2 + 1
    cx  = (icx - 1) * H
    cy  = (icy - 1) * H
    cz  = (icz - 1) * H

    do ix = 1, N
      x = (ix - 1) * H
      do iy = 1, N
        y = (iy - 1) * H
        do iz = 1, N
          z = (iz - 1) * H
          r2 = (x-cx)**2 + (y-cy)**2 + (z-cz)**2
          cb%field(ix,iy,iz) = exp(-ALPHA_G * r2)
        end do
      end do
    end do

    ! |grad rho| at the centre should be zero (stationary point)
    gnorm = cube_grad_norm(cb)
    call check("Gaussian: |grad| at centre < 1e-12", &
      gnorm%field(icx, icy, icz) < 1.0e-12_dp)

    ! integral of nabla^2 rho ~ 0 (Gauss theorem, density -> 0 at walls)
    lap  = cube_laplacian(cb)
    dvol = H**3
    int_lap = sum(lap%field) * dvol
    call check("Gaussian: integral nabla^2 rho < 0.01", abs(int_lap) < 0.01_dp)

    ! nabla^2 rho at centre: analytic value = -6*alpha + 4*alpha^2*r^2 |_{r=0}
    !   = -6 * alpha  (for isotropic Gaussian exp(-alpha*r^2))
    ! Wait: d2/dx2 exp(-alpha x^2) = -2*alpha*exp(-alpha*x^2)*(1 - 2*alpha*x^2)
    ! At centre (r=0): d2/dx2 = -2*alpha, and summing 3 axes: nabla^2 = -6*alpha
    call check("Gaussian: laplacian at centre ~ -6*alpha", &
      abs(lap%field(icx,icy,icz) - (-6.0_dp*ALPHA_G)) < 1.0e-4_dp)

    call cube_free(cb);  call cube_free(gnorm);  call cube_free(lap)
  end subroutine test_gaussian_properties

  ! ================================================================== !
  !  Test 6: second axis derivative  d2/ds_k^2
  ! ================================================================== !
  subroutine test_second_axis_deriv()
    integer(ip), parameter :: N = 12
    real(dp),    parameter :: H = 0.3_dp

    type(CubeFile) :: cb, d2
    integer(ip)    :: ix, iy, iz
    real(dp)       :: z, err

    call make_cube(cb, N, H)

    ! f(x,y,z) = z^2  ->  d2f/dz^2 = 2.0 everywhere (interior)
    do iz = 1, N
      z = cb%origin(3) + (iz - 1) * H
      cb%field(:, :, iz) = z**2
    end do

    d2 = cube_diff2_axis(cb, 3)

    err = 0.0_dp
    do ix = 1, N
      do iy = 1, N
        do iz = 4, N-3
          err = max(err, abs(d2%field(ix,iy,iz) - 2.0_dp))
        end do
      end do
    end do
    call check("d2z2/dz^2 = 2 at interior (< 1e-11)", err < 1.0e-11_dp)

    call cube_free(cb);  call cube_free(d2)
  end subroutine test_second_axis_deriv

  ! ================================================================== !
  !  Test 7: result metadata is correct
  ! ================================================================== !
  subroutine test_metadata()
    integer(ip), parameter :: N = 6
    real(dp),    parameter :: H = 0.2_dp

    type(CubeFile) :: cb, dc, d2, gx, gy, gz, gn, lap

    call make_cube(cb, N, H)
    cb%field = 1.0_dp

    dc  = cube_diff_axis(cb, 1)
    d2  = cube_diff2_axis(cb, 2)
    call cube_gradient(cb, gx, gy, gz)
    gn  = cube_grad_norm(cb)
    lap = cube_laplacian(cb)

    ! All results should share the same geometry as the input
    call check("meta: diff_axis compatible",  cube_compatible(dc,  cb))
    call check("meta: diff2_axis compatible", cube_compatible(d2,  cb))
    call check("meta: gx compatible",         cube_compatible(gx,  cb))
    call check("meta: gy compatible",         cube_compatible(gy,  cb))
    call check("meta: gz compatible",         cube_compatible(gz,  cb))
    call check("meta: grad_norm compatible",  cube_compatible(gn,  cb))
    call check("meta: laplacian compatible",  cube_compatible(lap, cb))

    ! Derivative of constant = 0
    call check("const: df/ds_1 = 0", maxval(abs(dc%field))  < 1.0e-14_dp)
    call check("const: d2f/ds_2^2 = 0", maxval(abs(d2%field)) < 1.0e-14_dp)
    call check("const: |grad 1| = 0", maxval(abs(gn%field))  < 1.0e-14_dp)
    call check("const: nabla^2 1 = 0", maxval(abs(lap%field)) < 1.0e-14_dp)

    call cube_free(cb);   call cube_free(dc);  call cube_free(d2)
    call cube_free(gx);   call cube_free(gy);  call cube_free(gz)
    call cube_free(gn);   call cube_free(lap)
  end subroutine test_metadata

end program test_diff
