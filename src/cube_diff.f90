!> cube_diff -- derivatives for CubeFile scalar fields
!>
!> Public API
!> ----------
!>   cube_diff_axis(cb, iax)       -> CubeFile   d f/d s_iax  (Bohr^-1)
!>   cube_diff2_axis(cb, iax)      -> CubeFile   d2f/d s_iax2 (Bohr^-2)
!>   cube_gradient(cb, gx, gy, gz)              Cartesian gradient components
!>   cube_grad_norm(cb)            -> CubeFile   |grad f|     (Bohr^-1)
!>   cube_laplacian(cb)            -> CubeFile   nabla^2 f    (Bohr^-2)
!>
!> Method dispatch (per axis, based on cb%periodicity(iax))
!> --------------------------------------------------------
!>   PERIODIC_ALL  : spectral derivative via FFT (cube_fft)
!>                   -- exact for band-limited functions; machine-precision
!>                      accuracy for smooth periodic fields.
!>   PERIODIC_NONE : 6th-order centred finite differences with zero-padding
!>                   -- accurate at interior points; O(h^6) truncation error.
!>
!> All higher-level functions (cube_gradient, cube_laplacian, etc.) call
!> cube_diff_axis / cube_diff2_axis and therefore automatically use the
!> appropriate method for each axis without any additional code.
!>
!> Non-orthogonal grids
!> --------------------
!>   Axis derivatives are in grid coordinates; cube_gradient and
!>   cube_laplacian transform to Cartesian via the unit vectors
!>   u_hat_i = axes(:,i)/|axes(:,i)|.
!>   Laplacian cross-terms (non-orthogonal grids only) use 2nd-order FD.
module cube_diff
  use cube_kinds, only: dp, ip, zero, two_pi, stencil_d1, stencil_d2
  use cube_data,  only: CubeFile, cube_copy_header, cube_is_alloc, &
                        cube_step, cube_free, PERIODIC_ALL, PERIODIC_NONE
  use cube_fft,   only: fft1d, FFT_FORWARD, FFT_BACKWARD
  implicit none
  private

  public :: cube_diff_axis    ! d f/d s_iax          (Bohr^-1)
  public :: cube_diff2_axis   ! d2f/d s_iax^2        (Bohr^-2)
  public :: cube_gradient     ! Cartesian (gx, gy, gz)
  public :: cube_grad_norm    ! |grad f|              (Bohr^-1)
  public :: cube_laplacian    ! nabla^2 f             (Bohr^-2)

contains

  ! ================================================================== !
  !  Private: FD helpers (pure -- callable from inner loops)
  ! ================================================================== !

  !> Field value at (ix, iy, iz) with per-axis boundary conditions.
  !>   Periodic  -> index wraps modulo N.
  !>   Aperiodic -> returns 0 for out-of-range index (zero-padding).
  pure real(dp) function gv(f, ix, iy, iz, nx, ny, nz, per)
    real(dp),    intent(in) :: f(:,:,:)
    integer(ip), intent(in) :: ix, iy, iz, nx, ny, nz, per(3)
    integer(ip)             :: jx, jy, jz

    if (per(1) == PERIODIC_ALL) then
      jx = modulo(ix - 1, nx) + 1
    else
      if (ix < 1 .or. ix > nx) then; gv = zero; return; end if
      jx = ix
    end if
    if (per(2) == PERIODIC_ALL) then
      jy = modulo(iy - 1, ny) + 1
    else
      if (iy < 1 .or. iy > ny) then; gv = zero; return; end if
      jy = iy
    end if
    if (per(3) == PERIODIC_ALL) then
      jz = modulo(iz - 1, nz) + 1
    else
      if (iz < 1 .or. iz > nz) then; gv = zero; return; end if
      jz = iz
    end if

    gv = f(jx, jy, jz)
  end function gv

  ! ------------------------------------------------------------------ !

  !> 6th-order centred first derivative at (ix,iy,iz) along axis iax.
  !> Returns d f / d s_iax  in Bohr^-1.
  pure real(dp) function fd1(f, ix, iy, iz, nx, ny, nz, per, iax, h)
    real(dp),    intent(in) :: f(:,:,:), h
    integer(ip), intent(in) :: ix, iy, iz, nx, ny, nz, per(3), iax
    integer(ip) :: d(3), j
    real(dp)    :: s

    d = 0;  d(iax) = 1
    s = zero
    do j = 1, 3
      s = s + stencil_d1(j) * ( &
        gv(f, ix+j*d(1), iy+j*d(2), iz+j*d(3), nx, ny, nz, per) - &
        gv(f, ix-j*d(1), iy-j*d(2), iz-j*d(3), nx, ny, nz, per) )
    end do
    fd1 = s / h
  end function fd1

  ! ------------------------------------------------------------------ !

  !> 6th-order centred second derivative at (ix,iy,iz) along axis iax.
  !> Returns d2f / d s_iax^2  in Bohr^-2.
  pure real(dp) function fd2(f, ix, iy, iz, nx, ny, nz, per, iax, h)
    real(dp),    intent(in) :: f(:,:,:), h
    integer(ip), intent(in) :: ix, iy, iz, nx, ny, nz, per(3), iax
    integer(ip) :: d(3), j
    real(dp)    :: s

    d = 0;  d(iax) = 1
    s = stencil_d2(0) * f(ix, iy, iz)
    do j = 1, 3
      s = s + stencil_d2(j) * ( &
        gv(f, ix+j*d(1), iy+j*d(2), iz+j*d(3), nx, ny, nz, per) + &
        gv(f, ix-j*d(1), iy-j*d(2), iz-j*d(3), nx, ny, nz, per) )
    end do
    fd2 = s / (h * h)
  end function fd2

  ! ------------------------------------------------------------------ !

  !> 2nd-order centred mixed partial d2f/(d s_iax d s_jax)  (Bohr^-2).
  !> Uses the 4-point cross stencil; non-zero only for non-orthogonal grids.
  pure real(dp) function fd_mixed(f, ix, iy, iz, nx, ny, nz, per, iax, jax, hi, hj)
    real(dp),    intent(in) :: f(:,:,:), hi, hj
    integer(ip), intent(in) :: ix, iy, iz, nx, ny, nz, per(3), iax, jax
    integer(ip) :: di(3), dj(3)

    di = 0;  di(iax) = 1
    dj = 0;  dj(jax) = 1

    fd_mixed = ( &
      gv(f, ix+di(1)+dj(1), iy+di(2)+dj(2), iz+di(3)+dj(3), nx, ny, nz, per) &
    - gv(f, ix+di(1)-dj(1), iy+di(2)-dj(2), iz+di(3)-dj(3), nx, ny, nz, per) &
    - gv(f, ix-di(1)+dj(1), iy-di(2)+dj(2), iz-di(3)+dj(3), nx, ny, nz, per) &
    + gv(f, ix-di(1)-dj(1), iy-di(2)-dj(2), iz-di(3)-dj(3), nx, ny, nz, per) &
    ) / (4.0_dp * hi * hj)
  end function fd_mixed

  ! ================================================================== !
  !  Private: FD loops over all voxels
  ! ================================================================== !

  subroutine fd_d1_all(cb, iax, out)
    type(CubeFile), intent(in)  :: cb
    integer(ip),    intent(in)  :: iax
    real(dp),       intent(out) :: out(:,:,:)
    integer(ip) :: ix, iy, iz
    real(dp)    :: h

    h = cube_step(cb, iax)
    do ix = 1, cb%nx
      do iy = 1, cb%ny
        do iz = 1, cb%nz
          out(ix,iy,iz) = fd1(cb%field, ix, iy, iz, &
            cb%nx, cb%ny, cb%nz, cb%periodicity, iax, h)
        end do
      end do
    end do
  end subroutine fd_d1_all

  ! ------------------------------------------------------------------ !

  subroutine fd_d2_all(cb, iax, out)
    type(CubeFile), intent(in)  :: cb
    integer(ip),    intent(in)  :: iax
    real(dp),       intent(out) :: out(:,:,:)
    integer(ip) :: ix, iy, iz
    real(dp)    :: h

    h = cube_step(cb, iax)
    do ix = 1, cb%nx
      do iy = 1, cb%ny
        do iz = 1, cb%nz
          out(ix,iy,iz) = fd2(cb%field, ix, iy, iz, &
            cb%nx, cb%ny, cb%nz, cb%periodicity, iax, h)
        end do
      end do
    end do
  end subroutine fd_d2_all

  ! ================================================================== !
  !  Private: spectral multiplier (operates on FFT coefficients)
  ! ================================================================== !

  !> Multiply FFT coefficients x(0:n-1) by the spectral operator for a
  !> derivative of order `order` (1 or 2) along an axis of step h.
  !>
  !>   order = 1 : multiply by  i * k_phys      (first derivative)
  !>   order = 2 : multiply by -k_phys^2         (second derivative)
  !>
  !> where k_phys = 2*pi * freq(k) / (N * h) and freq(k) is the
  !> centred DFT frequency.
  !>
  !> Nyquist (k = N/2, even N): zeroed for order=1 to enforce real result;
  !> retained (as -k^2) for order=2.
  pure subroutine apply_spectral_mult(x, n, h, order)
    integer(ip), intent(in)    :: n, order
    real(dp),    intent(in)    :: h
    complex(dp), intent(inout) :: x(0:n-1)

    integer(ip) :: k
    real(dp)    :: freq, kphys

    do k = 0, n - 1
      ! Centred frequency: 0, 1, ..., N/2-1, [N/2,] -(N/2-1), ..., -1
      if (k <= (n - 1) / 2) then
        freq = real(k, dp)
      else if (order == 1 .and. 2 * k == n) then   ! Nyquist for even N
        freq = 0.0_dp                               ! zero out: avoids imaginary drift
      else
        freq = real(k - n, dp)
      end if

      kphys = two_pi * freq / (real(n, dp) * h)

      if (order == 1) then
        x(k) = x(k) * cmplx(0.0_dp, kphys, kind=dp)     ! x i k
      else
        x(k) = x(k) * cmplx(-kphys**2, 0.0_dp, kind=dp) ! x (-k^2)
      end if
    end do
  end subroutine apply_spectral_mult

  ! ================================================================== !
  !  Private: spectral derivative along one axis (FFT-based)
  ! ================================================================== !

  !> Compute d^order f / d s_iax^order via FFT along axis iax.
  !>
  !> Processes N_other pencils of length N_iax in sequence; each pencil:
  !>   1. Extract complex pencil from field slice.
  !>   2. Forward FFT.
  !>   3. Multiply by spectral operator (apply_spectral_mult).
  !>   4. Backward FFT + normalise by N.
  !>   5. Store real part back.
  subroutine spectral_d(cb, iax, order, out)
    type(CubeFile), intent(in)  :: cb
    integer(ip),    intent(in)  :: iax, order
    real(dp),       intent(out) :: out(:,:,:)

    integer(ip)              :: N_ax, ix, iy, iz
    real(dp)                 :: h
    complex(dp), allocatable :: pen(:)

    select case (iax)
    case (1);  N_ax = cb%nx
    case (2);  N_ax = cb%ny
    case (3);  N_ax = cb%nz
    case default; N_ax = 0   ! iax out of range: guard against gfortran -Wmaybe-uninitialized
    end select
    h = cube_step(cb, iax)

    allocate(pen(0:N_ax - 1))

    select case (iax)
    case (1)   ! field(:,iy,iz) is contiguous in Fortran (ix fastest)
      do iz = 1, cb%nz
        do iy = 1, cb%ny
          pen = cmplx(cb%field(:, iy, iz), zero, kind=dp)
          call fft1d(pen, N_ax, FFT_FORWARD)
          call apply_spectral_mult(pen, N_ax, h, order)
          call fft1d(pen, N_ax, FFT_BACKWARD)
          out(:, iy, iz) = real(pen, kind=dp) / N_ax
        end do
      end do

    case (2)   ! field(ix,:,iz) -- stride nx
      do iz = 1, cb%nz
        do ix = 1, cb%nx
          pen = cmplx(cb%field(ix, :, iz), zero, kind=dp)
          call fft1d(pen, N_ax, FFT_FORWARD)
          call apply_spectral_mult(pen, N_ax, h, order)
          call fft1d(pen, N_ax, FFT_BACKWARD)
          out(ix, :, iz) = real(pen, kind=dp) / N_ax
        end do
      end do

    case (3)   ! field(ix,iy,:) -- stride nx*ny
      do iy = 1, cb%ny
        do ix = 1, cb%nx
          pen = cmplx(cb%field(ix, iy, :), zero, kind=dp)
          call fft1d(pen, N_ax, FFT_FORWARD)
          call apply_spectral_mult(pen, N_ax, h, order)
          call fft1d(pen, N_ax, FFT_BACKWARD)
          out(ix, iy, :) = real(pen, kind=dp) / N_ax
        end do
      end do
    end select

    deallocate(pen)
  end subroutine spectral_d

  ! ================================================================== !
  !  Public: axis derivatives
  ! ================================================================== !

  !> First derivative along grid axis iax (1, 2, or 3).
  !>
  !> Returns field = d f/d s_iax in Bohr^-1, where s_iax is arc length
  !> along the iax-th grid axis.  For orthogonal grids: axis 1->d/dx, etc.
  !>
  !>   PERIODIC_ALL  : FFT spectral derivative (machine-precision accuracy)
  !>   PERIODIC_NONE : 6th-order centred FD with zero-padding BC
  function cube_diff_axis(cb, iax) result(out)
    type(CubeFile), intent(in)  :: cb
    integer(ip),    intent(in)  :: iax
    type(CubeFile)              :: out

    if (.not. cube_is_alloc(cb)) &
      error stop "cube_diff_axis: field not allocated"
    if (iax < 1 .or. iax > 3) &
      error stop "cube_diff_axis: iax must be 1, 2, or 3"

    call cube_copy_header(cb, out)
    select case (iax)
    case (1);  out%title(1) = "df/ds_1 (axis 1)"
    case (2);  out%title(1) = "df/ds_2 (axis 2)"
    case (3);  out%title(1) = "df/ds_3 (axis 3)"
    end select

    if (cb%periodicity(iax) == PERIODIC_ALL) then
      call spectral_d(cb, iax, 1, out%field)
      out%title(2) = "spectral FFT (periodic BC)"
    else
      call fd_d1_all(cb, iax, out%field)
      out%title(2) = "6th-order centred FD (zero-padding BC)"
    end if
  end function cube_diff_axis

  ! ------------------------------------------------------------------ !

  !> Second derivative along grid axis iax.
  !> Returns d2f/ds_iax^2 in Bohr^-2.
  !>
  !>   PERIODIC_ALL  : FFT spectral (multiply by -k^2)
  !>   PERIODIC_NONE : 6th-order centred FD
  function cube_diff2_axis(cb, iax) result(out)
    type(CubeFile), intent(in)  :: cb
    integer(ip),    intent(in)  :: iax
    type(CubeFile)              :: out

    if (.not. cube_is_alloc(cb)) &
      error stop "cube_diff2_axis: field not allocated"
    if (iax < 1 .or. iax > 3) &
      error stop "cube_diff2_axis: iax must be 1, 2, or 3"

    call cube_copy_header(cb, out)
    select case (iax)
    case (1);  out%title(1) = "d2f/ds_1^2 (axis 1)"
    case (2);  out%title(1) = "d2f/ds_2^2 (axis 2)"
    case (3);  out%title(1) = "d2f/ds_3^2 (axis 3)"
    end select

    if (cb%periodicity(iax) == PERIODIC_ALL) then
      call spectral_d(cb, iax, 2, out%field)
      out%title(2) = "spectral FFT (periodic BC)"
    else
      call fd_d2_all(cb, iax, out%field)
      out%title(2) = "6th-order centred FD (zero-padding BC)"
    end if
  end function cube_diff2_axis

  ! ================================================================== !
  !  Public: Cartesian quantities  (unchanged -- use updated axis funcs)
  ! ================================================================== !

  !> Cartesian gradient components (gx, gy, gz).
  !>
  !> Uses cube_diff_axis for each axis (spectral or FD per periodicity).
  !> For non-orthogonal grids: grad f = sum_i (df/ds_i) * u_hat_i
  subroutine cube_gradient(cb, gx, gy, gz)
    type(CubeFile), intent(in)  :: cb
    type(CubeFile), intent(out) :: gx, gy, gz

    type(CubeFile) :: g1, g2, g3
    real(dp)       :: h1, h2, h3, u1(3), u2(3), u3(3)

    if (.not. cube_is_alloc(cb)) &
      error stop "cube_gradient: field not allocated"

    h1 = cube_step(cb, 1);  u1 = cb%axes(:, 1) / h1
    h2 = cube_step(cb, 2);  u2 = cb%axes(:, 2) / h2
    h3 = cube_step(cb, 3);  u3 = cb%axes(:, 3) / h3

    g1 = cube_diff_axis(cb, 1)
    g2 = cube_diff_axis(cb, 2)
    g3 = cube_diff_axis(cb, 3)

    call cube_copy_header(cb, gx)
    call cube_copy_header(cb, gy)
    call cube_copy_header(cb, gz)
    gx%title(1) = "gradient x (Cartesian)";  gx%title(2) = ""
    gy%title(1) = "gradient y (Cartesian)";  gy%title(2) = ""
    gz%title(1) = "gradient z (Cartesian)";  gz%title(2) = ""

    ! grad f = sum_i g_i * u_hat_i
    gx%field = u1(1)*g1%field + u2(1)*g2%field + u3(1)*g3%field
    gy%field = u1(2)*g1%field + u2(2)*g2%field + u3(2)*g3%field
    gz%field = u1(3)*g1%field + u2(3)*g2%field + u3(3)*g3%field

    call cube_free(g1);  call cube_free(g2);  call cube_free(g3)
  end subroutine cube_gradient

  ! ------------------------------------------------------------------ !

  !> |grad f| in Bohr^-1.
  function cube_grad_norm(cb) result(out)
    type(CubeFile), intent(in) :: cb
    type(CubeFile)             :: out, gx, gy, gz

    call cube_gradient(cb, gx, gy, gz)
    call cube_copy_header(cb, out)
    out%title(1) = "|grad f|  (Bohr^-1)"
    out%title(2) = ""
    out%field = sqrt(gx%field**2 + gy%field**2 + gz%field**2)

    call cube_free(gx);  call cube_free(gy);  call cube_free(gz)
  end function cube_grad_norm

  ! ------------------------------------------------------------------ !

  !> Laplacian: nabla^2 f = sum_k d2f/ds_k^2
  !>                      + 2 sum_{k<l} (u_hat_k.u_hat_l) d2f/ds_k ds_l
  !>
  !> Diagonal terms: cube_diff2_axis (spectral or 6th-order FD per axis).
  !> Cross-terms   : 2nd-order 4-point FD (non-zero for non-orthogonal grids
  !>                 only; zero for all standard orthogonal cube files).
  function cube_laplacian(cb) result(out)
    type(CubeFile), intent(in) :: cb
    type(CubeFile)             :: out

    real(dp)    :: h(3), udot(3,3), cosangle
    integer(ip) :: k, l, ix, iy, iz

    if (.not. cube_is_alloc(cb)) &
      error stop "cube_laplacian: field not allocated"

    do k = 1, 3
      h(k) = cube_step(cb, k)
    end do
    do k = 1, 3
      do l = 1, 3
        udot(k,l) = dot_product(cb%axes(:,k), cb%axes(:,l)) / (h(k) * h(l))
      end do
    end do

    call cube_copy_header(cb, out)
    out%title(1) = "Laplacian nabla^2 f  (Bohr^-2)"
    out%title(2) = "spectral FFT (periodic axes) / 6th-order FD (aperiodic)"
    out%field    = zero

    ! -- Diagonal: sum_k d2f/ds_k^2 (method chosen per axis) --
    do k = 1, 3
      block
        type(CubeFile) :: d2k
        d2k = cube_diff2_axis(cb, k)
        out%field = out%field + d2k%field
        call cube_free(d2k)
      end block
    end do

    ! -- Cross-terms (non-zero only for non-orthogonal grids) --
    do k = 1, 3
      do l = k + 1, 3
        cosangle = udot(k, l)
        if (abs(cosangle) < 1.0e-12_dp) cycle
        do ix = 1, cb%nx
          do iy = 1, cb%ny
            do iz = 1, cb%nz
              out%field(ix,iy,iz) = out%field(ix,iy,iz) &
                + 2.0_dp * cosangle &
                * fd_mixed(cb%field, ix, iy, iz, &
                           cb%nx, cb%ny, cb%nz, cb%periodicity, &
                           k, l, h(k), h(l))
            end do
          end do
        end do
      end do
    end do
  end function cube_laplacian

end module cube_diff
