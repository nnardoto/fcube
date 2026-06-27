!> cube_fft -- pure Fortran 1D complex fast Fourier transform
!>
!> Public API
!> ----------
!>   call fft1d(x, n, sign)
!>
!> where x(0:n-1) is complex(dp), n is the transform length, and sign is
!> FFT_FORWARD (-1) or FFT_BACKWARD (+1).
!>
!> Convention (identical to FFTW):
!>   Forward  (sign=-1): X(k) = sum_{n=0}^{N-1} x(n) exp(-2 pi i k n / N)
!>   Backward (sign=+1): x(n) = sum_{k=0}^{N-1} X(k) exp(+2 pi i k n / N)
!>   Normalisation: the caller is responsible; divide by N after backward.
!>
!> Algorithm selection:
!>   N = 2^m  : Cooley-Tukey radix-2 in-place FFT  O(N log N)
!>   otherwise: direct DFT                           O(N^2)
!>              (acceptable for typical axis sizes; accurate for any N)
module cube_fft
  use cube_kinds, only: dp, ip, two_pi
  implicit none
  private

  !> Sign constant for forward transform: X(k) = sum x(n) exp(-2 pi i k n/N)
  integer(ip), parameter, public :: FFT_FORWARD  = -1_ip

  !> Sign constant for backward (inverse) transform (not normalised).
  integer(ip), parameter, public :: FFT_BACKWARD = +1_ip

  public :: fft1d

contains

  ! ================================================================== !
  !  Public entry point
  ! ================================================================== !

  !> In-place 1D complex DFT.
  !>
  !> x(0:n-1) is overwritten on return.
  !> The result is NOT normalised: divide by n after a backward transform.
  subroutine fft1d(x, n, sign)
    integer(ip), intent(in)    :: n, sign
    complex(dp), intent(inout) :: x(0:n-1)

    if (n <= 1) return

    if (is_pow2(n)) then
      call fft_radix2(x, n, sign)
    else
      call dft_direct(x, n, sign)
    end if
  end subroutine fft1d

  ! ================================================================== !
  !  Internal helpers
  ! ================================================================== !

  !> True iff n is a positive power of 2.
  pure logical function is_pow2(n)
    integer(ip), intent(in) :: n
    is_pow2 = (n > 0) .and. (iand(n, n - 1) == 0)
  end function is_pow2

  ! ------------------------------------------------------------------ !

  !> Cooley-Tukey radix-2 in-place FFT (n must be a power of 2).
  !>
  !> Step 1: bit-reversal permutation (0-based index).
  !> Step 2: log2(N) butterfly stages with twiddle factor
  !>         wlen = exp(sign * 2 pi i / len).
  !> For sign = FFT_FORWARD (-1): wlen = exp(-2 pi i / len)  [standard DFT].
  subroutine fft_radix2(x, n, sign)
    integer(ip), intent(in)    :: n, sign
    complex(dp), intent(inout) :: x(0:n-1)

    integer(ip) :: i, j, k, half, len
    complex(dp) :: u, v, w, wlen
    real(dp)    :: angle

    ! -- Bit-reversal permutation (0-based) ----------------------------
    j = 0
    do i = 1, n - 1
      k = n / 2
      do while (j >= k)
        j = j - k
        k = k / 2
      end do
      j = j + k
      if (i < j) then
        u = x(i);  x(i) = x(j);  x(j) = u
      end if
    end do

    ! -- Butterfly stages ----------------------------------------------
    len = 2
    do while (len <= n)
      half  = len / 2
      angle = real(sign, dp) * two_pi / real(len, dp)
      wlen  = cmplx(cos(angle), sin(angle), kind=dp)
      k = 0
      do while (k < n)
        w = (1.0_dp, 0.0_dp)
        do i = 0, half - 1
          u = x(k + i)
          v = x(k + i + half) * w
          x(k + i)        = u + v
          x(k + i + half) = u - v
          w = w * wlen
        end do
        k = k + len
      end do
      len = len * 2
    end do
  end subroutine fft_radix2

  ! ------------------------------------------------------------------ !

  !> Direct DFT, O(N^2).  Used when N is not a power of 2.
  !> Accurate for any N; the bottleneck is N^2 complex multiplications.
  subroutine dft_direct(x, n, sign)
    integer(ip), intent(in)    :: n, sign
    complex(dp), intent(inout) :: x(0:n-1)

    complex(dp) :: y(0:n-1), tw
    integer(ip) :: k, j
    real(dp)    :: angle

    do k = 0, n - 1
      y(k) = (0.0_dp, 0.0_dp)
      do j = 0, n - 1
        angle = real(sign, dp) * two_pi * real(k * j, dp) / real(n, dp)
        tw    = cmplx(cos(angle), sin(angle), kind=dp)
        y(k)  = y(k) + x(j) * tw
      end do
    end do
    x = y
  end subroutine dft_direct

end module cube_fft
