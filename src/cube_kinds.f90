!> cube_kinds — precision kinds and physical constants
!>
!> All numeric literals in the library are written as KIND-qualified constants
!> imported from here.  Never write 1.0d0 or kind(0.0d0) elsewhere.
!>
!> Unit convention (matches the .cube file format):
!>   Lengths  : Bohr (atomic units)
!>   Energies : Hartree (atomic units)
!>   Charge   : electron charge (e)
!>
!> Conversion factors follow CODATA 2018 (NIST SP 961, 2019).
module cube_kinds
  implicit none
  private

  ! ------------------------------------------------------------------ !
  !  Floating-point kinds
  ! ------------------------------------------------------------------ !

  !> Double precision: 15 significant digits, exponent range ±307.
  !> This is the default throughout the library.
  integer, parameter, public :: dp = selected_real_kind(15, 307)

  !> Single precision: available for output or memory-constrained paths.
  integer, parameter, public :: sp = selected_real_kind(6, 37)

  !> Default integer kind (standard 4-byte on every relevant platform).
  integer, parameter, public :: ip = selected_int_kind(9)

  !> Large integer kind for flat array indices on grids > 2^31 points.
  integer, parameter, public :: lp = selected_int_kind(18)

  ! ------------------------------------------------------------------ !
  !  Mathematical constants
  ! ------------------------------------------------------------------ !

  real(dp), parameter, public :: pi    = 3.14159265358979323846_dp
  real(dp), parameter, public :: two_pi = 2.0_dp * pi
  real(dp), parameter, public :: four_pi = 4.0_dp * pi
  real(dp), parameter, public :: sqrt2 = 1.41421356237309504880_dp
  real(dp), parameter, public :: third = 1.0_dp / 3.0_dp
  real(dp), parameter, public :: zero  = 0.0_dp
  real(dp), parameter, public :: one   = 1.0_dp

  ! ------------------------------------------------------------------ !
  !  Physical constants (CODATA 2018)
  ! ------------------------------------------------------------------ !

  !> Bohr radius in Ångström
  real(dp), parameter, public :: bohr_to_ang  = 0.529177210903_dp

  !> Ångström to Bohr
  real(dp), parameter, public :: ang_to_bohr  = 1.0_dp / bohr_to_ang

  !> Hartree in eV
  real(dp), parameter, public :: hartree_to_ev = 27.211386245988_dp

  !> eV to Hartree
  real(dp), parameter, public :: ev_to_hartree = 1.0_dp / hartree_to_ev

  !> Hartree in kJ/mol
  real(dp), parameter, public :: hartree_to_kjmol = 2625.4996394799_dp

  !> Hartree in kcal/mol
  real(dp), parameter, public :: hartree_to_kcalmol = 627.50947406_dp

  !> Atomic mass unit in kg
  real(dp), parameter, public :: amu_to_kg = 1.66053906660e-27_dp

  !> Elementary charge in Coulombs
  real(dp), parameter, public :: elem_charge = 1.602176634e-19_dp

  !> Vacuum permittivity in SI (F/m)
  real(dp), parameter, public :: eps0_si = 8.8541878128e-12_dp

  ! ------------------------------------------------------------------ !
  !  Stencil coefficients (stored here so cube_diff has a single source
  !  of truth; the caller multiplies by 1/h^n for the right derivative)
  ! ------------------------------------------------------------------ !

  !> 6th-order centred first-derivative stencil coefficients
  !> f'(x) ≈ (c1*(f(x+h)-f(x-h)) + c2*(f(x+2h)-f(x-2h))
  !>         + c3*(f(x+3h)-f(x-3h))) / h
  !> Indices: stencil_d1(1)=c1, stencil_d1(2)=c2, stencil_d1(3)=c3
  real(dp), parameter, public :: stencil_d1(3) = [ &
       3.0_dp / 4.0_dp,   &  ! c1
      -3.0_dp / 20.0_dp,  &  ! c2
       1.0_dp / 60.0_dp   ]  ! c3

  !> 6th-order centred second-derivative stencil coefficients
  !> f''(x) ≈ (c0*f(x) + c1*(f(x+h)+f(x-h)) + c2*(f(x+2h)+f(x-2h))
  !>          + c3*(f(x+3h)+f(x-3h))) / h^2
  !> Indices: stencil_d2(0)=c0, stencil_d2(1)=c1, ..., stencil_d2(3)=c3
  real(dp), parameter, public :: stencil_d2(0:3) = [ &
      -49.0_dp / 18.0_dp,  &  ! c0 (centre)
       3.0_dp / 2.0_dp,    &  ! c1
      -3.0_dp / 20.0_dp,   &  ! c2
       1.0_dp / 90.0_dp    ]  ! c3

  ! ------------------------------------------------------------------ !
  !  Tiny numbers — guards against division by zero in ELF / NCI
  ! ------------------------------------------------------------------ !

  !> Density floor: below this, rho is treated as zero (a.u.)
  real(dp), parameter, public :: rho_tol = 1.0e-10_dp

  !> General numerical epsilon for ratio guards
  real(dp), parameter, public :: eps_safe = 1.0e-30_dp

end module cube_kinds