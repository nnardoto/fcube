!> fcube_example — demonstrates fcube library with a synthetic H₂ system
!>
!> We build two Gaussian 1s electron densities centred on the two hydrogen
!> atoms of H₂ (equilibrium bond length 1.4 Bohr) and exercise the full
!> arithmetic API:
!>
!>   ρ     = ρ_A + ρ_B           total density    → integrates to ≈ 2 e⁻
!>   σ     = ρ_A − ρ_B           difference density
!>   ρ²    = ρ   * ρ             element-wise square (∝ pair density)
!>   |σ|   = cube_abs(σ)         absolute difference
!>   √ρ    = cube_sqrt(ρ)        square root of density
!>
!> After the arithmetic the cubes are written to /tmp and ρ is read back
!> to verify the round-trip fidelity of cube_io.
!>
!> Grid:  40×40×40 voxels over a 20×20×20 Bohr³ box (h = 0.5 Bohr)
!> Basis: STO-1G exponent α = 0.4166 Bohr⁻² for hydrogen 1s
program fcube_example
  use cube_kinds
  use cube_data
  use cube_arith
  use cube_io
  use cube_diff
  implicit none

  ! ── Grid parameters ──────────────────────────────────────────────── !
  integer(ip), parameter :: NX = 40, NY = 40, NZ = 40
  real(dp),    parameter :: LBOX   = 10.0_dp                   ! half-extent (Bohr)
  real(dp),    parameter :: H_STEP = 2.0_dp * LBOX / NZ        ! voxel edge (Bohr)

  ! ── Gaussian exponent (STO-1G hydrogen) ──────────────────────────── !
  real(dp), parameter :: ALPHA = 0.4166_dp    ! Bohr⁻²

  ! ── H₂ geometry (R_e = 1.4 Bohr, bond along z) ──────────────────── !
  real(dp), parameter :: RE     = 1.4_dp
  real(dp), parameter :: RA(3) = [ 0.0_dp, 0.0_dp, -RE / 2.0_dp ]
  real(dp), parameter :: RB(3) = [ 0.0_dp, 0.0_dp,  RE / 2.0_dp ]

  ! ── Working cubes ────────────────────────────────────────────────── !
  type(CubeFile) :: rhoA, rhoB          ! atom-centred densities
  type(CubeFile) :: rho                 ! total density  (rhoA + rhoB)
  type(CubeFile) :: sigma               ! difference     (rhoA - rhoB)
  type(CubeFile) :: rho_sq              ! element-wise ρ²
  type(CubeFile) :: abs_sigma           ! |σ|
  type(CubeFile) :: sqrt_rho            ! √ρ

  real(dp) :: dvol
  real(dp) :: int_rho, int_sigma, int_abs_sigma

  ! ── Banner ───────────────────────────────────────────────────────── !
  print '(a)', ""
  print '(a)', "+--------------------------------------------------+"
  print '(a)', "|      fcube -- H2 electron density example        |"
  print '(a)', "+--------------------------------------------------+"
  print '(a)', ""
  print '(a,i0,a,i0,a,i0,a)', "  Grid   : ", NX, " × ", NY, " × ", NZ, " voxels"
  print '(a,f6.3,a)',          "  Step   : ", H_STEP, " Bohr"
  print '(a,f6.3,a,f6.3,a)', &
    "  Box    : ±", LBOX, " Bohr  (", 2.0_dp*LBOX, " Bohr per axis)"
  print '(a,f6.4,a)',          "  α      : ", ALPHA, " Bohr⁻² (STO-1G H)"
  print '(a,3f7.3,a,3f7.3,a)', &
    "  H_A   : (", RA, ")  H_B: (", RB, ") Bohr"
  print '(a)', ""

  ! ──────────────────────────────────────────────────────────────────── !
  !  Step 1 — build atom-centred Gaussian densities
  ! ──────────────────────────────────────────────────────────────────── !
  print '(a)', "  [1] Building atom-centred Gaussian densities …"
  call build_gaussian(rhoA, "H2: density on atom A", RA)
  call build_gaussian(rhoB, "H2: density on atom B", RB)
  print '(a,es11.4)', "      rhoA peak = ", cube_max_val(rhoA)
  print '(a,es11.4)', "      rhoB peak = ", cube_max_val(rhoB)
  print '(a)', ""

  ! ──────────────────────────────────────────────────────────────────── !
  !  Step 2 — arithmetic operations via overloaded operators
  ! ──────────────────────────────────────────────────────────────────── !
  print '(a)', "  [2] Arithmetic operations …"
  rho       = rhoA + rhoB       ! total density
  sigma     = rhoA - rhoB       ! difference density (A minus B)
  rho_sq    = rho  * rho        ! element-wise square
  abs_sigma = cube_abs(sigma)   ! |σ|
  sqrt_rho  = cube_sqrt(rho)    ! √ρ  (element-wise)
  print '(a)', "      ρ = ρ_A + ρ_B   ✓"
  print '(a)', "      σ = ρ_A − ρ_B   ✓"
  print '(a)', "      ρ² = ρ * ρ      ✓"
  print '(a)', "      |σ|             ✓"
  print '(a)', "      √ρ              ✓"
  print '(a)', ""

  ! ──────────────────────────────────────────────────────────────────── !
  !  Step 3 — scalar operations
  ! ──────────────────────────────────────────────────────────────────── !
  print '(a)', "  [3] Scalar operations …"
  block
    type(CubeFile) :: rho_half, rho_shifted, two_rho

    rho_half    = 0.5_dp * rho           ! scale by scalar
    rho_shifted = rho + 1.0e-6_dp        ! shift by constant
    two_rho     = rho * 2.0_dp           ! multiply by scalar (rhs)

    print '(a,es11.4,a,es11.4)', &
      "      0.5*ρ   : min=", cube_min_val(rho_half),    &
      "  max=", cube_max_val(rho_half)
    print '(a,es11.4,a,es11.4)', &
      "      ρ+1e-6  : min=", cube_min_val(rho_shifted), &
      "  max=", cube_max_val(rho_shifted)
    print '(a,es11.4,a,es11.4)', &
      "      ρ*2     : min=", cube_min_val(two_rho),     &
      "  max=", cube_max_val(two_rho)

    call cube_free(rho_half)
    call cube_free(rho_shifted)
    call cube_free(two_rho)
  end block
  print '(a)', ""

  ! ──────────────────────────────────────────────────────────────────── !
  !  Step 4 — numerical integrals
  ! ──────────────────────────────────────────────────────────────────── !
  print '(a)', "  [4] Numerical integration (Σ f * ΔV) …"
  dvol = cube_voxel_vol(rho)

  int_rho       = sum(rho%field)       * dvol
  int_sigma     = sum(sigma%field)     * dvol
  int_abs_sigma = sum(abs_sigma%field) * dvol

  print '(a,f12.6,a)',  "      ∫ ρ   d³r = ", int_rho,       "  e⁻  (expected ≈ 2)"
  print '(a,f12.6,a)',  "      ∫ σ   d³r = ", int_sigma,      "  e⁻  (expected ≈ 0)"
  print '(a,f12.6)',    "      ∫ |σ| d³r = ", int_abs_sigma
  print '(a)', ""

  if (abs(int_rho - 2.0_dp) < 0.01_dp) then
    print '(a)', "      ✓ Electron count within 1% of 2.00"
  else
    print '(a)', "      ✗ WARNING: electron count deviates from 2.00"
  end if
  print '(a)', ""

  ! ──────────────────────────────────────────────────────────────────── !
  !  Step 5 — write output cubes
  ! ──────────────────────────────────────────────────────────────────── !
  print '(a)', "  [5] Writing .cube files to /tmp …"
  call write_cube("/tmp/fcube_rho.cube",       rho)
  call write_cube("/tmp/fcube_sigma.cube",     sigma)
  call write_cube("/tmp/fcube_rho_sq.cube",    rho_sq)
  call write_cube("/tmp/fcube_abs_sigma.cube", abs_sigma)
  call write_cube("/tmp/fcube_sqrt_rho.cube",  sqrt_rho)
  print '(a)', "      /tmp/fcube_rho.cube"
  print '(a)', "      /tmp/fcube_sigma.cube"
  print '(a)', "      /tmp/fcube_rho_sq.cube"
  print '(a)', "      /tmp/fcube_abs_sigma.cube"
  print '(a)', "      /tmp/fcube_sqrt_rho.cube"
  print '(a)', ""

  ! ──────────────────────────────────────────────────────────────────── !
  !  Step 6 — round-trip fidelity check
  ! ──────────────────────────────────────────────────────────────────── !
  print '(a)', "  [6] Round-trip check (write → read → compare) …"
  block
    type(CubeFile) :: rho_rt
    real(dp)       :: max_err

    call read_cube("/tmp/fcube_rho.cube", rho_rt, &
                   [PERIODIC_NONE, PERIODIC_NONE, PERIODIC_NONE])

    max_err = maxval(abs(rho_rt%field - rho%field))

    print '(a,es10.3)', "      max |field_rt − field| = ", max_err

    if (max_err < 1.0e-5_dp) then
      print '(a)', "      ✓ Round-trip fidelity OK (< 1e-5)"
    else
      print '(a)', "      ✗ Round-trip error exceeds threshold"
    end if

    call cube_free(rho_rt)
  end block
  print '(a)', ""
  ! ──────────────────────────────────────────────────────────────────── !
  !  Step 7 — Gradient and Laplacian of the total density
  ! ──────────────────────────────────────────────────────────────────── !
  print '(a)', "  [7] Derivatives of rho (cube_diff) ..."
  block
    type(CubeFile) :: gx, gy, gz, gnorm, lap
    real(dp) :: int_gx, int_gz, int_lap, dvol7

    call cube_gradient(rho, gx, gy, gz)
    gnorm = cube_grad_norm(rho)
    lap   = cube_laplacian(rho)

    dvol7  = cube_voxel_vol(rho)
    int_gx  = sum(gx%field)  * dvol7
    int_gz  = sum(gz%field)  * dvol7
    int_lap = sum(lap%field) * dvol7

    print '(a,es11.4)', "      |nabla rho| max  = ", cube_max_val(gnorm)
    print '(a,es11.4)', "      nabla^2 rho  min  = ", cube_min_val(lap)
    print '(a,f12.7)',  "      integral gx        = ", int_gx
    print '(a,f12.7)',  "      integral gz        = ", int_gz
    print '(a,f12.7)',  "      integral nabla^2   = ", int_lap
    print '(a)', "      (all integrals -> 0 by Gauss theorem)"
    print '(a)', ""

    call write_cube("/tmp/fcube_grad_norm.cube", gnorm)
    call write_cube("/tmp/fcube_laplacian.cube", lap)
    print '(a)', "      /tmp/fcube_grad_norm.cube"
    print '(a)', "      /tmp/fcube_laplacian.cube"
    print '(a)', ""

    call cube_free(gx);  call cube_free(gy);  call cube_free(gz)
    call cube_free(gnorm);  call cube_free(lap)
  end block

  print '(a)', "Done."
  print '(a)', ""

  ! ── Free all cubes ─────────────────────────────────────────────────── !
  call cube_free(rhoA);      call cube_free(rhoB)
  call cube_free(rho);       call cube_free(sigma)
  call cube_free(rho_sq);    call cube_free(abs_sigma)
  call cube_free(sqrt_rho)

contains

  !> Build a normalised 1s Gaussian electron density centred at position R.
  !>
  !>   ρ(r) = (α/π)^{3/2} exp(−α |r−R|²)
  !>
  !> This integrates to exactly 1 electron over all space; on our finite
  !> grid (±10 Bohr) the truncation error is < 10⁻¹⁷.
  !>
  !> Atom metadata: both H atoms are written to every cube so that
  !> visualisation tools (VESTA, VMD) can overlay the molecule geometry.
  subroutine build_gaussian(cb, title1, R)
    type(CubeFile),   intent(inout) :: cb
    character(len=*), intent(in)    :: title1
    real(dp),         intent(in)    :: R(3)

    integer(ip) :: ix, iy, iz
    real(dp)    :: x, y, z, rsq, norm_fac

    cb%title(1)    = title1
    cb%title(2)    = "alpha=0.4166 Bohr^-2 (STO-1G H); grid 40^3, h=0.5 Bohr"
    cb%nx          = NX;   cb%ny = NY;   cb%nz = NZ
    cb%natoms      = 2
    cb%origin      = [ -LBOX, -LBOX, -LBOX ]
    cb%axes        = zero
    cb%axes(1, 1)  = H_STEP
    cb%axes(2, 2)  = H_STEP
    cb%axes(3, 3)  = H_STEP
    cb%periodicity = [ PERIODIC_NONE, PERIODIC_NONE, PERIODIC_NONE ]

    call cube_alloc(cb)

    ! Atom records (both atoms in every cube for visualisation)
    cb%atoms(1)%z      = 1;   cb%atoms(1)%charge = 0.0_dp;   cb%atoms(1)%pos = RA
    cb%atoms(2)%z      = 1;   cb%atoms(2)%charge = 0.0_dp;   cb%atoms(2)%pos = RB

    ! Normalisation: ∫ (α/π)^{3/2} exp(−α r²) d³r = 1
    norm_fac = (ALPHA / pi) ** 1.5_dp

    do ix = 1, NX
      x = cb%origin(1) + (ix - 1) * H_STEP
      do iy = 1, NY
        y = cb%origin(2) + (iy - 1) * H_STEP
        do iz = 1, NZ
          z = cb%origin(3) + (iz - 1) * H_STEP
          rsq = (x - R(1))**2 + (y - R(2))**2 + (z - R(3))**2
          cb%field(ix, iy, iz) = norm_fac * exp(-ALPHA * rsq)
        end do
      end do
    end do
  end subroutine build_gaussian

end program fcube_example
