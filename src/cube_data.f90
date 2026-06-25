!> cube_data — central data type for a Gaussian .cube file
!>
!> All coordinates and lengths are in Bohr (atomic units), matching the
!> .cube file convention.  The step vectors axes(:,i) are NOT required to
!> be orthogonal; the full 3×3 matrix is stored so that non-orthogonal
!> (e.g. hexagonal) grids are supported transparently.
!>
!> Memory layout:  data(ix, iy, iz)
!>   ix is the slowest index (stride = ny*nz*8 bytes)
!>   iz is the fastest index (stride = 8 bytes)
!> This matches the .cube file scan order (z varies fastest on disk) so
!> reading and writing are sequential with no transposition needed.
!>
!> The `periodicity` field records which axes are periodic, enabling
!> cube_diff to switch automatically between the finite-difference (FD)
!> and FFT kernels.  A fully periodic grid uses FFT on all three axes;
!> a molecular (non-periodic) grid uses 6th-order FD throughout.
module cube_data
  use cube_kinds, only: dp, ip, lp, zero, one, bohr_to_ang
  implicit none
  private

  ! ------------------------------------------------------------------ !
  !  Periodicity flag constants (stored in CubeFile%periodicity)
  ! ------------------------------------------------------------------ !

  !> Axis is NOT periodic — use finite-difference stencils at borders.
  integer(ip), parameter, public :: PERIODIC_NONE = 0

  !> Axis IS periodic — use FFT-based derivatives on this axis.
  integer(ip), parameter, public :: PERIODIC_ALL  = 1

  ! ------------------------------------------------------------------ !
  !  Atom record
  ! ------------------------------------------------------------------ !

  !> Single atom: atomic number, partial charge (as written in .cube),
  !> and Cartesian position in Bohr.
  type, public :: Atom
    integer(ip) :: z       = 0       !! atomic number
    real(dp)    :: charge  = zero    !! partial charge (from header)
    real(dp)    :: pos(3)  = zero    !! position in Bohr
  end type Atom

  ! ------------------------------------------------------------------ !
  !  Main type
  ! ------------------------------------------------------------------ !

  type, public :: CubeFile

    ! -- Header --
    character(len=80) :: title(2) = ''  !! two free-text comment lines

    ! -- Grid geometry --
    integer(ip) :: nx = 0, ny = 0, nz = 0   !! number of voxels per axis
    real(dp)    :: origin(3) = zero          !! grid origin in Bohr

    !> Step vectors in Bohr: axes(:, 1) = step along x-axis, etc.
    !> For orthogonal grids axes(:,i) = h_i * e_i; for non-orthogonal
    !> grids the full 3×3 matrix is needed.
    real(dp) :: axes(3, 3) = zero

    !> Which axes are periodic: PERIODIC_NONE or PERIODIC_ALL per axis.
    !> Default: all periodic — override explicitly for molecular grids.
    !> Pass periodicity=[PERIODIC_NONE,PERIODIC_NONE,PERIODIC_NONE] to
    !> read_cube, or set the field directly after construction.
    integer(ip) :: periodicity(3) = PERIODIC_ALL

    ! -- Molecule --
    integer(ip) :: natoms = 0
    type(Atom), allocatable :: atoms(:)   !! (natoms)

    ! -- Volumetric data --
    !> Scalar field sampled on the grid: data(ix, iy, iz).
    !> iz is the fastest-varying index, matching .cube disk order.
    real(dp), allocatable :: field(:,:,:)   !! (nx, ny, nz)

    ! -- Optional: spin / multi-orbital .cube --
    !> Number of orbital/spin components beyond the first (0 = standard).
    integer(ip) :: ncomp = 0

    ! -- Provenance --
    character(len=256) :: source_file = ''  !! path this was read from, if any

  end type CubeFile

  ! ------------------------------------------------------------------ !
  !  Public procedures
  ! ------------------------------------------------------------------ !

  public :: cube_alloc       ! allocate atoms + field after header is set
  public :: cube_free        ! deallocate all allocatables
  public :: cube_is_alloc    ! check whether field is allocated
  public :: cube_compatible  ! check whether two grids share the same geometry
  public :: cube_voxel_vol   ! scalar volume of one voxel in Bohr^3
  public :: cube_total_vol   ! total grid volume in Bohr^3
  public :: cube_step        ! step length along one axis (norm of axes(:,i))
  public :: cube_copy_header ! copy geometry/atom data, leave field unset
  public :: cube_clone       ! deep copy of an entire CubeFile

contains

  ! ------------------------------------------------------------------ !

  !> Allocate `atoms(natoms)` and `field(nx,ny,nz)`.
  !> Call this after filling the scalar header fields (nx, ny, nz, natoms).
  subroutine cube_alloc(cb)
    type(CubeFile), intent(inout) :: cb

    if (cb%natoms < 0) &
      error stop "cube_alloc: natoms < 0"
    if (cb%nx < 1 .or. cb%ny < 1 .or. cb%nz < 1) &
      error stop "cube_alloc: grid dimensions must be >= 1"

    ! Reallocate atoms only if size changed (preserves existing atom data)
    if (allocated(cb%atoms)) then
      if (size(cb%atoms) /= cb%natoms) then
        deallocate(cb%atoms)
        allocate(cb%atoms(cb%natoms))
      end if
    else
      allocate(cb%atoms(cb%natoms))
    end if

    if (allocated(cb%field)) deallocate(cb%field)
    allocate(cb%field(cb%nx, cb%ny, cb%nz))
    cb%field = zero
  end subroutine cube_alloc

  ! ------------------------------------------------------------------ !

  !> Deallocate all allocatable components and reset to default state.
  subroutine cube_free(cb)
    type(CubeFile), intent(inout) :: cb

    if (allocated(cb%atoms)) deallocate(cb%atoms)
    if (allocated(cb%field)) deallocate(cb%field)

    cb%nx = 0;  cb%ny = 0;  cb%nz = 0
    cb%natoms    = 0
    cb%ncomp     = 0
    cb%origin    = zero
    cb%axes      = zero
    cb%periodicity = PERIODIC_NONE
    cb%title     = ''
    cb%source_file = ''
  end subroutine cube_free

  ! ------------------------------------------------------------------ !

  !> Returns .true. if field is allocated (i.e. cube_alloc has been called).
  pure logical function cube_is_alloc(cb)
    type(CubeFile), intent(in) :: cb
    cube_is_alloc = allocated(cb%field)
  end function cube_is_alloc

  ! ------------------------------------------------------------------ !

  !> Returns .true. if two CubeFiles have the same grid geometry
  !> (same nx/ny/nz, same origin and axes within tolerance `tol`).
  !> Does NOT compare field data.
  pure logical function cube_compatible(a, b, tol)
    type(CubeFile), intent(in)           :: a, b
    real(dp),       intent(in), optional :: tol

    real(dp) :: eps
    integer(ip) :: i, j

    eps = 1.0e-6_dp
    if (present(tol)) eps = tol

    cube_compatible = .false.

    if (a%nx /= b%nx .or. a%ny /= b%ny .or. a%nz /= b%nz) return

    do i = 1, 3
      if (abs(a%origin(i) - b%origin(i)) > eps) return
      do j = 1, 3
        if (abs(a%axes(j,i) - b%axes(j,i)) > eps) return
      end do
    end do

    cube_compatible = .true.
  end function cube_compatible

  ! ------------------------------------------------------------------ !

  !> Volume of a single voxel in Bohr^3.
  !> Computed as |det(axes)| — works for non-orthogonal grids.
  pure real(dp) function cube_voxel_vol(cb)
    type(CubeFile), intent(in) :: cb
    real(dp) :: a(3), b(3), c(3)

    a = cb%axes(:, 1)
    b = cb%axes(:, 2)
    c = cb%axes(:, 3)

    ! det via scalar triple product  a · (b × c)
    cube_voxel_vol = abs( &
      a(1)*(b(2)*c(3) - b(3)*c(2)) - &
      a(2)*(b(1)*c(3) - b(3)*c(1)) + &
      a(3)*(b(1)*c(2) - b(2)*c(1))   )
  end function cube_voxel_vol

  ! ------------------------------------------------------------------ !

  !> Total volume of the grid in Bohr^3: voxel_vol * nx * ny * nz.
  pure real(dp) function cube_total_vol(cb)
    type(CubeFile), intent(in) :: cb
    cube_total_vol = cube_voxel_vol(cb) * real(cb%nx, dp) &
                                        * real(cb%ny, dp) &
                                        * real(cb%nz, dp)
  end function cube_total_vol

  ! ------------------------------------------------------------------ !

  !> Euclidean norm of the step vector along axis `iax` (1, 2, or 3).
  !> For orthogonal grids this is simply h_x, h_y, or h_z.
  pure real(dp) function cube_step(cb, iax)
    type(CubeFile), intent(in) :: cb
    integer(ip),    intent(in) :: iax
    real(dp) :: v(3)

    v = cb%axes(:, iax)
    cube_step = sqrt(v(1)**2 + v(2)**2 + v(3)**2)
  end function cube_step

  ! ------------------------------------------------------------------ !

  !> Copy grid geometry (title, nx/ny/nz, origin, axes, periodicity,
  !> natoms, atoms) into `dst`, then allocate dst%field (zeroed).
  !> Useful for creating a result cube with the same grid as a source.
  subroutine cube_copy_header(src, dst)
    type(CubeFile), intent(in)    :: src
    type(CubeFile), intent(inout) :: dst

    dst%title       = src%title
    dst%nx          = src%nx
    dst%ny          = src%ny
    dst%nz          = src%nz
    dst%origin      = src%origin
    dst%axes        = src%axes
    dst%periodicity = src%periodicity
    dst%natoms      = src%natoms
    dst%ncomp       = src%ncomp

    ! Copy atoms before cube_alloc so alloc preserves them
    if (allocated(dst%atoms)) deallocate(dst%atoms)
    allocate(dst%atoms(src%natoms))
    if (src%natoms > 0) dst%atoms = src%atoms

    call cube_alloc(dst)   ! allocates field (zeroed); atoms already set
  end subroutine cube_copy_header

  ! ------------------------------------------------------------------ !

  !> Deep copy: clone all fields of `src` into a new `dst`.
  subroutine cube_clone(src, dst)
    type(CubeFile), intent(in)    :: src
    type(CubeFile), intent(inout) :: dst

    call cube_copy_header(src, dst)
    dst%source_file = src%source_file

    if (allocated(src%field)) then
      dst%field = src%field
    end if
  end subroutine cube_clone

end module cube_data