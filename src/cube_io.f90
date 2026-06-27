!> cube_io — read and write Gaussian .cube files
!>
!> Format reference:
!>   https://paulbourke.net/dataformats/cube/
!>   Gaussian 16 User's Reference, §C.1
!>
!> Reading strategy:
!>   - Header fields use explicit READ with format to preserve title strings.
!>   - Volumetric data uses list-directed READ (*) to tolerate the varying
!>     line breaks produced by different generators (Gaussian, ORCA, CP2K,
!>     VASP, Multiwfn).  List-directed I/O ignores line boundaries, so the
!>     6-values-per-line convention (and the extra newline after each iz-row)
!>     are all handled transparently.
!>
!> Writing strategy:
!>   - Header follows the Gaussian convention exactly (field widths, signs).
!>   - Volumetric data: 6 values per line, ES13.5 format, extra newline after
!>     each (ix,iy) block — maximising compatibility with downstream tools.
!>
!> Unit convention:
!>   n_axis > 0 → lengths already in Bohr  (standard; no conversion needed)
!>   n_axis < 0 → lengths in Ångström      (converted to Bohr on read;
!>                                           written back as Bohr with |n|)
!>
!> Periodicity:
!>   The .cube format carries no periodicity information.  The caller
!>   supplies it via the optional `periodicity` argument to read_cube:
!>
!>     call read_cube("bulk.cube",  cb)                       ! all periodic (default)
!>     call read_cube("mol.cube",   cb, [0,0,0])              ! fully aperiodic
!>     call read_cube("slab.cube",  cb, [1,1,0])              ! periodic in x,y only
!>
!>   PERIODIC_ALL (1) and PERIODIC_NONE (0) are exported from cube_data.
!>   When writing, the periodicity field of cb is preserved as metadata but
!>   not encoded in the file (the format has no slot for it).
!>
!> Not yet implemented:
!>   natoms < 0 (multi-field / orbital cube) — raises error stop with message.
module cube_io
  use cube_kinds, only: dp, ip, ang_to_bohr, zero
  use cube_data,  only: CubeFile, Atom, cube_alloc, cube_free, &
                        PERIODIC_ALL, PERIODIC_NONE
  implicit none
  private

  public :: read_cube
  public :: write_cube

  !> Values per line in the volumetric section when writing.
  integer(ip), parameter :: VALS_PER_LINE = 6

  !> Default periodicity applied when the argument is absent: all axes periodic.
  integer(ip), parameter :: DEFAULT_PERIODICITY(3) = PERIODIC_ALL

contains

  ! ================================================================== !
  !  READ
  ! ================================================================== !

  !> Read a .cube file from `filename` into `cb`.
  !>
  !> `periodicity(3)` — optional, default [1,1,1] (all periodic).
  !>   Each element is PERIODIC_ALL (1) or PERIODIC_NONE (0) for axes x,y,z.
  !>   Examples:
  !>     call read_cube("bulk.cube", cb)            ! [1,1,1]
  !>     call read_cube("mol.cube",  cb, [0,0,0])   ! fully aperiodic
  !>     call read_cube("slab.cube", cb, [1,1,0])   ! periodic xy, open z
  !>
  !> On any error, error stop is called with a descriptive message.
  subroutine read_cube(filename, cb, periodicity)
    character(len=*),  intent(in)           :: filename
    type(CubeFile),    intent(inout)        :: cb
    integer(ip),       intent(in), optional :: periodicity(3)

    integer(ip) :: u, ios
    integer(ip) :: i, ix, iy, iz
    integer(ip) :: raw_natoms
    integer(ip) :: nx_raw, ny_raw, nz_raw
    logical     :: ang_x, ang_y, ang_z
    character(len=512) :: errmsg

    ! -- Validate periodicity argument before touching the file -------- !
    if (present(periodicity)) then
      if (any(periodicity /= PERIODIC_NONE .and. periodicity /= PERIODIC_ALL)) &
        error stop "read_cube: periodicity elements must be PERIODIC_NONE (0) " // &
                   "or PERIODIC_ALL (1)"
    end if

    ! -- Open ---------------------------------------------------------- !
    open(newunit=u, file=filename, status='old', action='read', &
         iostat=ios, iomsg=errmsg)
    if (ios /= 0) &
      error stop "read_cube: cannot open file: " // trim(errmsg)

    call cube_free(cb)
    cb%source_file = filename

    ! Store periodicity now (before cube_alloc which does not touch it).
    if (present(periodicity)) then
      cb%periodicity = periodicity
    else
      cb%periodicity = DEFAULT_PERIODICITY
    end if

    ! -- Lines 1-2: title ---------------------------------------------- !
    read(u, '(a)', iostat=ios, iomsg=errmsg) cb%title(1)
    if (ios /= 0) call io_error("reading title line 1", errmsg)

    read(u, '(a)', iostat=ios, iomsg=errmsg) cb%title(2)
    if (ios /= 0) call io_error("reading title line 2", errmsg)

    ! -- Line 3: natoms + origin --------------------------------------- !
    read(u, *, iostat=ios, iomsg=errmsg) raw_natoms, cb%origin
    if (ios /= 0) call io_error("reading natoms and origin", errmsg)

    if (raw_natoms < 0) then
      close(u)
      error stop "read_cube: multi-field cube files (natoms < 0) are not yet " // &
                 "implemented. The file contains multiple orbital/spin " // &
                 "components indicated by a negative atom count. " // &
                 "Support is planned for a future version."
    end if
    cb%natoms = raw_natoms

    ! -- Lines 4-6: grid axes ----------------------------------------- !
    ! n < 0 means the step length is in Ångström → convert to Bohr.
    read(u, *, iostat=ios, iomsg=errmsg) nx_raw, cb%axes(:, 1)
    if (ios /= 0) call io_error("reading x-axis", errmsg)

    read(u, *, iostat=ios, iomsg=errmsg) ny_raw, cb%axes(:, 2)
    if (ios /= 0) call io_error("reading y-axis", errmsg)

    read(u, *, iostat=ios, iomsg=errmsg) nz_raw, cb%axes(:, 3)
    if (ios /= 0) call io_error("reading z-axis", errmsg)

    ang_x = (nx_raw < 0);  cb%nx = abs(nx_raw)
    ang_y = (ny_raw < 0);  cb%ny = abs(ny_raw)
    ang_z = (nz_raw < 0);  cb%nz = abs(nz_raw)

    if (cb%nx < 1 .or. cb%ny < 1 .or. cb%nz < 1) &
      error stop "read_cube: grid dimensions must be >= 1"

    if (ang_x) cb%axes(:, 1) = cb%axes(:, 1) * ang_to_bohr
    if (ang_y) cb%axes(:, 2) = cb%axes(:, 2) * ang_to_bohr
    if (ang_z) cb%axes(:, 3) = cb%axes(:, 3) * ang_to_bohr

    ! -- Allocate ------------------------------------------------------ !
    call cube_alloc(cb)

    ! -- Lines 7..7+natoms-1: atom records ----------------------------- !
    do i = 1, cb%natoms
      read(u, *, iostat=ios, iomsg=errmsg) &
        cb%atoms(i)%z, cb%atoms(i)%charge, cb%atoms(i)%pos
      if (ios /= 0) call io_error("reading atom record", errmsg)
    end do

    ! -- Volumetric data ----------------------------------------------- !
    ! List-directed: ignores line boundaries, tolerates any generator.
    read(u, *, iostat=ios, iomsg=errmsg) &
      (((cb%field(ix, iy, iz), iz=1,cb%nz), iy=1,cb%ny), ix=1,cb%nx)
    if (ios /= 0) call io_error("reading volumetric data", errmsg)

    close(u)
  end subroutine read_cube

  ! ================================================================== !
  !  WRITE
  ! ================================================================== !

  !> Write `cb` to `filename` in standard Gaussian .cube format.
  !> All lengths are written in Bohr (n_axis > 0).
  !> Note: periodicity metadata is NOT encoded in the file (the format
  !> has no slot for it); it lives only in the CubeFile struct.
  subroutine write_cube(filename, cb)
    character(len=*), intent(in) :: filename
    type(CubeFile),   intent(in) :: cb

    integer(ip) :: u, ios, i, ix, iy, iz, col
    character(len=512) :: errmsg

    if (.not. allocated(cb%field)) &
      error stop "write_cube: field not allocated"

    open(newunit=u, file=filename, status='replace', action='write', &
         iostat=ios, iomsg=errmsg)
    if (ios /= 0) &
      error stop "write_cube: cannot open file for writing: " // trim(errmsg)

    ! -- Lines 1-2: title ---------------------------------------------- !
    write(u, '(a)') trim(cb%title(1))
    write(u, '(a)') trim(cb%title(2))

    ! -- Line 3: natoms + origin --------------------------------------- !
    write(u, '(i5,3f12.6)') cb%natoms, cb%origin

    ! -- Lines 4-6: axes (Bohr, so n_axis > 0) ------------------------ !
    write(u, '(i5,3f12.6)') cb%nx, cb%axes(:, 1)
    write(u, '(i5,3f12.6)') cb%ny, cb%axes(:, 2)
    write(u, '(i5,3f12.6)') cb%nz, cb%axes(:, 3)

    ! -- Atom records -------------------------------------------------- !
    do i = 1, cb%natoms
      write(u, '(i5,4f12.6)') &
        cb%atoms(i)%z, cb%atoms(i)%charge, cb%atoms(i)%pos
    end do

    ! -- Volumetric data ----------------------------------------------- !
    ! 6 values per line; mandatory extra newline after each (ix,iy) strip.
    do ix = 1, cb%nx
      do iy = 1, cb%ny
        col = 0
        do iz = 1, cb%nz
          col = col + 1
          write(u, '(es13.5)', advance='no') cb%field(ix, iy, iz)
          if (col == VALS_PER_LINE) then
            write(u, '()')
            col = 0
          end if
        end do
        if (col > 0) write(u, '()')   ! flush partial line
      end do
    end do

    close(u)
  end subroutine write_cube

  ! ================================================================== !
  !  Internal helpers
  ! ================================================================== !

  subroutine io_error(context, msg)
    character(len=*), intent(in) :: context, msg
    error stop "read_cube: error while " // context // ": " // trim(msg)
  end subroutine io_error

end module cube_io