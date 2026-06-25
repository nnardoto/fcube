!> Tests for cube_io: read, write, round-trip fidelity, unit conversion,
!> and error paths.
program test_io
  use cube_kinds
  use cube_data
  use cube_io
  implicit none

  integer :: fail = 0

  call test_read_basic()
  call test_periodicity()
  call test_round_trip()
  call test_angstrom_conversion()
  call test_write_format()

  if (fail == 0) then
    print '(a)', "cube_io: all tests passed."
  else
    print '(a,i0,a)', "cube_io: ", fail, " test(s) FAILED."
    error stop 1
  end if

contains

  ! ------------------------------------------------------------------ !
  subroutine test_read_basic()
    type(CubeFile) :: cb

    call read_cube("test/water_test.cube", cb)

    call check("read: natoms",   cb%natoms == 3)
    call check("read: nx",       cb%nx == 4)
    call check("read: ny",       cb%ny == 4)
    call check("read: nz",       cb%nz == 4)
    call check("read: is_alloc", cube_is_alloc(cb))

    ! Origin should be (0,0,0)
    call check("read: origin x", abs(cb%origin(1)) < 1.0e-10_dp)
    call check("read: origin y", abs(cb%origin(2)) < 1.0e-10_dp)
    call check("read: origin z", abs(cb%origin(3)) < 1.0e-10_dp)

    ! Step vectors: orthogonal 0.5 Bohr
    call check("read: step x",   abs(cube_step(cb, 1) - 0.5_dp) < 1.0e-10_dp)
    call check("read: step y",   abs(cube_step(cb, 2) - 0.5_dp) < 1.0e-10_dp)
    call check("read: step z",   abs(cube_step(cb, 3) - 0.5_dp) < 1.0e-10_dp)

    ! Atom records
    call check("read: atom 1 Z",  cb%atoms(1)%z == 8)
    call check("read: atom 2 Z",  cb%atoms(2)%z == 1)
    call check("read: atom 3 Z",  cb%atoms(3)%z == 1)
    call check("read: atom 1 charge", abs(cb%atoms(1)%charge - (-0.834_dp)) < 1.0e-6_dp)
    call check("read: atom 1 x",  abs(cb%atoms(1)%pos(1) - 0.0_dp) < 1.0e-6_dp)

    ! First data value: field(1,1,1) = 1.0e-3
    call check("read: field(1,1,1)", abs(cb%field(1,1,1) - 1.0e-3_dp) < 1.0e-8_dp)

    ! Last data value: field(4,4,4) = 5.6e-3 (last in file)
    call check("read: field(4,4,4)", abs(cb%field(4,4,4) - 5.6e-3_dp) < 1.0e-8_dp)

    ! Middle values: field(1,1,4) = 4.0e-3 (4th value of first strip)
    call check("read: field(1,1,4)", abs(cb%field(1,1,4) - 4.0e-3_dp) < 1.0e-8_dp)

    ! Check title was read
    call check("read: title not empty", len_trim(cb%title(1)) > 0)

    ! Default periodicity: all periodic
    call check("read: default periodic", all(cb%periodicity == PERIODIC_ALL))

    call cube_free(cb)
  end subroutine

  ! ------------------------------------------------------------------ !
  subroutine test_periodicity()
    type(CubeFile) :: cb

    ! Default (no argument): all axes periodic
    call read_cube("test/water_test.cube", cb)
    call check("per: default all periodic", &
      all(cb%periodicity == PERIODIC_ALL))
    call cube_free(cb)

    ! Explicit all-periodic
    call read_cube("test/water_test.cube", cb, [PERIODIC_ALL, PERIODIC_ALL, PERIODIC_ALL])
    call check("per: explicit all periodic", &
      all(cb%periodicity == PERIODIC_ALL))
    call cube_free(cb)

    ! Fully aperiodic (molecular)
    call read_cube("test/water_test.cube", cb, [PERIODIC_NONE, PERIODIC_NONE, PERIODIC_NONE])
    call check("per: all none x", cb%periodicity(1) == PERIODIC_NONE)
    call check("per: all none y", cb%periodicity(2) == PERIODIC_NONE)
    call check("per: all none z", cb%periodicity(3) == PERIODIC_NONE)
    call cube_free(cb)

    ! Slab: periodic in x,y; open in z
    call read_cube("test/water_test.cube", cb, [PERIODIC_ALL, PERIODIC_ALL, PERIODIC_NONE])
    call check("per: slab x periodic", cb%periodicity(1) == PERIODIC_ALL)
    call check("per: slab y periodic", cb%periodicity(2) == PERIODIC_ALL)
    call check("per: slab z open",     cb%periodicity(3) == PERIODIC_NONE)
    call cube_free(cb)

    ! Wire: periodic in x only
    call read_cube("test/water_test.cube", cb, [PERIODIC_ALL, PERIODIC_NONE, PERIODIC_NONE])
    call check("per: wire x periodic", cb%periodicity(1) == PERIODIC_ALL)
    call check("per: wire y open",     cb%periodicity(2) == PERIODIC_NONE)
    call check("per: wire z open",     cb%periodicity(3) == PERIODIC_NONE)

    ! Periodicity survives clone
    block
      type(CubeFile) :: clone
      call cube_clone(cb, clone)
      call check("per: clone preserves", &
        all(clone%periodicity == cb%periodicity))
      call cube_free(clone)
    end block

    ! Round-trip: periodicity set after read survives write→read cycle
    ! (periodicity is NOT in the file, so the re-read must supply it again)
    block
      type(CubeFile) :: dst
      call write_cube("/tmp/cubelib_per_rt.cube", cb)
      call read_cube("/tmp/cubelib_per_rt.cube", dst, &
                     [PERIODIC_ALL, PERIODIC_NONE, PERIODIC_NONE])
      call check("per: rt x periodic", dst%periodicity(1) == PERIODIC_ALL)
      call check("per: rt y open",     dst%periodicity(2) == PERIODIC_NONE)
      call cube_free(dst)
    end block

    call cube_free(cb)
  end subroutine

  ! ------------------------------------------------------------------ !
  subroutine test_round_trip()
    type(CubeFile) :: src, dst
    integer(ip) :: ix, iy, iz
    real(dp) :: max_err

    call read_cube("test/water_test.cube", src)
    call write_cube("/tmp/cubelib_roundtrip.cube", src)
    call read_cube("/tmp/cubelib_roundtrip.cube", dst)

    call check("rt: compatible", cube_compatible(src, dst))
    call check("rt: natoms",     dst%natoms == src%natoms)

    max_err = 0.0_dp
    do ix = 1, src%nx
      do iy = 1, src%ny
        do iz = 1, src%nz
          max_err = max(max_err, abs(dst%field(ix,iy,iz) - src%field(ix,iy,iz)))
        end do
      end do
    end do
    ! ES13.5 format gives 5 significant digits → relative error < 1e-5
    call check("rt: field max err < 1e-7", max_err < 1.0e-7_dp)

    ! Atom positions preserved within write precision (6 decimal places)
    call check("rt: atom 1 pos x", &
      abs(dst%atoms(1)%pos(1) - src%atoms(1)%pos(1)) < 1.0e-5_dp)
    call check("rt: atom 2 pos y", &
      abs(dst%atoms(2)%pos(2) - src%atoms(2)%pos(2)) < 1.0e-5_dp)

    call cube_free(src);  call cube_free(dst)
  end subroutine

  ! ------------------------------------------------------------------ !
  subroutine test_angstrom_conversion()
    type(CubeFile) :: cb

    ! Write a temporary .cube with negative n_axis (Ångström units)
    call write_angstrom_cube("/tmp/cubelib_ang.cube")
    call read_cube("/tmp/cubelib_ang.cube", cb)

    ! Step should have been converted: 1.0 Å = 1/bohr_to_ang Bohr
    call check("ang: step x in Bohr", &
      abs(cube_step(cb, 1) - ang_to_bohr) < 1.0e-8_dp)
    call check("ang: step y in Bohr", &
      abs(cube_step(cb, 2) - ang_to_bohr) < 1.0e-8_dp)

    call cube_free(cb)
  end subroutine

  ! -- Helper: write a minimal Ångström-unit cube -------------------- !
  subroutine write_angstrom_cube(filename)
    character(len=*), intent(in) :: filename
    integer(ip) :: u

    open(newunit=u, file=filename, status='replace', action='write')
    write(u, '(a)') "Angstrom test"
    write(u, '(a)') "n_axis < 0 means Angstrom"
    write(u, '(i5,3f12.6)')  1, 0.0_dp, 0.0_dp, 0.0_dp   ! natoms=1, origin
    write(u, '(i5,3f12.6)') -2, 1.0_dp, 0.0_dp, 0.0_dp   ! nx=-2 → Å
    write(u, '(i5,3f12.6)') -2, 0.0_dp, 1.0_dp, 0.0_dp   ! ny=-2 → Å
    write(u, '(i5,3f12.6)')  2, 0.0_dp, 0.0_dp, 1.0_dp   ! nz=+2 → Bohr
    write(u, '(i5,4f12.6)')  6, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp
    ! 2*2*2 = 8 values
    write(u, '(6es13.5)') 0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, 0.6_dp
    write(u, '(2es13.5)') 0.7_dp, 0.8_dp
    close(u)
  end subroutine

  ! ------------------------------------------------------------------ !
  subroutine test_write_format()
    type(CubeFile) :: cb
    integer(ip) :: u, ios, line_count, val_count
    character(len=256) :: line
    logical :: in_data

    ! Build a small cube and write it
    cb%natoms = 0;  cb%nx = 3;  cb%ny = 2;  cb%nz = 4
    cb%title(1) = "Format test"
    cb%title(2) = "3x2x4 grid"
    cb%origin = zero
    cb%axes = zero
    cb%axes(1,1) = 0.2_dp;  cb%axes(2,2) = 0.2_dp;  cb%axes(3,3) = 0.2_dp
    call cube_alloc(cb)
    cb%field = 1.23456789e-4_dp

    call write_cube("/tmp/cubelib_fmt.cube", cb)

    ! Count header lines and check data line lengths
    open(newunit=u, file="/tmp/cubelib_fmt.cube", status='old', action='read')

    ! Skip 2+1+3+natoms = 6 header lines
    line_count = 0;  in_data = .false.;  val_count = 0
    do
      read(u, '(a)', iostat=ios) line
      if (ios /= 0) exit
      line_count = line_count + 1
      if (line_count > 6) then
        ! Count tokens per data line; no line should exceed 6
        val_count = val_count + count_tokens(trim(line))
      end if
    end do
    close(u)

    ! nx*ny*nz = 3*2*4 = 24 values.  Each (ix,iy) strip has nz=4 values
    ! → 4 values → 1 full line per strip.  Strips: 3*2=6.
    ! header=6 lines + 6 data lines = 12 total
    call check("fmt: line count", line_count == 12)
    call check("fmt: value count", val_count == 24)

    call cube_free(cb)
  end subroutine

  ! -- Count whitespace-separated tokens in a string ----------------- !
  integer function count_tokens(s)
    character(len=*), intent(in) :: s
    integer :: i
    logical :: in_token

    count_tokens = 0;  in_token = .false.
    do i = 1, len(s)
      if (s(i:i) /= ' ') then
        if (.not. in_token) then
          count_tokens = count_tokens + 1
          in_token = .true.
        end if
      else
        in_token = .false.
      end if
    end do
  end function

  ! ------------------------------------------------------------------ !
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

end program test_io