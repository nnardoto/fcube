# fcube

A Fortran 2018 library for reading, writing, and operating on
[Gaussian `.cube`](https://paulbourke.net/dataformats/cube/) volumetric
data files.

Typical use cases: electron density analysis, ELF, NCI/RDG, electrostatic
potentials, and any scalar field on a 3-D grid produced by Gaussian, ORCA,
CP2K, VASP, or Multiwfn.

---

## Status

Work in progress — foundation complete, analysis modules upcoming.

| Module | Description | Status |
|---|---|---|
| `cube_kinds` | Precision kinds, physical constants, stencil coefficients | done |
| `cube_data` | `CubeFile` type, allocation, geometry helpers | done |
| `cube_io` | Read/write `.cube` files | done |
| `cube_arith` | Arithmetic operators and field functions | done |
| `cube_integrate` | Volumetric integration | planned |
| `cube_diff` | Gradient, Laplacian, Hessian (6th-order FD + FFT) | planned |
| `cube_elf` | Electron Localization Function | planned |
| `cube_nci` | Non-Covalent Interactions / RDG | planned |
| `cube_lol` | Localized Orbital Locator and related | planned |

---

## Requirements

- A Fortran 2018 compiler: `gfortran ≥ 9`, `ifort`, `ifx`, or `nagfor`
- GNU `make`
- No external dependencies

---

## Building

```bash
# Release build (default)
make

# Debug build — bounds checking, backtraces
make BUILD=debug

# Intel compiler
make FC=ifort

# Run all tests
make test

# Clean object files and test binaries
make clean

# Remove everything including the library archive
make distclean
```

The build system writes all generated files (`*.o`, `*.mod`,
`libfcube.a`, test binaries) to `build/`. Source directories are
never modified.

---

## Quick start

```fortran
program example
  use cube_kinds, only: dp, PERIODIC_ALL, PERIODIC_NONE
  use cube_data,  only: CubeFile, cube_free
  use cube_io,    only: read_cube, write_cube
  use cube_arith, only: operator(-)
  implicit none

  type(CubeFile) :: rho_alpha, rho_beta, spin_density

  ! Read two spin-density cubes (periodic bulk)
  call read_cube("rho_alpha.cube", rho_alpha)
  call read_cube("rho_beta.cube",  rho_beta)

  ! Spin density = α − β  (operator overloading, grid checked automatically)
  spin_density = rho_alpha - rho_beta

  call write_cube("spin_density.cube", spin_density)

  call cube_free(rho_alpha)
  call cube_free(rho_beta)
  call cube_free(spin_density)
end program
```

Compile against the library:

```bash
gfortran -std=f2018 -Ibuild example.f90 -Lbuild -lfcube -o example
```

---

## Periodicity

The `.cube` format carries no periodicity metadata. Supply it at read
time via the optional `periodicity` argument:

```fortran
use cube_data, only: PERIODIC_ALL, PERIODIC_NONE

! Fully periodic (default — most production use cases)
call read_cube("bulk.cube",  cb)

! Molecular / cluster (non-periodic)
call read_cube("mol.cube",   cb, [PERIODIC_NONE, PERIODIC_NONE, PERIODIC_NONE])

! Slab: periodic in x and y, open in z
call read_cube("slab.cube",  cb, [PERIODIC_ALL, PERIODIC_ALL, PERIODIC_NONE])

! Wire: periodic in x only
call read_cube("wire.cube",  cb, [PERIODIC_ALL, PERIODIC_NONE, PERIODIC_NONE])
```

Periodicity is stored in `cb%periodicity(3)` and propagates through
`cube_clone` and arithmetic operations. It is **not** written to disk
(the format has no slot for it), so the caller must supply it again
when re-reading a file.

---

## Module reference

### `cube_kinds`

Precision kinds and physical constants (CODATA 2018). Import from here
rather than writing literals directly.

```fortran
use cube_kinds, only: dp, sp, ip, lp          ! kind parameters
use cube_kinds, only: pi, two_pi, bohr_to_ang  ! constants
use cube_kinds, only: stencil_d1, stencil_d2   ! 6th-order FD coefficients
use cube_kinds, only: rho_tol, eps_safe        ! numerical guards
```

### `cube_data`

The `CubeFile` derived type and its lifecycle routines.

```fortran
type(CubeFile) :: cb

! Key fields
cb%nx, cb%ny, cb%nz          ! grid dimensions
cb%origin(3)                 ! grid origin in Bohr
cb%axes(3,3)                 ! step vectors in Bohr (columns)
cb%periodicity(3)            ! per-axis: PERIODIC_ALL or PERIODIC_NONE
cb%natoms                    ! number of atoms
cb%atoms(i)%z                ! atomic number
cb%atoms(i)%pos(3)           ! position in Bohr
cb%field(nx,ny,nz)           ! scalar field (iz fastest, matching .cube order)

! Routines
call cube_alloc(cb)          ! allocate atoms + field after setting dimensions
call cube_free(cb)           ! deallocate and reset
call cube_copy_header(a, b)  ! copy geometry into b, allocate b%field zeroed
call cube_clone(a, b)        ! deep copy including field data
cube_is_alloc(cb)            ! .true. if field is allocated
cube_compatible(a, b)        ! .true. if grids match within tolerance
cube_voxel_vol(cb)           ! voxel volume in Bohr³ (det of axes matrix)
cube_total_vol(cb)           ! total grid volume in Bohr³
cube_step(cb, iax)           ! step length along axis iax (1, 2, or 3)
```

### `cube_io`

```fortran
! Read — periodicity defaults to all-periodic
call read_cube(filename, cb)
call read_cube(filename, cb, periodicity=[1,1,0])

! Write — all lengths in Bohr, 6 values per line (Gaussian convention)
call write_cube(filename, cb)
```

Reads list-directed data so it tolerates the varying line-break
conventions of different generators (Gaussian, ORCA, CP2K, VASP,
Multiwfn). Lengths in Ångström (`n_axis < 0`) are converted to Bohr
on read and written back as Bohr.

Files with `natoms < 0` (multi-field / orbital cubes) are not yet
supported and produce a descriptive `error stop`.

### `cube_arith`

Operator overloading for `CubeFile` scalar fields. All binary operators
require compatible grids; an incompatible pair calls `error stop` with
a message identifying which dimension or vector differs.

```fortran
c = a + b          ! element-wise addition
c = a - b          ! element-wise subtraction
c = a * b          ! element-wise product
c = a / b          ! element-wise division (zero denominator → zero)

c = a + 2.0_dp     ! scalar offset
c = 3.0_dp * a     ! scalar scale
c = -a             ! negation

c = cube_abs(a)    ! |field|
c = cube_sqrt(a)   ! sqrt(field), negative voxels clamped to zero

call cube_apply(a, f)        ! apply pure function f(x) in-place
x = cube_max_val(a)          ! maximum voxel value
x = cube_min_val(a)          ! minimum voxel value
```

---

## Design notes

**Unit convention.** All internal values are in atomic units (Bohr,
Hartree, electron charge). Conversion constants live in `cube_kinds`;
no silent conversions happen inside the library.

**Memory layout.** `field(ix, iy, iz)` — `iz` is the fastest index,
matching the `.cube` disk order (z varies fastest). Read and write are
sequential with no transposition.

**Numerical differentiation** (upcoming `cube_diff`). Non-periodic axes
use 6th-order centred finite differences (7-point stencil). Periodic
axes use FFT-based derivatives, exact within the sampling. The
per-axis `periodicity` vector drives the choice automatically.

**Error handling.** All errors call `error stop` with a descriptive
message. There are no silent fallbacks or iostat outputs.

---

## Repository layout

```
fcube/
├── Makefile
├── README.md
├── src/
│   ├── cube_kinds.f90
│   ├── cube_data.f90
│   ├── cube_io.f90
│   ├── cube_arith.f90
│   └── analysis/          ← ELF, NCI, LOL (upcoming)
├── test/
│   ├── test_kinds.f90
│   ├── test_data.f90
│   ├── test_io.f90
│   ├── test_arith.f90
│   └── water_test.cube
└── build/                 ← generated, exclude from version control
```

---

## License

To be decided.
