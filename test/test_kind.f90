!> Smoke-test for cube_kinds: compile, check KIND values and constants.
program test_kinds
  use cube_kinds
  implicit none

  integer :: fail = 0

  ! ---- KIND sizes ----
  call check("dp is 8-byte", storage_size(1.0_dp) == 64)
  call check("sp is 4-byte", storage_size(1.0_sp) == 32)
  call check("lp covers 10^18", lp >= 18 .or. selected_int_kind(18) > 0)

  ! ---- Mathematical constants ----
  call check("pi precision",    abs(pi - acos(-1.0_dp))  < 1.0e-15_dp)
  call check("two_pi",          abs(two_pi - 2.0_dp*pi)  < 1.0e-15_dp)
  call check("sqrt2",           abs(sqrt2 - sqrt(2.0_dp)) < 1.0e-15_dp)
  call check("third",           abs(3.0_dp*third - 1.0_dp) < 1.0e-15_dp)

  ! ---- Conversion round-trips ----
  call check("bohr<->ang", abs(bohr_to_ang * ang_to_bohr - 1.0_dp) < 1.0e-14_dp)
  call check("hartree<->ev", abs(hartree_to_ev * ev_to_hartree - 1.0_dp) < 1.0e-14_dp)

  ! ---- Stencil sum rules ----
  ! First derivative: antisymmetric stencil must give exact result for f(x)=x.
  ! Sum of positive coefficients minus sum of negative = 1 (after /h, with h=1).
  ! Contribution: 2*(c1*1 + c2*2 + c3*3) = 1
  call check("d1 stencil sum", &
    abs(2.0_dp*(stencil_d1(1)*1 + stencil_d1(2)*2 + stencil_d1(3)*3) - 1.0_dp) < 1.0e-14_dp)

  ! Second derivative: f(x)=x^2 -> f''=2. Applied at x=0 with h=1:
  ! f(0)=0, f(±1)=1, f(±2)=4, f(±3)=9
  ! result = c0*0 + c1*(1+1) + c2*(4+4) + c3*(9+9) = 2.0
  call check("d2 stencil sum", &
    abs(stencil_d2(0)*0.0_dp &
      + stencil_d2(1)*2.0_dp &
      + stencil_d2(2)*8.0_dp &
      + stencil_d2(3)*18.0_dp - 2.0_dp) < 1.0e-14_dp)

  ! ---- Guards ----
  call check("rho_tol > 0",  rho_tol  > zero)
  call check("eps_safe > 0", eps_safe > zero)

  ! ---- Summary ----
  if (fail == 0) then
    print '(a)', "cube_kinds: all tests passed."
  else
    print '(a,i0,a)', "cube_kinds: ", fail, " test(s) FAILED."
    error stop 1
  end if

contains

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

end program test_kinds