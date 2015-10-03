! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution

! Global variables
! Note: This is a program-independent module
module globvars
  implicit none

  ! Double precision type / formatting
  integer, parameter :: dp = kind(0.d0)
  real(dp), parameter :: eps_dp = epsilon(1.0_dp)
  character(*), parameter :: dp_format = "(es25.16e3)"
  character(*), parameter :: dp_format_raw = "es25.16e3"

  ! Numerical constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: sqrt_pi = sqrt(pi)
  real(dp), parameter :: tau = 2.0_dp * pi
  real(dp), parameter :: e = exp(1.0_dp)
  complex(dp), parameter :: j = (0.0_dp, 1.0_dp)

end module globvars
