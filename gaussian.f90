! Copyright (c) 2015 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file at the top-level directory of this distribution.

! 2-dimensional time-dependent Gaussian wavepacket propagation module
module gaussian
  ! We assume that psi(x,0) is a Gaussian wavepacket such that:
  ! psi(x,0) = (1 / pi sig_x^2)^(1/4) exp(-1/2 ( (x-x0) / sig_x)^2) exp(i p0_x (x - x0) / hb)

  use globvars
  use progvars, only: hbar, m, En, phi, sig_x, sig_y, x0, y0

  implicit none

  private

  public :: gaussian_init
  public :: gaussian_cleanup
  public :: gaussian_xyt
  public :: gaussian_px
  public :: gaussian_py

  ! Precalculated constants
  real(dp) :: p0_x, p0_y
  real(dp) :: p0_x_m, p0_y_m
  real(dp) :: p0_x_2m, p0_y_2m

  real(dp) :: sig_x2, sig_y2

  real(dp) :: sig_p_x, sig_p_y
  real(dp) :: sig_p_x2, sig_p_y2
  real(dp) :: norm_p_x, norm_p_y

  complex(dp) :: j_hb, jhb_m

contains

  subroutine gaussian_init()

    ! Precalculate some constants
    p0_x = sqrt(2.0_dp * m * En) * cos(phi)
    p0_y = sqrt(2.0_dp * m * En) * sin(phi)

    p0_x_m = p0_x / m
    p0_x_2m = p0_x_m / 2.0_dp
    p0_y_m = p0_y / m
    p0_y_2m = p0_y_m / 2.0_dp

    sig_x2 = sig_x**2
    sig_y2 = sig_y**2

    j_hb = j / hbar
    jhb_m = j * hbar / m

    !write(*,*) p0_x, p0_y

    sig_p_x = hbar / sig_x
    sig_p_y = hbar / sig_y

    sig_p_x2 = sig_p_x**2
    sig_p_y2 = sig_p_y**2

    norm_p_x = sqrt(1.0_dp / (sig_x * sqrt_pi))
    norm_p_y = sqrt(1.0_dp / (sig_y * sqrt_pi))

  end subroutine gaussian_init

  subroutine gaussian_cleanup()
  end subroutine gaussian_cleanup

  real(dp) function gaussian_px(p_x) result(val)
    real(dp), intent(in) :: p_x

    real(dp) :: exp_p_x
    ! Taken from Sakurai, p. 58

    exp_p_x = exp(- 0.5_dp * ( (p_x - p0_x) / sig_p_x2)**2 )

    val = norm_p_x * exp_p_x
  end function gaussian_px

  real(dp) function gaussian_py(p_y) result(val)
    real(dp), intent(in) :: p_y

    real(dp) :: exp_p_y
    ! Taken from Sakurai, p. 58

    exp_p_y = exp(- 0.5_dp * ( (p_y - p0_y) / sig_p_y2)**2 )

    val = norm_p_y * exp_p_y
  end function gaussian_py

  ! This gives the two-dimensional momentum space wfunc, but we want the 1D ones?

  complex(dp) function gaussian_xyt(x, y, t) result(val)
    real(dp), intent(in) :: x
    real(dp), intent(in) :: y
    real(dp), intent(in) :: t

    complex(dp) :: norm_x, norm_x2
    complex(dp) :: norm_y, norm_y2

    complex(dp) :: exp_x1, exp_x2
    complex(dp) :: exp_y1, exp_y2
    complex(dp) :: exp_x, exp_y

    ! Taken from Shankar, p. 154.

    norm_x2 = sqrt_pi * (sig_x + jhb_m / sig_x * t)
    norm_x = 1 / sqrt(norm_x2)

    norm_y2 = sqrt_pi * (sig_y + jhb_m / sig_y * t)
    norm_y = 1 / sqrt(norm_y2)

    exp_x1 = - ((x - x0) - p0_x_m * t)**2 / &
         (2.0_dp * sig_x2 * (1 + jhb_m / sig_x2 * t))
    exp_y1 = - ((y - y0) - p0_y_m * t)**2 / &
         (2.0_dp * sig_y2 * (1 + jhb_m / sig_y2 * t))

    exp_x2 = j_hb * p0_x * ((x - x0) - p0_x_2m * t)
    exp_y2 = j_hb * p0_y * ((y - y0) - p0_y_2m * t)

    exp_x = exp(exp_x1 + exp_x2)
    exp_y = exp(exp_y1 + exp_y2)

    val = norm_x * norm_y * exp_x * exp_y

  end function gaussian_xyt

end module gaussian
