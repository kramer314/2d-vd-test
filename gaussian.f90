! 2-Dimensional Gaussian
module gaussian

  use globvars
  use progvars, only: hbar, m, En, phi, sig_x, sig_y, x0, y0

  implicit none
  
  private

  public :: gaussian_init
  public :: gaussian_cleanup
  public :: gaussian_xyt

  ! Precalculated constants
  real(dp) :: p0_x, p0_y
  real(dp) :: p0_x_m, p0_y_m
  real(dp) :: p0_x_2m, p0_y_2m

  real(dp) :: sig_x2, sig_y2

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
    
  end subroutine gaussian_init

  subroutine gaussian_cleanup()
  end subroutine gaussian_cleanup

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

    exp_x1 = - (x - p0_x_m * t)**2 / &
         (2.0_dp * sig_x2 * (1 + jhb_m / sig_x2 * t))
    exp_y1 = - (y - p0_y_m * t)**2 / &
         (2.0_dp * sig_y2 * (1 + jhb_m / sig_y2 * t))
    
    exp_x2 = j_hb * p0_x * (x - p0_x_2m * t)
    exp_y2 = j_hb * p0_y * (y - p0_y_2m * t)

    exp_x = exp(exp_x1 + exp_x2)
    exp_y = exp(exp_y1 + exp_y2)

    val = norm_x * norm_y * exp_x * exp_y
    
  end function gaussian_xyt
  
end module gaussian
