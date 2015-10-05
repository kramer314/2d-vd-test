! Run parameters
module params
  use progvars

  implicit none

  private

  public :: params_init
  public :: params_cleanup

contains
  subroutine params_init()
    ! Units
    hbar = 1.0_dp

    ! Particle parameters
    m = 1.0_dp

    ! Gaussian parameters
    En = 10.0_dp
    phi = pi / 4.0_dp ! Equal momentum splitting
    !phi = 0.0_dp
    sig_x = 1.0_dp
    sig_y = 1.0_dp
    x0 = 0.0_dp
    y0 = 0.0_dp

    ! Numerical grid parameters
    x_min = -10.0_dp
    x_max = 50_dp
    nx = int(1e3)

    y_min = -10.0_dp
    y_max = 50.0_dp
    ny = int(1e3)

    ! External ("non-interaction") part of the grid
    nxl_external = 100
    nxr_external = 100
    nyl_external = 100
    nyr_external = 100

    ! # of grid points used as virtual detector
    ! (placed right outside the interaction region)
    vd_nxl = 1
    vd_nxr = 5
    vd_nyl = 1
    vd_nyr = 5

    ! Virtual detector binning parameters
    vd_kx_min = 0.0_dp
    vd_kx_max = 10.0_dp
    vd_nkx = int(1e2)

    vd_ky_min = 0.0_dp
    vd_ky_max = 10.0_dp
    vd_nky = int(1e2)

    ! Time parameters
    t_min = 0.0_dp
    t_max = 20.0_dp
    nt = int(5e2)

    ! Output parameters
    output_dir = "./output/"
    nkx_fname = "kxt.dat"
    nky_fname = "kyt.dat"
    output_kx = .true.
    output_ky = .true.

    ! Print filters
    print_mod_t = 10
    print_mod_kx = 10
    print_mod_ky = 10
  end subroutine params_init

  subroutine params_cleanup()
  end subroutine params_cleanup
end module params
