! Copyright (c) 2015 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file at the top-level directory of this distribution.

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
    h = 2.0_dp * pi * hbar

    ! Particle parameters
    m = 1.0_dp

    ! Gaussian parameters
    En = 12.5_dp ! p = sqrt(2mE) = 5
    phi = pi / 4.0_dp ! Equal momentum splitting
    phi = 0.0_dp
    sig_x = 1.0_dp
    sig_y = 1.0_dp
    x0 = 0.0_dp
    y0 = 0.0_dp

    ! Numerical grid parameters
    x_min = -10.0_dp
    x_max = 25_dp
    nx = int(3e3)

    y_min = -10.0_dp
    y_max = 25.0_dp
    ny = int(1e3)

    ! External ("non-interaction") part of the grid
    nxl_external = 100
    nxr_external = 100
    nyl_external = 100
    nyr_external = 100

    ! # of grid points used as virtual detector
    ! (placed right outside the interaction region)
    vd_nxl = 5
    vd_nxr = 20
    vd_nyl = 5
    vd_nyr = 5

    ! Virtual detector binning parameters
    vd_px_min = 0.0_dp
    vd_px_max = 10.0_dp
    vd_npx = int(1e2)

    vd_py_min = 0.0_dp
    vd_py_max = 10.0_dp
    vd_npy = int(1e2)

    ! Time parameters
    t_min = 0.0_dp
    t_max = 15_dp
    nt = int(1e2)

    ! Residual analysis parameters
    resid_x_eps = 1.0e-6_dp
    resid_y_eps = 1.0e-6_dp

    ! Output parameters
    output_dir = "./output/"
    params_fname = "vd_params.dat"
    vd_px_fname = "vd_px.dat"
    vd_py_fname = "vd_py.dat"
    log_fname = "log.dat"
    vd_stats_fname = "vd_stats.dat"
    output_px = .true.
    output_py = .true.

    ! Print filters
    print_mod_t = 10
    print_mod_px = 10
    print_mod_py = 10
  end subroutine params_init

  subroutine params_cleanup()
  end subroutine params_cleanup
end module params
