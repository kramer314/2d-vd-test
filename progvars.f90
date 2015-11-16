! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution

! Program-specific variables
module progvars
  use globvars

  implicit none

  ! Units
  real(dp) :: hbar
  real(dp) :: h

  ! Particle parameters
  real(dp) :: m

  ! Gaussian parameters
  ! Total energy (En)
  real(dp) :: En
  ! Momentum splitting - p_x = p cos(phi), p_y = p sin(phi)
  real(dp) :: phi
  ! Gaussian x, y standard deviations
  real(dp) :: sig_x, sig_y
  ! Gaussian x, y initial central points
  real(dp) :: x0, y0

  ! Grid parameters
  real(dp) :: x_min, x_max, dx
  real(dp) :: y_min, y_max, dy
  integer :: nx, ny

  ! Virtual detector grid parameters
  ! Left / right number of grid points outside region of interest
  integer :: nxl_external, nxr_external, nyl_external, nyr_external
  ! Left / right number of virtual detector grid points in external grid
  integer :: vd_nxl, vd_nxr, vd_nyl, vd_nyr

  ! Useful VD grid indices (initialized by the VD module)
  integer :: vd_xl_min, vd_xl_max, vd_xr_min, vd_xr_max
  integer :: vd_yl_min, vd_yl_max, vd_yr_min, vd_yr_max

  ! Virtual detector momentum bin parameters
  real(dp) :: vd_px_min, vd_px_max, vd_dpx
  real(dp) :: vd_py_min, vd_py_max, vd_dpy
  integer :: vd_npx, vd_npy

  ! Time grid parameters
  real(dp) :: t_min, t_max, dt
  integer :: nt

  ! Residual analysis parameters
  real(dp) :: resid_x_eps, resid_y_eps

  ! Output parameters
  character(120) :: output_dir
  character(120) :: log_fname, params_fname
  character(120) :: vd_px_fname, vd_py_fname, vd_stats_fname
  logical :: output_px, output_py
  integer :: print_mod_t, print_mod_px, print_mod_py

  ! Arrays
  real(dp), allocatable :: x_range(:), y_range(:), t_range(:)
  real(dp), allocatable :: vd_px_arr(:), vd_py_arr(:)
  real(dp), allocatable :: npx_arr(:), npy_arr(:)
  complex(dp), allocatable :: psi_arr(:,:)

  ! Theoretical result / residual analysis arrays
  real(dp), allocatable :: theor_npx_arr(:), theor_npy_arr(:)
  real(dp), allocatable :: resid_npx_arr(:), resid_npy_arr(:)
  logical, allocatable :: resid_npx_mask(:), resid_npy_mask(:)

end module progvars
