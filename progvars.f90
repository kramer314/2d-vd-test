! Copyright (c) 2015 Alex Kramer <alkramer@phys.ksu.edu>
! See the LICENSE.txt file in the top-level directory of this distribution

! Program-specific variables
module progvars
  use globvars

  implicit none

  ! Units
  real(dp) :: hbar

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
  integer(dp) :: nx, ny

  ! Virtual detector grid parameters
  ! Left / right number of grid points outside region of interest
  integer(dp) :: nxl_external, nxr_external, nyl_external, nyr_external
  ! Left / right number of virtual detector grid points in external grid
  integer(dp) :: vd_nxl, vd_nxr, vd_nyl, vd_nyr

  ! Useful VD grid indices (initialized by the VD module)
  integer(dp) :: vd_xl_min, vd_xl_max, vd_xr_min, vd_xr_max
  integer(dp) :: vd_yl_min, vd_yl_max, vd_yr_min, vd_yr_max
 
  ! Virtual detector momentum bin parameters
  real(dp) :: vd_px_min, vd_px_max, vd_dpx
  real(dp) :: vd_py_min, vd_py_max, vd_dpy
  integer(dp) :: vd_npx, vd_npy

  ! Time grid parameters
  real(dp) :: t_min, t_max, dt
  integer(dp) :: nt

  ! Output parameters
  character(120) :: output_dir
  character(120) :: npx_fname, npy_fname
  logical :: output_px, output_py
  integer(dp) :: print_mod_t, print_mod_px, print_mod_py

  ! Arrays
  real(dp), allocatable :: x_range(:), y_range(:), t_range(:)
  real(dp), allocatable :: vd_px_arr(:), vd_py_arr(:)
  real(dp), allocatable :: npx_arr(:), npy_arr(:)
  complex(dp), allocatable :: psi_arr(:,:)

end module progvars
