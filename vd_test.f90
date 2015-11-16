! Copyright (c) 2015 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file at the top-level directory of this distribution.

program vd_test
  use progvars
  use propagate, only: propagate_psi
  use setup, only: setup_init, setup_cleanup
  use vd, only: vd_update, vd_normalize
  use gaussian, only: gaussian_px, gaussian_py, gaussian_calc_theor_npxy
  use stats, only: stats_residuals, stats_mean, stats_variance, stats_median, &
       stats_stdev, stats_mean_sq_err, stats_mean_abs_err
  use output, only: output_vd_results, output_vd_resid_stats, output_log

  implicit none

  integer :: i_t

  write(*,*) "Initializing program"
  call setup_init()

  ! Print out grid setup
  write(*,*) "Full Grid: x_min; x_max:", x_min, x_max
  write(*,*) "Full Grid: x min; x max:",  x_min, x_max
  write(*,*) "Full Grid: y_min; y_max:", y_min, y_max
  write(*,*)

  write(*,*) "Internal Grid: x_min, x_max", x_range(nxl_external), x_range(nx - nxr_external)
  write(*,*) "Internal Grid: y_min; y_max", y_range(nyl_external), y_range(ny - nyr_external)
  write(*,*)

  write(*,*) "Left Virtual X Grid: xl_min, xl_max", x_range(vd_xl_min), x_range(vd_xl_max)
  write(*,*) "Right Virtual X Grid: xr_min, xr_max", x_range(vd_xr_min), x_range(vd_xr_max)
  write(*,*) "Left Virtual Y Grid: yr_min, yr_max", y_range(vd_yl_min), x_range(vd_yl_max)
  write(*,*) "Right Virtual Y Grid: yr_min, yr_max", y_range(vd_yr_min), x_range(vd_yr_max)
  write(*,*)

  ! Propagate in time and run virtual detector
  write(*,*) "Beginning time propagation"
  write(*,*) "Format: [current step] [total steps]"

  do i_t = 1, nt
     call propagate_psi(psi_arr, i_t)

     call vd_update(psi_arr, npx_arr, npy_arr)

     if (mod(i_t, print_mod_t) .eq. 0) then
        write(*,*) i_t, nt
     end if
  end do

  call vd_normalize(npx_arr, npy_arr)

  call gaussian_calc_theor_npxy(vd_px_arr, vd_py_arr, &
       theor_npx_arr, theor_npy_arr)

  call output_vd_results()
  call output_vd_resid_stats()

  call setup_cleanup()

end program vd_test
