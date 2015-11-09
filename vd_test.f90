! Copyright (c) 2015 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file at the top-level directory of this distribution.

program vd_test
  use progvars
  use propagate, only: propagate_psi
  use setup, only: setup_init, setup_cleanup
  use vd, only: vd_update, vd_normalize
  use gaussian, only: gaussian_px, gaussian_py
  use stats, only: stats_residuals, stats_mean, stats_variance, stats_median, &
       stats_stdev, stats_mean_sq_err, stats_mean_abs_err

  implicit none

  integer :: i_t
  integer :: i_px, i_py
  real(dp) :: px, py

  call setup_init()

  ! Print out grid setup
  write(*,*) "Full Grid: x min; x max:", x_min, x_max
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
  write(*,*)
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

  ! Calculate theoretical results and calculate errors
  do i_px = 1, size(theor_npx_arr)
     px = vd_px_arr(i_px)
     theor_npx_arr(i_px) = abs(gaussian_px(px))**2
  end do
  do i_py = 1, size(theor_npy_arr)
     py = vd_py_arr(i_py)
     theor_npy_arr(i_py) = abs(gaussian_py(py))**2
  end do

  call stats_residuals(npx_arr, theor_npx_arr, resid_npx_arr)
  call stats_residuals(npy_arr, theor_npy_arr, resid_npy_arr)

  ! Output results
  open(98, file="vd_x.dat")
  open(99, file="vd_y.dat")

  do i_px = 1, size(npx_arr)
     px = vd_px_arr(i_px)
     write(98,*) px, npx_arr(i_px), theor_npx_arr(i_px), resid_npx_arr(i_px)
  end do
  do i_py = 1, size(npy_arr)
     py = vd_py_arr(i_py)
     write(99, *) py, npy_arr(i_py), theor_npy_arr(i_py), resid_npy_arr(i_py)
  end do

  ! Residual analysis
  resid_npx_mask(:) = resid_npx_arr(:) .ge. resid_x_eps
  resid_npy_mask(:) = resid_npy_arr(:) .ge. resid_y_eps

  write(*,*)
  write(*,*) "Residual analysis for detected/theoretical momentum distributions"
  write(*,*) "Format: [Label] [x value] [y value]"
  write(*,*)
  write(*,*) "Sum of residuals:", sum(resid_npx_arr), sum(resid_npy_arr)
  write(*,*) "Sum of squared residauls:", sum(resid_npx_arr(:)**2), &
       sum(resid_npy_arr(:)**2)

  write(*,*)
  write(*,*) "Residual thresholds:", resid_x_eps, resid_y_eps
  write(*,*) "Number of above-threshold residuals", &
       count(resid_npx_mask), count(resid_npy_mask)

  write(*,*)
  write(*,*) "Mean of above-threshold residuals:", &
       stats_mean(resid_npx_arr, mask=resid_npx_mask), &
       stats_mean(resid_npy_arr, mask=resid_npy_mask)
  write(*,*) "Median of above-threshold residuals:", &
       stats_median(resid_npx_arr, mask=resid_npx_mask), &
       stats_median(resid_npy_arr, mask=resid_npy_mask)
  write(*,*) "Variance of above-threshold residuals:", &
       stats_variance(resid_npx_arr, mask=resid_npx_mask), &
       stats_variance(resid_npy_arr, mask=resid_npy_mask)
  write(*,*) "Mean squared error of above-threshold residuals:", &
       stats_mean_sq_err(resid_npx_arr, mask=resid_npx_mask), &
       stats_mean_sq_err(resid_npy_arr, mask=resid_npy_mask)
  write(*,*) "Mean absolute error of above-threshold residuals:", &
       stats_mean_abs_err(resid_npx_arr, mask=resid_npx_mask), &
       stats_mean_abs_err(resid_npy_arr, mask=resid_npy_mask)
  write(*,*) "Maximum residuals:", maxval(resid_npx_arr), maxval(resid_npy_arr)
  write(*,*) "Minimum residuals:", minval(resid_npx_arr), minval(resid_npy_arr)

  call setup_cleanup()

end program vd_test
