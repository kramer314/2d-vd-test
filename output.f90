! Copyright (c) 2015 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file in the top-level directory of this distribution.

module output
  use progvars
  use files, only: files_ensure_dir
  use stats, only: stats_residuals, stats_mean, stats_variance, stats_median, &
       stats_stdev, stats_mean_sq_err, stats_mean_abs_err

  implicit none

  private

  ! Generic module functions
  public :: output_init
  public :: output_cleanup

  public :: output_log
  public :: output_vd_results
  public :: output_vd_resid_stats
  public :: output_params

  ! Output file unit numbers
  integer, parameter :: log_unit = 99
  integer, parameter :: params_unit = 98
  integer, parameter :: vd_px_unit = 89
  integer, parameter :: vd_py_unit = 87
  integer, parameter :: vd_stats_unit = 86

  character(10) :: date, time

contains
  ! Module initializaiton
  subroutine output_init
    call files_ensure_dir(output_dir)

    open(unit=log_unit, file=trim(output_dir)//trim(log_fname))
    open(unit=params_unit, file=trim(output_dir)//trim(params_fname))
    open(unit=vd_px_unit, file=trim(output_dir)//trim(vd_px_fname))
    open(unit=vd_py_unit, file=trim(output_dir)//trim(vd_py_fname))
    open(unit=vd_stats_unit, file=trim(output_dir)//trim(vd_stats_fname))
  end subroutine output_init

  subroutine output_log(log_text, fmt)
    character(*), intent(in) :: log_text
    character(*), intent(in), optional :: fmt

    call date_and_time(date=date, time=time)
    write(log_unit, *) date, time
    
    if (present(fmt)) then
       write(log_unit, fmt) log_text
    else
       write(log_unit, *) log_text
    end if
    
  end subroutine output_log

  subroutine output_params()
  end subroutine output_params
  
  subroutine output_vd_results()

    integer :: i_px, i_py
    real(dp) :: px, py

    ! Residual analysis
    call stats_residuals(npx_arr, theor_npx_arr, resid_npx_arr)
    call stats_residuals(npy_arr, theor_npy_arr, resid_npy_arr)

    ! Output vd results
    do i_px = 1, size(npx_arr)
       px = vd_px_arr(i_px)
       write(vd_px_unit, *) px, npx_arr(i_px), theor_npx_arr(i_px), resid_npx_arr(i_px)
    end do
    do i_py = 1, size(npy_arr)
       py = vd_py_arr(i_py)
       write(vd_py_unit, *) py, npy_arr(i_py), theor_npy_arr(i_py), resid_npy_arr(i_py)
    end do

  end subroutine output_vd_results

  subroutine output_vd_resid_stats()

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

    write(*,*)
    write(*,*) "Variance of above-threshold residuals:", &
         stats_variance(resid_npx_arr, mask=resid_npx_mask), &
         stats_variance(resid_npy_arr, mask=resid_npy_mask)
    write(*,*) "Mean squared error of above-threshold residuals:", &
         stats_mean_sq_err(resid_npx_arr, mask=resid_npx_mask), &
         stats_mean_sq_err(resid_npy_arr, mask=resid_npy_mask)
    write(*,*) "Mean absolute error of above-threshold residuals:", &
         stats_mean_abs_err(resid_npx_arr, mask=resid_npx_mask), &
         stats_mean_abs_err(resid_npy_arr, mask=resid_npy_mask)

    write(*,*)
    write(*,*) "Maximum residuals:", maxval(resid_npx_arr), maxval(resid_npy_arr)
    write(*,*) "Minimum residuals:", minval(resid_npx_arr), minval(resid_npy_arr)

  end subroutine output_vd_resid_stats

  subroutine output_cleanup
    close(unit=log_unit)
    close(unit=params_unit)
    close(unit=vd_px_unit)
    close(unit=vd_py_unit)
    close(unit=vd_stats_unit)
  end subroutine output_cleanup

end module output
