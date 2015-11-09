! Copyright (c) 2015 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file at the top-level directory of this distribution.

! Basic statistics functionality
module stats
  use globvars
  use arrays, only: arrays_quicksort

  implicit none

  private

  public :: stats_residuals
  public :: stats_mean
  public :: stats_median
  public :: stats_mean_sq_err
  public :: stats_mean_abs_err
  public :: stats_variance
  public :: stats_stdev

contains
  ! Calculate residuals
  subroutine stats_residuals(data_arr, theor_arr, resid_arr)
    real(dp), intent(in) :: data_arr(:), theor_arr(:)
    real(dp), intent(inout) :: resid_arr(:)

    resid_arr(:) = abs(data_arr(:) - theor_arr(:))

  end subroutine stats_residuals

  ! Calculate mean
  real(dp) function stats_mean(data_arr, mask) result(val)
    real(dp), intent(in) :: data_arr(:)
    logical, intent(in), optional :: mask(:)

    integer :: n

    if (present(mask)) then
       n = count(mask)
       val = sum(data_arr, mask=mask) / n
    else
       n = size(data_arr)
       val = sum(data_arr) / n
    end if

  end function stats_mean

  ! Calculate median
  real(dp) function stats_median(data_arr, mask) result(val)

    real(dp), intent(in) :: data_arr(:)
    logical, intent(in), optional :: mask(:)

    real(dp), allocatable :: data_work_arr(:)
    integer :: data_size

    if (present(mask)) then
       data_size = count(mask)
       allocate(data_work_arr(data_size))
       data_work_arr = pack(data_arr, mask)
    else
       data_size = size(data_arr)
       allocate(data_work_arr(data_size))
       data_work_arr(:) = data_arr(:)
    end if

    call arrays_quicksort(data_work_arr, 1, size(data_work_arr))

    if (mod(data_size, 2) .eq. 0) then
       val = 0.5_dp * (data_work_arr(data_size / 2) + &
            data_work_arr(data_size / 2 + 1))
    else
       val = data_work_arr( (data_size - 1) / 2 + 1 )
    end if

    deallocate(data_work_arr)

  end function stats_median

  ! Calculate mean absolute error from residuals
  real(dp) function stats_mean_abs_err(resid_arr, mask) result(val)
    real(dp), intent(in) :: resid_arr(:)
    logical, intent(in), optional :: mask(:)

    integer :: n

    if (present(mask)) then
       n = count(mask)
       val = sum(resid_arr, mask=mask) / n
    else
       n = size(resid_arr)
       val = sum(resid_arr) / n
    end if

  end function stats_mean_abs_err

  ! Calculate mean squared error from residuals
  real(dp) function stats_mean_sq_err(resid_arr, mask) result(val)
    real(dp), intent(in) :: resid_arr(:)
    logical, intent(in), optional :: mask(:)


    integer :: n

    if (present(mask)) then
       n = count(mask)
       val = sum(resid_arr**2, mask=mask) / n
    else
       n = size(resid_arr)
       val = sum(resid_arr**2) / n
    end if

  end function stats_mean_sq_err

  ! Calculate variance
  real(dp) function stats_variance(data_arr, mask) result(val)
    real(dp), intent(in) :: data_arr(:)
    logical, intent(in), optional :: mask(:)

    real(dp) :: mean
    integer :: n


    if (present(mask)) then
       n = count(mask)
       mean = stats_mean(data_arr, mask=mask)
       ! Biased sample variance
       val = sum(data_arr**2, mask=mask) / n - mean**2
       ! Convert to unbiased sample variance
       val = n / (n - 1.0_dp) * val
    else
       n = size(data_arr)
       mean = stats_mean(data_arr)
       val = sum((data_arr - mean)**2) / (n - 1)
    end if

  end function stats_variance

  ! Compute standard deviation
  real(dp) function stats_stdev(data_arr, mask) result(val)
    real(dp), intent(in) :: data_arr(:)
    logical, intent(in), optional :: mask(:)

    if (present(mask)) then
       val = sqrt(stats_variance(data_arr, mask=mask))
    else
       val = sqrt(stats_variance(data_arr))
    end if

  end function stats_stdev
end module stats
