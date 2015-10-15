! Copyright (c) 2015 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file at the top-level directory of this distribution.

! Basic statistics functionality
module stats
  use globvars

  implicit none

  private

  public :: stats_residuals
  public :: stats_mean
  public :: stats_variance

contains
  ! Calculate residuals
  subroutine stats_residuals(data_arr, theor_arr, resid_arr)
    real(dp), intent(in) :: data_arr(:), theor_arr(:)
    real(dp), intent(inout) :: resid_arr(:)

    resid_arr(:) = abs(data_arr(:) - theor_arr(:))

  end subroutine stats_residuals

  ! Calculate mean
  real(dp) function stats_mean(data_arr) result(val)
    real(dp), intent(in) :: data_arr(:)

    integer(dp) :: n

    n = size(data_arr)
    val = sum(data_arr) / n
  end function stats_mean

  ! Calculate variance
  real(dp) function stats_variance(data_arr) result(val)
    real(dp), intent(in) :: data_arr(:)

    real(dp) :: mean
    integer(dp) :: n

    mean = stats_mean(data_arr)
    n = size(data_arr)

    val = sum((data_arr(:) - mean)**2) / (n - 1)

  end function stats_variance
end module stats
