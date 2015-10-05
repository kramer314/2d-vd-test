module propagate
  use progvars
  use gaussian, only: gaussian_xyt

  implicit none

  private

  public :: propagate_psi

contains

  subroutine propagate_psi(psi_arr, i_t)
    ! For testing, only fill in the VD grid and adjacent grid points (for
    ! derivatives)

    complex(dp), intent(inout) :: psi_arr(:,:)
    integer(dp), intent(in) :: i_t

    real(dp) :: x, y, t
    integer(dp) :: i_x, i_y

    t = t_range(i_t)

    ! "left" end of VD grid
    do i_x = vd_xl_min - 1, vd_xl_max + 1
       x = x_range(i_x)

       do i_y = vd_yl_min - 1, vd_yr_max + 1
          y = y_range(i_y)
          psi_arr(i_x, i_y) = gaussian_xyt(x, y, t)
       end do

    end do

    ! "right" end of VD grid
    do i_x = vd_xr_min - 1, vd_xr_max + 1
       x = x_range(i_x)

       do i_y = vd_yl_min - 1, vd_yr_max + 1
          y = y_range(i_y)
          psi_arr(i_x, i_y) = gaussian_xyt(x, y, t)

       end do

    end do

    ! rest of VD grid
    do i_x = vd_xl_max + 1, vd_xr_min - 1
       x = x_range(i_x)

       ! unfilled portion of "bottom" end of VD grid
       do i_y = vd_yl_min - 1, vd_yl_max + 1
          y = y_range(i_y)
          psi_arr(i_x, i_y) = gaussian_xyt(x, y, t)
       end do

       ! unfilled portion of "top" end of VD grid
       do i_y = vd_yr_min - 1, vd_yr_max + 1
          y = y_range(i_y)
          psi_arr(i_x, i_y) = gaussian_xyt(x, y, t)
       end do

    end do

  end subroutine propagate_psi

end module propagate
