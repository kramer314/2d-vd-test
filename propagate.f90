module propagate
  use progvars
  use gaussian, only: gaussian_xyt

  implicit none

  private

  public :: propagate_psi

contains

  subroutine propagate_psi(psi_arr, i_t)
    ! Only fill in the virtual detector grid

    complex(dp), intent(inout) :: psi_arr(:,:)
    integer(dp), intent(in) :: i_t

    real(dp) :: x, y, t
    integer(dp) :: i_x, i_y

    t = t_range(i_t)
    
    ! "left" end of VD grid
    do i_x = vd_xl_min, vd_xl_max
       x = x_range(i_x)

       do i_y = vd_yl_min, vd_yr_max
          y = y_range(i_y)
          psi_arr(i_x, i_y) = gaussian_xyt(x, y, t)
       end do

    end do

    ! "right" end of VD grid
    do i_x = vd_xr_min, vd_xr_max
       x = x_range(i_x)

       do i_y = vd_yl_min, vd_yr_max
          y = y_range(i_y)
          psi_arr(i_x, i_y) = gaussian_xyt(x, y, t)

       end do

    end do

    ! rest of VD grid
    do i_x = vd_xl_max, vd_xr_min
       x = x_range(i_x)

       ! unfilled portion of "bottom" end of VD grid
       do i_y = vd_yl_min, vd_yl_max
          y = y_range(i_y)
          psi_arr(i_x, i_y) = gaussian_xyt(x, y, t)
       end do

       ! unfilled portion of "top" end of VD grid
       do i_y = vd_yr_min, vd_yr_max
          y = y_range(i_y)
          psi_arr(i_x, i_y) = gaussian_xyt(x, y, t)
       end do

    end do

  end subroutine propagate_psi

end module propagate
