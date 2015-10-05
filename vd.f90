! Copyright (c) 2015 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file at the top-level directory of this distribution.

! 2-dimensional virtual detector module
module vd
  use progvars
  use numerics, only: numerics_linspace, numerics_cmplx_phase, numerics_d1, &
       numerics_trapz

  implicit none

  private

  public :: vd_update
  public :: vd_normalize
  public :: vd_init
  public :: vd_cleanup

  ! Work arrays
  complex(dp), allocatable :: phi_arr(:,:), px_arr(:,:), py_arr(:,:)
  real(dp), allocatable :: mag_arr(:,:), jx_arr(:,:), jy_arr(:,:)

contains

  subroutine vd_init()

    logical :: sane

    ! Check internal / external / virtual detector grid limits
    sane = .true.

    if (nx .le. nxl_external + nxr_external) then
       sane = .false.
    else if (nxl_external .le. vd_nxl) then
       sane = .false.
    else if (nxr_external .le. vd_nxr) then
       sane = .false.
    end if

    if (.not. sane) then
       write(*,*) "Error in x-grid configuration, exiting."
       call exit(0)
    end if

    if (ny .le. nyl_external + nyr_external) then
       sane = .false.
    else if (nyl_external .le. vd_nyl) then
       sane = .false.
    else if (nyr_external .le. vd_nyr) then
       sane = .false.
    end if

    if (.not. sane) then
       write(*,*) "Error in y-grid configuration, exiting."
       call exit(0)
    end if

    ! Calculate edges of VD grid
    vd_xl_min = nxl_external - vd_nxl
    vd_xl_max = nxl_external
    vd_xr_min = nx - nxr_external
    vd_xr_max = nx - nxr_external + vd_nxr

    vd_yl_min = nyl_external - vd_nyl
    vd_yl_max = nyl_external
    vd_yr_min = ny - nyr_external
    vd_yr_max = ny - nyr_external + vd_nyr

    ! Allocate work arrays
    allocate(phi_arr(nx, ny))
    allocate(px_arr(nx, ny))
    allocate(py_arr(nx, ny))
    allocate(mag_arr(nx, ny))
    allocate(jx_arr(nx, ny))
    allocate(jy_arr(nx, ny))

  end subroutine vd_init

  subroutine vd_normalize(npx_arr, npy_arr)
    real(dp), intent(inout) :: npx_arr(:), npy_arr(:)

    real(dp) :: npx_norm, npy_norm

    npx_norm = 1.0_dp / numerics_trapz(npx_arr, vd_dpx)
    npy_norm = 1.0_dp / numerics_trapz(npy_arr, vd_dpy)

    npx_arr(:) = npx_norm * npx_arr(:)
    npy_arr(:) = npy_norm * npy_arr(:)
  end subroutine vd_normalize

  subroutine vd_update(psi_arr, npx_arr, npy_arr)
    complex(dp), intent(in) :: psi_arr(:,:)
    real(dp) :: npx_arr(:), npy_arr(:)

    call vd_fill_arrays(psi_arr)
    call vd_calc_pj()
    call vd_bin(npx_arr, npy_arr)

  end subroutine vd_update

  subroutine vd_bin(npx_arr, npy_arr)
    integer(dp) :: i_x, i_y

    real(dp), intent(inout) :: npx_arr(:), npy_arr(:)

    ! Go along the entirety of the y-edge (bottom to top)
    do i_y = vd_yl_min, vd_yr_max

       ! Accumulate counts along the left side
       do i_x = vd_xl_min, vd_xl_max
          call accumulate_counts(i_x, i_y, px_arr, jx_arr, npx_arr, &
               vd_px_arr, vd_dpx)
          call accumulate_counts(i_x, i_y, py_arr, jy_arr, npy_arr, &
               vd_py_arr, vd_dpy)
       end do

       ! Accumulate counts along the right side
       do i_x = vd_xr_min, vd_xr_max
          call accumulate_counts(i_x, i_y, px_arr, jx_arr, npx_arr, &
               vd_px_arr, vd_dpx)
          call accumulate_counts(i_x, i_y, py_arr, jy_arr, npy_arr, &
               vd_py_arr, vd_dpy)
       end do

    end do

    ! Go along the entirety of the x-edge (left to right)
    do i_x = vd_xl_min, vd_xr_max

       ! Accumulate counts along the bottom side
       do i_y = vd_yl_min, vd_yl_max
          call accumulate_counts(i_x, i_y, px_arr, jx_arr, npx_arr, &
               vd_px_arr, vd_dpx)
          call accumulate_counts(i_x, i_y, py_arr, jy_arr, npy_arr, &
               vd_py_arr, vd_dpy)
       end do

       ! Accumulate counts along the top side
       do i_y = vd_yr_min, vd_yr_max
          call accumulate_counts(i_x, i_y, px_arr, jx_arr, npx_arr, &
               vd_px_arr, vd_dpx)
          call accumulate_counts(i_x, i_y, py_arr, jy_arr, npy_arr, &
               vd_py_arr, vd_dpy)
       end do
    end do
!
  contains
    subroutine accumulate_counts(i_x, i_y, p_arr, j_arr, count_arr,&
         p_bin_arr, dp_bin)
      integer(dp), intent(in) :: i_x, i_y
      real(dp), intent(inout) :: count_arr(:)
      complex(dp), intent(in) :: p_arr(:,:)
      real(dp), intent(in) :: j_arr(:,:)
      real(dp), intent(in) :: p_bin_arr(:), dp_bin

      real(dp) :: p_bin_min
      integer(dp) :: i_p
      real(dp) :: p

      p = real(p_arr(i_x, i_y))

      ! Mid-point histogram
      p_bin_min = p_bin_arr(1) - dp_bin / 2

      ! Calculate bin index - much more efficient than brute force
      i_p = ceiling((p - p_bin_min) / dp_bin)

      if (i_p .ge. 1 .and. i_p .le. size(p_bin_arr)) then
         count_arr(i_p) = count_arr(i_p) + dp_bin * abs(j_arr(i_x, i_y)) * dt
      end if

    end subroutine accumulate_counts
  end subroutine vd_bin

  subroutine vd_calc_pj()

    integer(dp) :: i_x, i_y

    ! Iterate over each y index and calculate x-component of p, j
    do i_y = vd_yl_min - 1, vd_yl_max + 1
       call calc_along_x(vd_xl_min - 1, vd_xr_max + 1, i_y)
    end do
    do i_y = vd_yr_min - 1, vd_yr_max + 1
       call calc_along_x(vd_xl_min - 1, vd_xr_max + 1, i_y)
    end do
    do i_y = vd_yl_max + 1, vd_yr_min - 1
       call calc_along_x(vd_xl_min - 1, vd_xl_max + 1, i_y)
       call calc_along_x(vd_xr_min - 1, vd_xr_max + 1, i_y)
    end do

    ! Iterate over each x index and calculate y-component of p, j
    do i_x = vd_xl_min - 1, vd_xl_max + 1
       call calc_along_y(i_x, vd_yl_min - 1, vd_yr_max + 1)
    end do
    do i_x = vd_xr_min - 1, vd_xr_max + 1
       call calc_along_y(i_x, vd_yl_min - 1, vd_yr_max + 1)
    end do
    do i_x = vd_xl_max + 1, vd_xr_min - 1
       call calc_along_y(i_x, vd_yl_min - 1, vd_yl_max + 1)
       call calc_along_y(i_x, vd_yr_min - 1, vd_yr_max + 1)
    end do

  contains

    subroutine calc_along_x(i_x_min, i_x_max, i_y)
      integer(dp) :: i_x_min, i_x_max, i_y

      call numerics_d1(phi_arr(i_x_min:i_x_max, i_y), &
           px_arr(i_x_min:i_x_max, i_y), dx)
      jx_arr(i_x_min:i_x_max, i_y) = real(mag_arr(i_x_min:i_x_max, i_y) / m * &
           px_arr(i_x_min:i_x_max, i_y))
    end subroutine calc_along_x

    subroutine calc_along_y(i_x, i_y_min, i_y_max)
      integer(dp) :: i_x, i_y_min, i_y_max

      call numerics_d1(phi_arr(i_x, i_y_min:i_y_max), &
           py_arr(i_x, i_y_min:i_y_max), dy)
      jy_arr(i_x, i_y_min:i_y_max) = real(mag_arr(i_x, i_y_min:i_y_max) / m * &
           py_arr(i_x, i_y_min:i_y_max))
    end subroutine calc_along_y

  end subroutine vd_calc_pj

  subroutine vd_fill_arrays(psi_arr)

    complex(dp), intent(in) :: psi_arr(:,:)
    complex(dp) :: z
    integer(dp) :: i_x, i_y

    ! "left" end of VD grid
    do i_x = vd_xl_min - 1, vd_xl_max + 1

       do i_y = vd_yl_min - 1, vd_yr_max + 1
          call fill_by_index(i_x, i_y)
       end do

    end do

    ! "right" end of VD grid
    do i_x = vd_xr_min - 1, vd_xr_max + 1

       do i_y = vd_yl_min - 1, vd_yr_max + 1
          call fill_by_index(i_x, i_y)
       end do

    end do

    ! rest of VD grid
    do i_x = vd_xl_max + 1, vd_xr_min - 1

       ! unfilled portion of "bottom" end of VD grid
       do i_y = vd_yl_min - 1, vd_yl_max + 1
          call fill_by_index(i_x, i_y)
       end do

       ! unfilled portion of "top" end of VD grid
       do i_y = vd_yr_min - 1, vd_yr_max + 1
          call fill_by_index(i_x, i_y)
       end do

    end do

  contains

    subroutine fill_by_index(i_x, i_y)
      integer(dp) :: i_x, i_y

      z = psi_arr(i_x, i_y)

      ! multiply by hbar here since Psi ~ exp(i S / hbar)
      phi_arr(i_x, i_y) = numerics_cmplx_phase(z) * hbar

      mag_arr(i_x, i_y) = abs(z)**2
    end subroutine fill_by_index

  end subroutine vd_fill_arrays

  subroutine vd_cleanup()

    deallocate(phi_arr)
    deallocate(px_arr)
    deallocate(py_arr)
    deallocate(mag_arr)
    deallocate(jx_arr)
    deallocate(jy_arr)

  end subroutine vd_cleanup

end module vd
