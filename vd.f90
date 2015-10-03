module vd
  use progvars
  use numerics, only: numerics_linspace, numerics_cmplx_phase, numerics_d1

  ! TODO: Actually calculate from 1 outside of VD grid (for derivatives, etc..)
  ! TODO: Better binning procedure
  
  implicit none

  ! At each time, for each x, compute k = del(phi), and j = A^2 / mu * k
  !   bin k - iterate over k_arr and check whether it fits in the bin;
  !   if so, increment bin count * j

  private

  public :: vd_update
  public :: vd_init
  public :: vd_cleanup

  ! Work arrays

  complex(dp), allocatable :: phi_arr(:,:), k_arr(:,:)
  real(dp), allocatable :: mag_arr(:,:), j_arr(:,:)

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
    allocate(k_arr(nx, ny))
    allocate(mag_arr(nx, ny))
    allocate(j_arr(nx, ny))

  end subroutine vd_init

  subroutine vd_update(psi_arr, nkx_arr, nky_arr)
    complex(dp), intent(in) :: psi_arr(:,:)
    real(dp) :: nkx_arr(:), nky_arr(:)

    call vd_fill_arrays(psi_arr)
    call vd_calc_kj()
    call vd_bin(nkx_arr, nky_arr)

  end subroutine vd_update

  subroutine vd_bin(nkx_arr, nky_arr)
    integer(dp) :: i_x, i_y

    real(dp), intent(inout) :: nkx_arr(:), nky_arr(:)

    ! Go along the entirety of the y-edge (bottom to top)
    do i_y = vd_yl_min, vd_yr_max

       ! Accumulate counts along the left side
       do i_x = vd_xl_min, vd_xl_max
          call accumulate_counts(i_x, i_y, nkx_arr, vd_kx_arr, vd_dkx)
       end do

       ! Accumulate counts along the right side
       do i_x = vd_xr_min, vd_xr_max
          call accumulate_counts(i_x, i_y, nkx_arr, vd_kx_arr, vd_dkx)
       end do

    end do

    ! Go along the entirety of the x-edge (left to right)
    do i_x = vd_xl_min, vd_xr_max

       ! Accumulate counts along the bottom side
       do i_y = vd_yl_min, vd_yl_max
          call accumulate_counts(i_x, i_y, nky_arr, vd_ky_arr, vd_dky)
       end do

       ! Accumulate counts along the top side
       do i_y = vd_yr_min, vd_yr_max
          call accumulate_counts(i_x, i_y, nky_arr, vd_ky_arr, vd_dky)
       end do
    end do

  contains
    subroutine accumulate_counts(i_x, i_y, count_arr, k_bin_arr, dk_bin)
      integer(dp), intent(in) :: i_x, i_y
      real(dp), intent(inout) :: count_arr(:)
      real(dp), intent(in) :: k_bin_arr(:), dk_bin

      real(dp) :: k_bin_min
      integer(dp) :: i_k
      real(dp) :: k

      k = real(k_arr(i_x, i_y))

      ! Mid-point histogram
      k_bin_min = k_bin_arr(1) - dk_bin / 2
      
      i_k = ceiling((k - k_bin_min) / dk_bin)

      if (i_k .ge. 1 .and. i_k .le. size(k_bin_arr)) then
         count_arr(i_k) = count_arr(i_k) + dk_bin * abs(j_arr(i_x, i_y)) * dt
      end if

      ! Brute-force approach
      !do i_k = 1, size(k_bin_arr) - 1
      !  if (k .ge. k_bin_arr(i_k) .and. k .lt. k_bin_arr(i_k + 1)) then
      !      count_arr(i_k) = count_arr(i_k) + dk_bin * abs(j_arr(i_x, i_y)) * dt
      !   end if
      !end do
    end subroutine accumulate_counts
  end subroutine vd_bin

  subroutine vd_calc_kj()

    integer(dp) :: i_x, i_y

    ! Go along the entirety of the y-edge (bottom to top)
    do i_y = vd_yl_min, vd_yr_max

       ! Calculate k,j along the left side
       call calc_along_x(vd_xl_min, vd_xl_max, i_y)
       ! Calculate k,j along the right side
       call calc_along_x(vd_xr_min, vd_xr_max, i_y)

    end do

    ! Go along the entirety of the x-edge (left to right)
    do i_x = vd_xl_min, vd_xr_max

       ! Calculate k,j along the bottom side
       call calc_along_y(i_x, vd_yl_min, vd_yl_max)
       ! Calculate k,j along the top side
       call calc_along_y(i_x, vd_yr_min, vd_yr_max)

    end do

  contains

    subroutine calc_along_x(i_x_min, i_x_max, i_y)
      integer(dp) :: i_x_min, i_x_max, i_y

      call numerics_d1(phi_arr(i_x_min:i_x_max, i_y), &
           k_arr(i_x_min:i_x_max, i_y), dx)
      j_arr(i_x_min:i_x_max, i_y) = real(mag_arr(i_x_min:i_x_max, i_y) / m * &
           k_arr(i_x_min:i_x_max, i_y))
    end subroutine calc_along_x

    subroutine calc_along_y(i_x, i_y_min, i_y_max)
      integer(dp) :: i_x, i_y_min, i_y_max

      call numerics_d1(phi_arr(i_x, i_y_min:i_y_max), &
           k_arr(i_x, i_y_min:i_y_max), dy)
      j_arr(i_x, i_y_min:i_y_max) = real(mag_arr(i_x, i_y_min:i_y_max) / m * &
           k_arr(i_x, i_y_min:i_y_max))
    end subroutine calc_along_y

  end subroutine vd_calc_kj

  subroutine vd_fill_arrays(psi_arr)

    complex(dp), intent(in) :: psi_arr(:,:)
    complex(dp) :: z
    integer(dp) :: i_x, i_y

    ! "left" end of VD grid
    do i_x = vd_xl_min, vd_xl_max

       do i_y = vd_yl_min, vd_yr_max
          call fill_by_index(i_x, i_y)
       end do

    end do

    ! "right" end of VD grid
    do i_x = vd_xr_min, vd_xr_max

       do i_y = vd_yl_min, vd_yr_max
          call fill_by_index(i_x, i_y)
       end do

    end do

    ! rest of VD grid
    do i_x = vd_xl_max, vd_xr_min

       ! unfilled portion of "bottom" end of VD grid
       do i_y = vd_yl_min, vd_yl_max
          call fill_by_index(i_x, i_y)
       end do

       ! unfilled portion of "top" end of VD grid
       do i_y = vd_yr_min, vd_yr_max
          call fill_by_index(i_x, i_y)
       end do

    end do

    ! Psi ~ exp(i S / hbar)
    phi_arr = hbar * phi_arr

  contains

    subroutine fill_by_index(i_x, i_y)
      integer(dp) :: i_x, i_y

      z = psi_arr(i_x, i_y)
      phi_arr(i_x, i_y) = numerics_cmplx_phase(z)
      mag_arr(i_x, i_y) = abs(z)**2
    end subroutine fill_by_index

  end subroutine vd_fill_arrays

  subroutine vd_cleanup()

    deallocate(phi_arr)
    deallocate(k_arr)
    deallocate(mag_arr)
    deallocate(j_arr)

  end subroutine vd_cleanup

end module vd
