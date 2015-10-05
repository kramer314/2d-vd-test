! Copyright (c) 2015 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file at the top-level directory of this distribution.

program vd_test
  use progvars
  use propagate, only: propagate_psi
  use setup, only: setup_init, setup_cleanup
  use vd, only: vd_update
  use gaussian, only: gaussian_xyt
  
  implicit none

  integer(dp) :: i_t
  integer(dp) :: i_kx

  call setup_init()

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

  do i_t = 1, nt
     call propagate_psi(psi_arr, i_t)

     call vd_update(psi_arr, nkx_arr, nky_arr)

     if (mod(i_t, print_mod_t) .eq. 0) then
        write(*,*) i_t, nt
     !   call output_vd_bins()
     end if
  end do

  
  open(98, file="vd_x.dat")
  open(99, file="vd_y.dat")

  do i_kx = 1, size(nkx_arr)
     write(98,*) vd_kx_arr(i_kx), nkx_arr(i_kx)
  end do
  do i_kx = 1, size(nky_arr)
     write(99, *) vd_ky_arr(i_kx), nky_arr(i_kx)
  end do

  call setup_cleanup()

end program vd_test
