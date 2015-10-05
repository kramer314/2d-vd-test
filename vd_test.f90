! Copyright (c) 2015 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file at the top-level directory of this distribution.

program vd_test
  use progvars
  use propagate, only: propagate_psi
  use setup, only: setup_init, setup_cleanup
  use vd, only: vd_update, vd_normalize
  use gaussian, only: gaussian_xyt, gaussian_px, gaussian_py

  implicit none

  integer(dp) :: i_t
  integer(dp) :: i_px, i_py
  real(dp) :: px, py

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

     call vd_update(psi_arr, npx_arr, npy_arr)

     if (mod(i_t, print_mod_t) .eq. 0) then
        write(*,*) i_t, nt
     !   call output_vd_bins()
     end if
  end do

  call vd_normalize(npx_arr, npy_arr)


  open(98, file="vd_x.dat")
  open(99, file="vd_y.dat")

  do i_px = 1, size(npx_arr)
     px = vd_px_arr(i_px)
     write(98,*) px, npx_arr(i_px), abs(gaussian_px(px))**2
  end do
  do i_py = 1, size(npy_arr)
     py = vd_py_arr(i_py)
     write(99, *) py, npy_arr(i_py), abs(gaussian_py(py))**2
  end do

  call setup_cleanup()

end program vd_test
