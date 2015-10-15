! Copyright (c) 2015 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file at the top-level directory of this distribution.

! Setup module
module setup

  use progvars

  use params, only: params_init, params_cleanup
  use gaussian, only: gaussian_init, gaussian_cleanup
  use numerics, only: numerics_linspace
  use vd, only: vd_init, vd_cleanup

  implicit none

  private

  public :: setup_init
  public :: setup_cleanup

contains
  subroutine setup_init()
    call params_init()
    call gaussian_init()
    call vd_init()

    call allocate_arrays()
    call init_arrays()
  end subroutine setup_init

  subroutine setup_cleanup()
    call vd_cleanup()
    call gaussian_cleanup()
    call params_cleanup()

    call deallocate_arrays()
  end subroutine setup_cleanup

  subroutine init_arrays()
    ! initialize numerical grids
    call numerics_linspace(t_min, t_max, t_range, dt)
    call numerics_linspace(x_min, x_max, x_range, dx)
    call numerics_linspace(y_min, y_max, y_range, dy)

    ! initialize virtual detector grids
    call numerics_linspace(vd_px_min, vd_px_max, vd_px_arr, vd_dpx)
    call numerics_linspace(vd_py_min, vd_py_max, vd_py_arr, vd_dpy)

    ! initialize virtual detector counts
    npx_arr(:) = 0.0_dp
    npy_arr(:) = 0.0_dp
  end subroutine init_arrays

  subroutine allocate_arrays()
    allocate(psi_arr(nx, ny))
    allocate(t_range(nt))
    allocate(x_range(nx))
    allocate(y_range(ny))

    allocate(vd_px_arr(vd_npx))
    allocate(vd_py_arr(vd_npy))

    allocate(npx_arr(vd_npx))
    allocate(npy_arr(vd_npy))

    allocate(theor_npx_arr(vd_npx))
    allocate(theor_npy_arr(vd_npy))
    allocate(resid_npx_arr(vd_npx))
    allocate(resid_npy_arr(vd_npy))
  end subroutine allocate_arrays

  subroutine deallocate_arrays()
    deallocate(psi_arr)
    deallocate(t_range)
    deallocate(x_range)
    deallocate(y_range)

    deallocate(vd_px_arr)
    deallocate(vd_py_arr)

    deallocate(npx_arr)
    deallocate(npy_arr)

    deallocate(theor_npx_arr)
    deallocate(theor_npy_arr)
    deallocate(resid_npx_arr)
    deallocate(resid_npy_arr)
  end subroutine deallocate_arrays

end module setup
