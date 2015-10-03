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
    call numerics_linspace(vd_kx_min, vd_kx_max, vd_kx_arr, vd_dkx)
    call numerics_linspace(vd_ky_min, vd_ky_max, vd_ky_arr, vd_dky)
    
    ! initialize virtual detector counts
    nkx_arr = 0.0_dp
    nky_arr = 0.0_dp
  end subroutine init_arrays
  
  subroutine allocate_arrays()
    allocate(psi_arr(nx, ny))
    allocate(t_range(nt))
    allocate(x_range(nx))
    allocate(y_range(ny))

    allocate(vd_kx_arr(vd_nkx))
    allocate(vd_ky_arr(vd_nky))
    
    allocate(nkx_arr(vd_nkx))
    allocate(nky_arr(vd_nky))
  end subroutine allocate_arrays

  subroutine deallocate_arrays()
    deallocate(psi_arr)
    deallocate(t_range)
    deallocate(x_range)
    deallocate(y_range)

    deallocate(vd_kx_arr)
    deallocate(vd_ky_arr)
    
    deallocate(nkx_arr)
    deallocate(nky_arr)
  end subroutine deallocate_arrays
  
end module setup
