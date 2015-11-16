! Copyright (c) 2015 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file at the top-level directory of this distribution.

! Assorted file operations
! Note: This is a program-independent module
module files
  implicit none

  private

  public :: files_ensure_dir

contains

  ! Ensure that a directory exists on a UNIX system
  !
  ! dir :: directory location
  subroutine files_ensure_dir(dir)

    character(*), intent(in) :: dir

    call system("mkdir -p "//trim(dir))

  end subroutine files_ensure_dir

end module files
