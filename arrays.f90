module arrays
  use globvars
  implicit none

  private

  public :: arrays_quicksort

contains

  subroutine arrays_swap_idex(arr, i, j)
    ! Swap values of arr(i) and arr(j)

    real(dp), intent(inout) :: arr(:)
      integer, intent(in) :: i, j

      real(dp) :: temp

      temp = arr(i)
      arr(i) = arr(j)
      arr(j) = temp

  end subroutine arrays_swap_idex

  recursive subroutine arrays_quicksort(arr, i_min, i_max)
    ! Quicksort with Lomuto partitioning

    real(dp), intent(inout) :: arr(:)
    integer, intent(in) :: i_min, i_max

    integer :: part

    if (i_min < i_max) then
       call partition(arr, i_min, i_max, part)
       call arrays_quicksort(arr, i_min, part - 1)
       call arrays_quicksort(arr, part + 1, i_max)
    end if

  contains

    subroutine partition(arr, i_min, i_max, part)
      ! Lomuto partition scheme for quicksort

      real(dp), intent(inout) :: arr(:)
      integer, intent(in) :: i_min, i_max
      integer, intent(out) :: part

      integer :: i, j
      real(dp) :: pivot

      pivot = arr(i_max)

      i = i_min

      do j = i_min, i_max - 1
         if (arr(j) .le. pivot) then
            call arrays_swap_idex(arr, i, j)
            i = i + 1
         end if
      end do
      call arrays_swap_idex(arr, i, i_max)

      part = i
    end subroutine partition

  end subroutine arrays_quicksort
end module arrays
