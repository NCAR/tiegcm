module char_module
! some useful functions to find a string in a string array

  implicit none

  contains
!-----------------------------------------------------------------------
  pure function ismember(instance,group) result(flag)

    character(len=*),intent(in) :: instance
    character(len=*),dimension(:),intent(in) :: group
    logical :: flag

    integer :: i

    flag = .false.
    do i = 1,size(group)
      if (trim(instance) == trim(group(i))) then
        flag = .true.
        exit
      endif
    enddo

  end function ismember
!-----------------------------------------------------------------------
  pure function find_index(instance,group) result(idx)

    character(len=*),intent(in) :: instance
    character(len=*),dimension(:),intent(in) :: group
    integer :: idx

    do idx = 1,size(group)
      if (trim(instance) == trim(group(idx))) return
    enddo
    idx = 0

  end function find_index
!-----------------------------------------------------------------------
end module char_module
