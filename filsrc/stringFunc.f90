!************** Functions for handling strings ***************

character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str

character(len=20) function realtostr(r)
!   "Convert an real to string."
    real, intent(in) :: r
    write (realtostr,'(f20.3)') r
end function realtostr

subroutine intStrToLStr(lintStr,iniStr)
	implicit none
	integer :: lint
	character(len=100) :: iniStr,lintStr
	if (trim(adjustl(lintStr)) == 'T' .or. trim(adjustl(lintStr)) == 't'  &
		& .or. trim(adjustl(lintStr)) == 'True' .or. trim(adjustl(lintStr)) == 'true') then
		iniStr = 'T'
		return
	elseif (trim(adjustl(lintStr)) == 'F' .or. trim(adjustl(lintStr)) == 'f' &
		& .or. trim(adjustl(lintStr)) == 'False' .or. trim(adjustl(lintStr)) == 'false') then
		iniStr = 'F'
		return
	endif
	read(lintStr,*) lint
	if (lint == 0) then
		iniStr = 'F'
		return
	elseif (lint == 1) then
		iniStr = 'T'
		return
	else
		iniStr = 'F'
		return
	endif
	iniStr = 'F'
	return
end subroutine