subroutine fric2d_sub(nx,nz,name)

USE coredata, ONLY : iprec, rprec, read_data, write_data, dirname
IMPLICIT NONE
integer(iprec) :: ierr      = 0
integer(iprec) :: rec       = 1
integer(iprec) :: istart = 0
integer(iprec) :: nx,nz
logical :: relaxing = .TRUE.        ! relax or not?
character (LEN=100) :: fricinit = 'fric2d.init'
character (LEN=100) :: fricdata = 'fric2d.dat'
character (LEN=100) :: fricinp = 'fric.inp'
character (LEN=100) :: name
character (LEN=100) :: word
dirname = name
write(6,'(/A/)') 'Setup'

call allocate_core_data(nx,nz)

CALL Setupd()

word = trim(dirname)//'/'//trim(fricinit)

call read_data(word, 1, ierr)  ! get initial conditions

word = trim(dirname)//'/'//trim(fricdata)

call write_data(word, 1, ierr) 

istart = 1

if (ierr < 0) STOP 'unable to read fricinit with rec=1, aborting'
write(6,'(/A/)') '--------------------- Computing Force Balance --------------------------'
word = trim(dirname)//'/'//trim(fricdata)
call relaxdc(trim(dirname)//'/'//fricdata, istart, word, 0, 1, ierr)
if(ierr < 0) stop 'RELAXDC returned non-zero error flag, ABORTING'

call deallocate_core_data

end subroutine