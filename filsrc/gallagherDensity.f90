
subroutine gallagher(x,y,z,density)
!GCPM v2.4 wrapper for filament code
   use parametersBackground, only: kpindex
   use constants, only: Pi
   implicit none
   real,intent(in) :: x,y,z
   real,intent(out) :: density
   real*4 r,amlt,alatr,outn(4)
   real*4 akp
   integer itime(2)
   !write(*,*) x,y,z
   r = real(sqrt(x**2+y**2+z**2),4)
   if(x .ne. 0.0)then
      amlt = real(atan(z/sqrt(x**2+y**2)),4)
   else
      amlt = real(sign(1.0,z)*Pi/2.0,4)
   endif

   if(x .ne. 0.0)then
      alatr = real((24.0/(2*Pi))*(atan(y/x)+Pi),4)
   else
      alatr = real((24.0/(2*Pi))*(sign(1.0,y)*Pi/2.0+Pi),4)
   endif

   akp = kpindex
   itime(1) = 2001093
   itime(2) = 0
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! c  Input parameters:
! c
! c  itime integer*4   dimensions=2
! c     (1) = yeardoy, e.g. 2001093
! c     (2) = miliseconds of day
! c  r     real*4      dimension=1
! c     geocentric radial distance in Re
! c  amlt  real*4      dimension=1
! c     solar magnetic local time in hours
! c  alatr real*4      dimension=1
! c     solar magnetic latitude in radians
! c  akp      real*4      dimension=1
! c     planetary Kp index
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! c  Output parameters:
! c  outn  real*4      dimensions=4
! c        (1) = total electron density in 1/cm^3
! c        (2) = total hydrogen density in 1/cm^3
! c        (3) = total helium density in 1/cm^3
! c        (4) = total oxygen density in 1/cm^3
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   call gcpm_v24(itime,r,amlt,alatr,akp,outn)
   !write(*,*) r,amlt,alatr,akp,outn(2)
   density = outn(2)
   return
end subroutine gallagher