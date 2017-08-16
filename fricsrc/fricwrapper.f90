
subroutine fric2d_wrapper(bxin,bzin,pin,xin,zin,nx,nz,tilt_anglein,bxo,bzo,po,dirnam)
! routine to convert marina data into friction code format data
! 12/14 frt
! modified as subroutine for linking with tfc
! 10/15 ams
use coredata, only: iprec,rprec
implicit none
integer(iprec) :: nx,nz
integer(iprec) :: i,k
integer(iprec) :: iin,kin,rec_1
integer(iprec),parameter :: in_tot =0
real*8 :: xin(nx),zin(nz),bxin(nx,nz),bzin(nx,nz),pin(nx,nz),tilt_anglein,bxo(nx,nz),bzo(nx,nz),po(nx,nz)
real(rprec), dimension(nx) :: x,rx
real(rprec), dimension(nz) :: z,rz
real(rprec),dimension(nx,nz) :: bx,bz,vx,vz,p,rho,kx,kz,fcp
real(rprec) :: t,zeit,timeut,dm,tilt_angle,vis,gamma,ekave0,ekave
real(rprec),parameter :: bfactor =1.0e-9 ! convert to T
real(rprec),parameter :: pfactor =1.0e-9 ! convert to P
character*40::header
integer(iprec) :: record_number,iter_number
real(rprec) :: time_fric
character(len=100) :: dirnam


tilt_angle = real(tilt_anglein,rprec)


dm = 3.0e4*bfactor


! set to zero
vx=0.0; vz =0.0; rho=0.0; kx=0.0; kz=0.0

do i=1,nx
 do k=1,nz
   bx(i,k) = real(bxin(i,k),rprec)*bfactor
   bz(i,k) = real(bzin(i,k),rprec)*bfactor
    p(i,k) = real(pin(i,k),rprec)*pfactor
 end do
end do

do i=1,nx
  x(i) = real(xin(i),rprec)
enddo
do i=1,nz
  z(i) = real(zin(i),rprec)
enddo

! now compute rx and rz
do i=2,nx-1
 rx(i) = 1.0/(x(i+1)-x(i-1))
end do
 rx(1) = rx(2)
 rx(nx) = rx(nx-1)
do k=2,nz-1
 rz(k) = 1.0/(z(k+1)-z(k-1))
end do
 rz(1) = rz(2)
 rz(nz) = rz(nz-1)



! now output
open(unit=11,file=trim(dirnam)//'/'//'fric2d.init',status='unknown',form='unformatted')

zeit = 0.0; zeit=0.0;timeut =0.0; dm = 0.0; tilt_angle = 0.0; vis =0.0
gamma = 1.667 ; ekave0 = 0.0; ekave = 0.0

write (11) in_tot,zeit, timeut, nx,nz,dm,tilt_angle,x,z,rx,rz, &                         
                                rho,vx,vz,bx,bz,p, &                                      
                                vis,gamma,kx,kz,ekave0,ekave 

close(11)

!call friction code here
call fric2d_sub(nx,nz,dirnam)


rec_1 = 2
OPEN (UNIT=10, FILE = trim(dirnam)//'/'//'fric2d.dat', STATUS='old', form='unformatted')
! Skip to the starting record:

DO record_number = 1, rec_1 - 1
 READ (10) iter_number, time_fric, timeut,&
           nx,  nz, &
           dm, tilt_angle,&
           x, z, rx, rz, &
           rho,vx,vz,bx,bz,p, &
           vis, gamma, &
           kx, kz, ekave0, ekave
END DO

! Begin reading in data:

DO record_number = rec_1, rec_1

 WRITE (*,'(//A,I6.6)') 'Reading record number=', record_number
 READ (10) iter_number, time_fric, timeut,&
           nx,  nz, &
           dm, tilt_angle,&
           x,  z, rx, rz, &
           rho,vx,vz,bx,bz,p, &
           vis, gamma, &
           kx,  kz, ekave0, ekave

 WRITE (*,'(A,I8.8,2(A,G9.2),3(A,I4.4))') 'ITER=', iter_number, &
        '   TIME_FRIC=', time_fric, '  UT=', timeut, &
        '   NX=',nx,'  NZ=', nz
END DO

close(10)
do i = 1,nx
 do k = 1,nz
  bxo(i,k) = real(bx(i,k),8)/bfactor
  bzo(i,k) = real(bz(i,k),8)/bfactor
  po(i,k) = real(p(i,k),8)/pfactor
 enddo
enddo

end subroutine

