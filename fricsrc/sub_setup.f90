   module setupmod
   use coredata, only: iprec
   implicit none
   integer(iprec):: init_option

   character (LEN=100), PARAMETER :: setupinp = 'setup.inp'
   end module setupmod

   subroutine read_setup_inp ()
   use setupmod
   use coredata,only: symmetric
   use fricoptions
   implicit none

! magnetic field option
! 1. Pritchett 2010
! 2. harris
! 3. birn 1975

    init_option = initialization_option
    if(init_option==0)then
      write(6,*)' reading in data file '
    elseif(init_option==1)then
      write(6,*)' Pritchett initial condition option '
    elseif(init_option==2)then
      write(6,*)' Harris initial condition option '
    elseif(init_option==3)then
      write(6,*)' Birn 1975 initial condition option '
    else
      write(6,*)' undefined initial condition, stopping code'
      STOP
    end if      

    symmetric = NSsymmetric_option
    if(symmetric)then
            write(6,*)' N-S symmetry'
    else
            write(6,*)' NO N-S symmetry'
    end if
    return
    end subroutine read_setup_inp
!--------------------------------------------------

    subroutine setupd()
    use coredata
    use setupmod
    implicit none

    real(rprec) :: xmin,zmin,xmax,zmax

    real(rprec), PARAMETER :: timeut_local=-1.0
    integer(iprec) :: nx_neg,nx_pos
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(rprec) :: ssxp, rxp, ssyp, ryp, sszp, rzp
    real(rprec) :: dx, dy
    real(rprec) :: ximx
    real(rprec) :: rc, dz
    real(rprec) :: db1, db2, db3
    real(rprec) :: bb
    integer(iprec) ::  ix, iy, iz
    real(rprec) :: xstdoff
    real(rprec) mu0,alpha,apot,d
! finite difference computed gradients
    real(rprec) :: rxg(zahlx),rzg(zahlz)
    real(rprec) :: betx
    real(rprec) :: betz  ! , delz
    real(rprec) :: dmin
    real(rprec) :: fac1 
    integer(iprec) :: i,j,k
    integer(iprec) :: kstart
! read input file
    call read_setup_inp()
    if(init_option ==0)return ! will read in data file
! stuff not really used but initialized anyway 
    dm=30742.  
    tilt_angle=0.0
! setup the grid, see comments in the 3D friction code
    xstdoff = 12.0
    xmin = -40.0
    xmax = 0.0
    zmax =  1.25*xstdoff
    zmin =  0
    nx = zahlx ; nz = zahlz
     nx_pos = nx*xmax/(xmax-xmin)
!    nx_pos = nx*xmax/(xmax-xmin) + 4
    nx_neg = nx - nx_pos - 1

    write(6,*)' nx =',nx,' nx_neg =',nx_neg,' nx_pos=',nx_pos
    write(6,*)' nz =',nz
    betx = (xmax-3.35)/pi
    betx = 0.8*xmax/pi
    rc = .5*(nx_pos-1)
    dx = 1.0/(nx_pos-1)
      
   do ix=nx_neg+1,nx
      rx(ix) = rc/( xmax - pi*betx*cos(pi*(ix-nx_neg-1)*dx) )
      ssx(ix) = xmax*(ix-nx_neg-1)*dx - betx*sin(pi*(ix-nx_neg-1)*dx)
   end do
  
    ximx = (nx_neg-1.0)/(nx_pos-1.0)
    betx=(-xmin-3.35*ximx)/pi
      betx=(-xmin*0.8)/pi
   do ix=nx_neg+1,1,-1
      rx(ix) = -rc*ximx/( xmin + pi*betx*cos(pi*(nx_neg+1-ix)*dx/ximx) )
      ssx(ix) = xmin*(nx_neg+1-ix)*dx/ximx + betx*sin(pi*(nx_neg+1-ix)*dx/ximx)
   end do

! __________________________________________________________________

    betz = (zmax-3.0)/pi
    if(symmetric)then
            kstart = 2
            rc=.5*(nz-3)
            dz = 1.0/(nz-3)
    else
            if(mod(nz,2)==0)then
                 write(6,*)' nz must be even for non-symmetric grid'
                 stop
             end if             
             kstart = (nz+1)/2
             rc=.5*(kstart-3)
             dz = 1.0/(kstart-3)
    endif
    do iz=kstart,nz
       rz(iz) = rc/( zmax - pi*betz*cos(pi*(iz-kstart)*dz)) ! +2*pi*delz*cos(2*pi*(iz-2)*dz) )
       ssz(iz) = zmax*(iz-kstart)*dz - betz*sin(pi*(iz-kstart)*dz) ! +delz*sin(2*pi*(iz-2)*dz)
    end do
! reconstruct the bottom part of the grid    
    if(.not.symmetric)then
            do iz=1,kstart -1
                ssz(iz) = -ssz(nz-iz+1)
                rz(iz) = rz(nz-iz+1)
            end do
    end if

    write(*,*) 'ssx =',ssx
    write(*,*) 'ssz =',ssz
! numerical grid spacing calculation

    rxg(1) = rx(1)
    rxg(nx) = rx(nx)
    do ix=2,nx-1
        rxg(ix) = 1.0/(ssx(ix+1)-ssx(ix-1))
    end do

    rzg(1) = rz(1)
    rzg(nz) = rz(nz)
    do iz=2,nz-1
        rzg(iz) = 1.0/(ssz(iz+1)-ssz(iz-1))
    end do   
! grid setup finished
! setup initial condition

open(unit=11,file='setup_gridz.dat',status='unknown',form='formatted')
     do iz=1,nz
     write(11,*)iz,ssz(iz),rz(iz),rzg(iz)
     end do
     close(11)
  
     mu0 = 4e-7*pi
    if( init_option  == 1) then ! pritchett
        do k=1,zahlz
        do i=1,zahlx
        call  Pritchett_B(ssx(i),ssz(k),bx(i,k),bz(i,k),p(i,k)) 
!       write(*,*) 'x,z,bxz,p',ssx(i),ssz(k),bx(i,k),bz(i,k),p(i,k)
        end do 
        end do
    else if ( init_option  == 2) then  ! harris
        do k=1,zahlz
        do i=1,zahlx
         bx(i,k) = -40000.e-9*tanh(ssz(k)/3.)
         bz(i,k) = 0.
         p(i,k) = 6.3662e-4*(cosh(ssz(k)/3.))**(-2)    !nT
!        if((ssx(i)**2+ssz(k)**2)< 1) then
!        rho(i,k)=1e15
!        endif
!        if((ssx(i)**2+ssz(k)**2)< 2) then
!        p(i,k)=0.0
!       endif
        end do
        end do
     else if( init_option ==3 ) then  ! Birn
     d= 20
     alpha = 0.01
     do k=1,zahlz
       do i=1,zahlx
         fac1=1.e-9*dm
         apot=cos(pi*ssz(k)/(2*d))*exp(alpha*ssx(i))*fac1
         bx(i,k) =fac1*(pi/(2*d))*sin(pi*ssz(k)/(2*d))*exp(alpha*ssx(i)) 
         bz(i,k) =  alpha*apot
         p(i,k) = 1./(2*mu0)*((pi/(2*d))**2-alpha**2)*apot**2 
          
         if((ssx(i)**2+ssz(k)**2)< 4) then
         rho(i,k)=1.e15
         endif
      end do
     end do
         
    endif  
    
    do k=1,zahlz
       do i=1,zahlx
        if((ssx(i)**2+ssz(k)**2)<4.0) then
        rho(i,k)=-1.0
        else
        rho(i,k)=1.0
        endif 
   enddo
    enddo
    
    kx=0.0;kz=0.0
    vx=0.0;vz=0.0
    open(41,file='fric2d.init',form='UNFORMATTED',status='REPLACE')
    write(41) 0,0.0,timeut_local,nx,nz,dm,tilt_angle,ssx,ssz,rxg,rzg,&
    rho,vx,vz,bx,bz,p,vis,gamma,kx,kz,1.0,1.0

    close(41)


    end subroutine setupd


    subroutine Pritchett_B(x,z,bx,bz,p)
    use coredata, only: rprec,iprec,pi
! This is base on Pritchett and Coroniti 2010 JGR paper
    real(rprec):: x,z,bx,bz,p,x1,x2,x3,x4
    real(rprec):: epsmin,eps1,eps2
    real(rprec):: f,minusdfoverfL
    real(rprec):: n,kT,b0,n0,L


! The parameters are
       x1=-5.   ! the start of the Bz minimum
       x2=-15.  ! the end of the Bz minumum
       x3=-35.  ! the start of high Bz region
       x4=-40.  ! the right-side simulation boundary
       epsmin = 0.1 !Bz/Bo for the minumum region 
       eps1  = 0.1  !Bz/Bo for the left-side boundary-epsmin
       eps2  = 0.3  !Bz/Bo for the right-side boundary-epsmin
       L=4.   ! scale in z direction for current sheet half-thickness
       b0= 25.e-9  ! in data file B is in si
       kT=5000*1.6e-19   ! 5keV
!      n0=0.31042127*10*6. ! per m^3
       n0=0.31042127*1.0e6 ! per m^3
!  n0= b0(nT)**2/(8*1.6*pi*T(keV)) to satisfy the initial equalibrium. 
     if(x.ge.x1.and.x.lt.0.00) then
! note the grid in the paper is backward
      f= exp(1./L*(epsmin*(x-x1)-0.5*eps1*(x-x1)**2/x1))
      minusdfoverfL = epsmin+eps1*(1-x/x1)
     else if (x.ge.x2.and.x.lt.x1) then
      f= exp(1./L*(epsmin*(x-x1)))
      minusdfoverfL = epsmin     
     else if (x.ge.x3.and.x.lt.x2) then
      f= exp(1./L*(epsmin*(x-x1)+0.5*eps2*(x-x2)**2/(x3-x2)))
      minusdfoverfL = epsmin+eps2*(x-x2)/(x3-x2)
     else if (x.ge.x4.and.x.lt.x3) then
      f= exp(1./L*(epsmin*(x-x1)+0.5*eps2*(2*x-x2-x3)))
      minusdfoverfL = epsmin+ eps2
      else
       f =0.0 ! default set to zero
     endif
     bx = b0*f*tanh(f*z/L)
!    reversed sign because of backward grid for bx
!    what was plotted in figure 1 is actually -L*dfoverf
     bz = b0*minusdfoverfL- b0*z*tanh(f*z/L)*minusdfoverfL*f/L
     n=n0*f**2*(1/cosh(f*z/L)**2)
     p=n*kT


     end subroutine  pritchett_b
     
