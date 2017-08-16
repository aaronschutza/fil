subroutine filfop
! calculation of the force operator for the linearized EOM
! version 2 where divB is enforced 5/17 frt
use stateMod
use stateIndex
use constants
use dnams
implicit none
real :: Bx,Bz,pBack,dxBx,dxBz,dzBx,dzBz,ptotal,dxPtot,dzPtot
real,allocatable,dimension(:) :: s_Bx,s_Bz,s_pBack,s_dxBx,s_dxBz,s_dzBx
real,allocatable,dimension(:) :: s_dzBz,s_ptotal,s_dxPtot,s_dzPtot
real,allocatable,dimension(:) :: dbydy,dbds,d2bds2,d2b2ds2,kdotgradb
real,allocatable,dimension(:) :: dss ! ss array spacing
real,allocatable,dimension(:,:) :: khat,bhat
real,allocatable,dimension(:) :: kappa,ck,Bmag,Bc1,Bc2,Bc3,Bc4,cs,ca,fs,fa
real,allocatable,dimension(:) :: C1,C2,C3,C4,C5,C6,C7,C8,C9,C10
real,allocatable,dimension(:,:) :: fop
integer :: i,j
real :: ss,xx,zz
real :: psum

write(*,*) 'Calculating force operator'

allocate(khat(2,dimj))
allocate(bhat(2,dimj))

allocate(kappa(dimj))
allocate(Bmag(dimj))
allocate(ck(dimj))
allocate(cs(dimj))
allocate(ca(dimj))
allocate(fs(dimj))
allocate(fa(dimj))
allocate(Bc1(dimj))
allocate(Bc2(dimj))
allocate(Bc3(dimj))
allocate(Bc4(dimj))
allocate(dbydy(dimj))
allocate(kdotgradb(dimj))
allocate(dss(dimj))

allocate(C1(dimj))
allocate(C2(dimj))
allocate(C3(dimj))
allocate(C4(dimj))
allocate(C5(dimj))
allocate(C6(dimj))
allocate(C7(dimj))
allocate(C8(dimj))
allocate(C9(dimj))
allocate(C10(dimj))
allocate(fop(2*dimj,2*dimj))

allocate(s_Bx(dimj))
allocate(s_Bz(dimj))
allocate(s_pBack(dimj))
allocate(s_dxBx(dimj))
allocate(s_dxBz(dimj))
allocate(s_dzBx(dimj))
allocate(s_dzBz(dimj))
allocate(s_ptotal(dimj))
allocate(s_dxPtot(dimj))
allocate(s_dzPtot(dimj))
allocate(dbds(dimj))
allocate(d2b2ds2(dimj))
allocate(d2bds2(dimj))


do j=1,dimj
       s_Bx(j)=Bx(state(x,j),state(z,j),-1)
       s_Bz(j)=Bz(state(x,j),state(z,j),-1)
       s_pBack(j)=pBack(state(x,j),state(z,j),-1)
       s_dxBx(j)=dxBx(state(x,j),state(z,j),-1)
       s_dxBz(j)=dxBz(state(x,j),state(z,j),-1)
       s_dzBx(j)=dzBx(state(x,j),state(z,j),-1)
       s_dzBz(j)=dzBz(state(x,j),state(z,j),-1)
       s_ptotal(j)=ptotal(state(x,j),state(z,j),-1)
       s_dxPtot(j)=dxPtot(state(x,j),state(z,j),-1)
       s_dzPtot(j)=dzPtot(state(x,j),state(z,j),-1)
enddo
! test frt make constant along a field line - doesn't seem to do much
! psum = 0.0
!do j=1,dimj
!        psum = psum + state(p,j)
!end do
!state(p,:) = psum/(dimj)

do j=1,dimj
       Bmag(j) = sqrt(s_Bx(j)**2+s_Bz(j)**2)
enddo

do j=1,dimj
       bhat(x,j) = s_Bx(j)/Bmag(j)
       bhat(z,j) = s_Bz(j)/Bmag(j)
enddo
! compute grid spacing centered on grid j
do j=1,dimj
       if(j==1)then
              dss(j)=sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
       elseif(j==dimj)then
              dss(j)=sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2)
       else
              dss(j) = 0.5*(sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
              &           + sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2))
       endif
enddo

do j=1,dimj
       if(j==1)then
              ss=sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
              khat(x,j)=(bhat(x,j+1)-bhat(x,j))/ss
              khat(z,j)=(bhat(z,j+1)-bhat(z,j))/ss
       elseif(j==dimj)then
              ss=sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2)
              khat(x,j)=(bhat(x,j)-bhat(x,j-1))/ss
              khat(z,j)=(bhat(z,j)-bhat(z,j-1))/ss
       else
              ss = sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
              &  + sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
              khat(x,j)=(bhat(x,j+1)-bhat(x,j-1))/ss
              khat(z,j)=(bhat(z,j+1)-bhat(z,j-1))/ss
       endif
enddo
do j=1,dimj
       kappa(j) = sqrt(khat(x,j)**2+khat(z,j)**2)
       !normalize khat
       khat(x,j) = khat(x,j)/kappa(j)
       khat(z,j) = khat(z,j)/kappa(j)
enddo

! dbds
do j=1,dimj
       if(j==1)then
              ss=sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
              dbds(j) = (Bmag(j+1)-Bmag(j))/ss
       elseif(j==dimj)then
              ss=sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2)
              dbds(j) = (Bmag(j)-Bmag(j-1))/ss
       else
              ss = sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
              &  + sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
              dbds(j) = (Bmag(j+1)-Bmag(j-1))/ss
       endif
enddo
! kdotgradb 
do j=1,dimj
              kdotgradb(j) = (khat(x,j)*(bhat(x,j)*s_dxBx(j)+bhat(z,j)*s_dxBz(j))  &
              +               khat(z,j)*(bhat(x,j)*s_dzBx(j)+bhat(z,j)*s_dzBz(j)))
enddo
! d2b2ds2 & d2bds2
do j=1,dimj
       if(j==1)then
              d2b2ds2(j) = (2.0*Bmag(j)**2 -5*Bmag(j+1)**2 +4*Bmag(j+2)**2-Bmag(j+3)**2)/dss(j)**2
              d2bds2(j) = (2.0*Bmag(j) -5*Bmag(j+1) +4*Bmag(j+2)-Bmag(j+3))/dss(j)**2
       elseif(j==dimj)then
              d2b2ds2(j) = (2.0*Bmag(j)**2 -5*Bmag(j-1)**2 +4*Bmag(j-2)**2-Bmag(j-3)**2)/dss(j)**2
              d2bds2(j) = (2.0*Bmag(j) -5*Bmag(j-1) +4*Bmag(j-2)-Bmag(j-3))/dss(j)**2
       else
              ss = 0.5*(((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
                      + ((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2))
              d2b2ds2(j) = (Bmag(j+1)**2 -2.0*Bmag(j)**2 + Bmag(j-1)**2)/ss
              d2bds2(j) = (Bmag(j+1) -2.0*Bmag(j) + Bmag(j-1))/ss
       endif
enddo

! dbydy using divB
do j=1,dimj
       dbydy(j) = -(s_dxBx(j) + s_dzBz(j))
enddo

do j=1,dimj
       ck(j)= (1/Bmag(j))*(khat(x,j)*(khat(x,j)*s_dxBx(j)+khat(z,j)*s_dzBx(j)) &
              &                        +khat(z,j)*(khat(x,j)*s_dxBz(j)+khat(z,j)*s_dzBz(j)))
enddo

do j=1,dimj
       cs(j) = (cm/Re)*sqrt(gamma*state(p,j)/(nano*state(d,j)*mi*cmC))
       ca(j) = (cm/Re)*Bmag(j)/(nano*sqrt(state(d,j)*mi*mu0*cmC/nano))
       fs(j) = cs(j)**2/(cs(j)**2+ca(j)**2)
       fa(j) = 1.0 - fs(j)
enddo

do j=1,dimj
       Bc1(j) = Bmag(j)*cs(j)**2/(cs(j)**2+ca(j)**2)
enddo

do j=1,dimj
       if(j==1)then
              ss=sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
              Bc2(j) = (s_ptotal(j+1)-s_ptotal(j))/ss &
              &      * ((cm/Re)**2)*Bmag(j)/(nano*state(d,j)*mi*cmC*(cs(j)**2+ca(j)**2)) &
              &      - (Bmag(j+1)-Bmag(j))/ss
       elseif(j==dimj)then
              ss=sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2)
              Bc2(j) = (s_ptotal(j)-s_ptotal(j-1))/ss &
              &      * ((cm/Re)**2)*Bmag(j)/(nano*state(d,j)*mi*cmC*(cs(j)**2+ca(j)**2)) &
              &      - (Bmag(j)-Bmag(j-1))/ss
       else
              ss = sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
              &  + sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
              Bc2(j) = (s_ptotal(j+1)-s_ptotal(j-1))/ss &
              &      * ((cm/Re)**2)*Bmag(j)/(nano*state(d,j)*mi*cmC*(cs(j)**2+ca(j)**2)) &
              &      - (Bmag(j+1)-Bmag(j-1))/ss
       endif
enddo



do j=1,dimj
       Bc3(j) = (khat(x,j)*s_dxPtot(j)+khat(z,j)*s_dzPtot(j)) &
              & * ((cm/Re)**2)*Bmag(j)/(nano*state(d,j)*mi*cmC*(cs(j)**2+ca(j)**2)) &
              & - Bmag(j)*kappa(j)*(cs(j)**2/(cs(j)**2+ca(j)**2)) &
              & - khat(x,j)*(bhat(x,j)*s_dxBx(j)+bhat(z,j)*s_dxBz(j)) &
              & - khat(z,j)*(bhat(x,j)*s_dzBx(j)+bhat(z,j)*s_dzBz(j))
enddo

do j=1,dimj
       Bc4(j) = -ck(j)*Bmag(j)
enddo

do j=1,dimj
!       c1(j)=Bmag(j)*Bc1(j)
       c1(j)=Bmag(j)**2*fs(j)
enddo

do j=1,dimj
!       if(j==1)then
!              ss=sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
!              c2(j) = Bmag(j)*Bc2(j)+(Bc1(j+1)*Bmag(j+1)-Bc1(j)*Bmag(j))/ss
!       elseif(j==dimj)then
!              ss=sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2)
!              c2(j) = Bmag(j)*Bc2(j)+(Bc1(j)*Bmag(j)-Bc1(j-1)*Bmag(j-1))/ss
!       else
!              ss = sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
!              &  + sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
!              c2(j) = Bmag(j)*Bc2(j)+(Bc1(j+1)*Bmag(j+1)-Bc1(j-1)*Bmag(j-1))/ss
!       endif
        c2(j) = Bmag(j)*dbds(j)*fs(j)*(1.0-2.0*fa(j))
enddo


do j=1,dimj
!       if(j==1)then
!              ss=sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
!              c3(j) = (Bc2(j+1)*Bmag(j+1)-Bc2(j)*Bmag(j))/ss
!       elseif(j==dimj)then
!              ss=sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2)
!              c3(j) = (Bc2(j)*Bmag(j)-Bc2(j-1)*Bmag(j-1))/ss
!       else
!              ss = sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
!              &  + sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
!              c3(j) = (Bc2(j+1)*Bmag(j+1)-Bc2(j-1)*Bmag(j-1))/ss
!       endif
        c3(j) =fs(j)*(-0.5*d2b2ds2(j) + 2.0*fa(j)*dbds(j)**2)
enddo

do j=1,dimj
!       c4(j)= -(Bmag(j)**2)*kappa(j)+Bc3(j)*Bmag(j)+Bmag(j) &
!              &  *(khat(x,j)*(bhat(x,j)*s_dxBx(j)+bhat(z,j)*s_dxBz(j)) &
!              &  + khat(z,j)*(bhat(x,j)*s_dzBx(j)+bhat(z,j)*s_dzBz(j)))
        c4(j) = -2*kappa(j)*Bmag(j)**2*fs(j)
enddo

do j=1,dimj
       if(j==1)then
              ss=sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
!              c5(j) = -kappa(j)*Bc4(j)*Bmag(j)  &
!              &     + (Bc3(j+1)*Bmag(j+1)-Bc3(j)*Bmag(j))/ss          &
!              &  + Bc4(j)*(khat(x,j)*(bhat(x,j)*s_dxBx(j)+bhat(z,j)*s_dxBz(j))  &
!              &  +         khat(z,j)*(bhat(x,j)*s_dzBx(j)+bhat(z,j)*s_dzBz(j)))
              c5(j) = -Bmag(j)*(kdotgradb(j+1)-kdotgradb(j))/ss+(kappa(j+1)-kappa(j))/ss*Bmag(j)**2*(fa(j)-fs(j))+&
                      kappa(j)*dbds(j)*Bmag(j)*(1.0-4.0*fs(j)**2)+dbydy(j)*(kdotgradb(j)-kappa(j)*Bmag(j))
       elseif(j==dimj)then
              ss=sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2)
!              c5(j) = -kappa(j)*Bc4(j)*Bmag(j)  &
!              &     + (Bc3(j)*Bmag(j)-Bc3(j-1)*Bmag(j-1))/ss          &
!             &  + Bc4(j)*(khat(x,j)*(bhat(x,j)*s_dxBx(j)+bhat(z,j)*s_dxBz(j))  &
!             &  +         khat(z,j)*(bhat(x,j)*s_dzBx(j)+bhat(z,j)*s_dzBz(j)))
              c5(j) = -Bmag(j)*(kdotgradb(j)-kdotgradb(j-1))/ss+(kappa(j)-kappa(j-1))/ss*Bmag(j)**2*(fa(j)-fs(j))+&
                      kappa(j)*dbds(j)*Bmag(j)*(1.0-4.0*fs(j)**2)+dbydy(j)*(kdotgradb(j)-kappa(j)*Bmag(j))
       else
              ss = sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
              &  + sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
!             c5(j) = -kappa(j)*Bc4(j)*Bmag(j)  &
!             &     + (Bc3(j+1)*Bmag(j+1)-Bc3(j-1)*Bmag(j-1))/ss          &
!             &  + Bc4(j)*(khat(x,j)*(bhat(x,j)*s_dxBx(j)+bhat(z,j)*s_dxBz(j))  &
!             &  +         khat(z,j)*(bhat(x,j)*s_dzBx(j)+bhat(z,j)*s_dzBz(j)))
              c5(j) = -Bmag(j)*(kdotgradb(j+1)-kdotgradb(j-1))/ss+(kappa(j+1)-kappa(j-1))/ss*Bmag(j)**2*(fa(j)-fs(j))+&
                      kappa(j)*Bmag(j)*dbds(j)*(1.0-4.0*fs(j)**2)+dbydy(j)*(kdotgradb(j)-kappa(j)*Bmag(j))
       endif
enddo

do j=1,dimj
!       c6(j) = 2.0*Bc1(j)*Bmag(j)*kappa(j)
       c6(j) = 2.0*Bmag(j)**2*kappa(j)*fs(j)
enddo

do j=1,dimj
!       c7(j) = 2.0*Bc2(j)*Bmag(j)*kappa(j)
       c7(j) = -fs(j)*kappa(j)*2.0*Bmag(j)*dbds(j)
enddo

do j=1,dimj
       c8(j) = Bmag(j)**2
enddo

do j=1,dimj
!       if(j==1)then
!              ss=sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
!              c9(j) = Bmag(j)*(Bmag(j+1)-Bmag(j))/ss
!       elseif(j==dimj)then
!              ss=sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2)
!              c9(j) = Bmag(j)*(Bmag(j)-Bmag(j-1))/ss
!       else
!              ss = sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
!              &  + sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
!              c9(j) = Bmag(j)*(Bmag(j+1)-Bmag(j-1))/ss
              c9(j) = Bmag(j)*dbds(j)
!       endif
enddo

do j=1,dimj
       if(j==1)then
              ss=sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
!              c10(j) = 2.0*Bmag(j)*Bc3(j)*kappa(j)-Bc4(j)**2 &
!              &      + Bmag(j)*(Bc4(j+1)-Bc4(j))/ss
             c10(j) = 2.0*(Bmag(j)*kappa(j))**2*(fa(j)-fs(j)) -2.0*Bmag(j)*kappa(j)*kdotgradb(j) - &
                      dbds(j)**2+Bmag(j)*d2b2ds2(j)-2.0*dbydy(j)*dbds(j)+Bmag(j)*(dbydy(j+1)-dbydy(j))/ss  - &
                      dbydy(j)**2  
       elseif(j==dimj)then
              ss=sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2)
!              c10(j) = 2.0*Bmag(j)*Bc3(j)*kappa(j)-Bc4(j)**2 &
!              &      + Bmag(j)*(Bc4(j)-Bc4(j-1))/ss
             c10(j) = 2.0*(Bmag(j)*kappa(j))**2*(fa(j)-fs(j)) -2.0*Bmag(j)*kappa(j)*kdotgradb(j) - &
                      dbds(j)**2+Bmag(j)*d2b2ds2(j)-2.0*dbydy(j)*dbds(j)+Bmag(j)*(dbydy(j)-dbydy(j-1))/ss  - &
                      dbydy(j)**2  
       else
              ss = sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
              &  + sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
!              c10(j) = 2.0*Bmag(j)*Bc3(j)*kappa(j)-Bc4(j)**2 &
!              &      + Bmag(j)*(Bc4(j+1)-Bc4(j-1))/ss

             c10(j) = 2.0*(Bmag(j)*kappa(j))**2*(fa(j)-fs(j)) -2.0*Bmag(j)*kappa(j)*kdotgradb(j) - &
                      dbds(j)**2+Bmag(j)*d2bds2(j)-2.0*dbydy(j)*dbds(j)+Bmag(j)*(dbydy(j+1)-dbydy(j-1))/ss  - &
                      dbydy(j)**2                      
       endif
enddo

fop = 0.0

!upper left quadrant
do j=2,dimj-1
       ss = 0.5*sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
       &  + 0.5*sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
       fop(j,j+1) = c1(j)/(ss**2)+c2(j)/(2.0*ss)
       fop(j,j) = -2.0*c1(j)/(ss**2) +c3(j)
       fop(j,j-1) = c1(j)/(ss**2)-c2(j)/(2.0*ss)
enddo

!upper right quadrant
do j=2,dimj-1
       ss = 0.5*sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
       &  + 0.5*sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
       fop(j,j+1+dimj) = c4(j)/(2.0*ss)
       fop(j,j+dimj) = 0.0 +c5(j)
       fop(j,j-1+dimj) = -c4(j)/(2.0*ss)
enddo

!lower left quadrant
do j=2,dimj-1
       ss = 0.5*sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
       &  + 0.5*sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
       fop(j+dimj,j+1) =  c6(j)/(2.0*ss)
       fop(j+dimj,j) = 0.0 +c7(j)
       fop(j+dimj,j-1) = -c6(j)/(2.0*ss)
enddo

!lower right quadrant
do j=2,dimj-1
       ss = 0.5*sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
       &  + 0.5*sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
       fop(j+dimj,j+1+dimj) = c8(j)/(ss**2)+c9(j)/(2.0*ss)
       fop(j+dimj,j+dimj) = -2.0*c8(j)/(ss**2) +c10(j)
       fop(j+dimj,j-1+dimj) = c8(j)/(ss**2)-c9(j)/(2.0*ss)
enddo

!open(100,file=trim(datapath)//'.debug.dat',status='unknown',recl=500)

!open(101,file=trim(datapath)//'.debug1.dat',status='unknown',recl=500)

!open(102,file=trim(datapath)//'.debug2.dat',status='unknown',recl=500)

!open(103,file=trim(datapath)//'.debug3.dat',status='unknown',recl=500)

!write(6,*)'dimj =',dimj
!write(100,*)'j, state(x,j), state(z,j), state(p,j), state(d,j), s_Bx(j), s_Bz(j), s_ptotal(j), Bmag(j)'
!do j=2,dimj-1
!write(100,*)j,state(x,j),state(z,j), state(p,j), state(d,j), s_Bx(j), s_Bz(j), s_ptotal(j), Bmag(j)
!enddo
!write(101,*)'j,khat(x,j), khat(z,j), cs(j), ca(j), Bc1(j), Bc2(j), Bc3(j), Bc4(j)'
!do j=2,dimj-1
!write(101,*)j,khat(x,j),khat(z,j),cs(j),ca(j),Bc1(j),Bc2(j),Bc3(j),Bc4(j)
!enddo
!write(102,*)'j, c1(j), c2(j), c3(j), c4(j), c5(j), c6(j), c7(j), c8(j), c9(j), c10(j)'
!do j=2,dimj-1
!write(102,*)j,c1(j), c2(j), c3(j), c4(j), c5(j), c6(j), c7(j), c8(j), c9(j), c10(j)
!enddo


!write(103,*)'j, s_dxBx(j), s_dxBz(j), s_dzBx(j), s_dzBz(j), s_dxPtot(j), s_dzPtot(j) '
!do j=2,dimj-1
!write(103,*)j,s_dxBx(j), s_dxBz(j), s_dzBx(j), s_dzBz(j), s_dxPtot(j), s_dzPtot(j)
!enddo

!close(100)
!close(101)
!close(102)
!close(103)

write(6,*)' writing fop file..',datapath

open(unit=111, file=trim(datapath)//'.fop', status="replace")
!do j = 2,dimj-1
!       write(111,*) fop(2:dimj-1,j),fop(dimj+2:2*dimj-1,j)
!enddo
!do j = dimj+2,2*dimj-1
!       write(111,*) fop(2:dimj-1,j),fop(dimj+2:2*dimj-1,j)
!enddo
do j = 2,dimj-1!
       write(111,*) fop(j,2:dimj-1),fop(j,dimj+2:2*dimj-1)
enddo
do j = dimj+2,2*dimj-1
       write(111,*) fop(j,2:dimj-1),fop(j,dimj+2:2*dimj-1)
enddo

close(111)

deallocate(khat)
deallocate(bhat)

deallocate(kappa)
deallocate(Bmag)
deallocate(ck)
deallocate(cs)
deallocate(ca)
deallocate(fs)
deallocate(fa)
deallocate(Bc1)
deallocate(Bc2)
deallocate(Bc3)
deallocate(Bc4)

deallocate(C1)
deallocate(C2)
deallocate(C3)
deallocate(C4)
deallocate(C5)
deallocate(C6)
deallocate(C7)
deallocate(C8)
deallocate(C9)
deallocate(C10)
deallocate(fop)

deallocate(s_Bx)
deallocate(s_Bz)
deallocate(s_pBack)
deallocate(s_dxBx)
deallocate(s_dxBz)
deallocate(s_dzBx)
deallocate(s_dzBz)
deallocate(s_ptotal)
deallocate(s_dxPtot)
deallocate(s_dzPtot)
deallocate(dbydy)
deallocate(dbds)
deallocate(d2b2ds2)
deallocate(d2bds2)
deallocate(kdotgradb)
deallocate(dss)

end subroutine filfop
