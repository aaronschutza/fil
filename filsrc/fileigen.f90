subroutine fileigen
!eigenvalue and eigenvector solver for the linearized EOM
! version 2 where divB is enforced 5/17 frt
! version 3 option 11/17 frt
use stateMod
use stateIndex
use constants
use dnams
use initCond, only: boundary
use boundmod, only: iopen
use parametersBackground, only: useGrid,wallBound
implicit none
real :: Bx,Bz,pBack,dxBx,dxBz,dzBx,dzBz,ptotal,dxPtot,dzPtot,dxp,dzp
real,allocatable,dimension(:) :: s_Bx,s_Bz,s_pBack,s_dxBx,s_dxBz,s_dzBx
real,allocatable,dimension(:) :: s_dzBz,s_ptotal,s_dxPtot,s_dzPtot,s_dxP,s_dzP
real,allocatable,dimension(:) :: dbydy,dbds,d2bds2,d2b2ds2,kdotgradb,db2ds,kdotgradp,dkappads,dkdotgradpds,dbkdotgradbds
real,allocatable,dimension(:) :: dss ! ss array spacing
!real,allocatable,dimension(:) :: s0,s1,x1,z1 ! for variable grid spacingA
real,allocatable,dimension(:) :: s0
real,allocatable,dimension(:,:) :: khat,bhat,kvec
real,allocatable,dimension(:) :: kappa,ck,Bmag,Bc1,Bc2,Bc3,Bc4,cs,ca,fs,fa
real,allocatable,dimension(:) :: C,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10
real,allocatable,dimension(:) :: f3s,g3s ! for 3rd approach
real,allocatable,dimension(:,:) :: fop
real,allocatable,dimension(:) :: dsdi,d2sdi2
integer :: i,j
real :: ss,xx,zz
real :: psum
real :: interp1d
real :: alpha,delta,geom
!integer :: jdim
integer :: jmid,idir,isopen,jend
integer,parameter :: method=1
integer,parameter :: c5method=4
!real, parameter :: fac = 0.0 !80 ! grid stretch parameter < 1
real :: a0
!integer,parameter :: iopen=1 ! open bc
real :: x0,y0,z0,dx,dz,ds,eps,xe,xequator,eps0
!real :: dsdi, d2sdi2

write(*,*) 'Calculating force operator'
write(*,*)' C5 method =',c5method

allocate(khat(2,dimj))
allocate(kvec(2,dimj))
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
allocate(kdotgradp(dimj))
allocate(dss(dimj))
allocate(s0(dimj))
!allocate(s1(dimj))
!allocate(x1(dimj))
!allocate(z1(dimj))

allocate(C1(dimj))
allocate(C(dimj))
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
allocate(f3s(dimj))
allocate(g3s(dimj))
allocate(dsdi(dimj))
allocate(d2sdi2(dimj))

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
allocate(s_dxP(dimj))
allocate(s_dzP(dimj))
allocate(dbds(dimj))
allocate(db2ds(dimj))
allocate(d2b2ds2(dimj))
allocate(d2bds2(dimj))
allocate(dkappads(dimj))
allocate(dkdotgradpds(dimj))
allocate(dbkdotgradbds(dimj))

      if(iopen == 0)then
              write(6,*)' Closed BC for eigenmode analysis'
      else
              write(6,*)' Open BC for eigenmode analysis'
      end if

!--------------------------------------------
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
       s_dxP(j)=dxP(state(x,j),state(z,j),-1)
       s_dzP(j)=dzP(state(x,j),state(z,j),-1)
enddo
s0(1) = 0.0
do j=2,dimj
  s0(j) = s0(j-1) + sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2)
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
              kvec(x,j)=(bhat(x,j+1)-bhat(x,j))/ss
              kvec(z,j)=(bhat(z,j+1)-bhat(z,j))/ss
       elseif(j==dimj)then
              ss=sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2)
              kvec(x,j)=(bhat(x,j)-bhat(x,j-1))/ss
              kvec(z,j)=(bhat(z,j)-bhat(z,j-1))/ss
       else
              ss = sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
              &  + sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
              kvec(x,j)=(bhat(x,j+1)-bhat(x,j-1))/ss
              kvec(z,j)=(bhat(z,j+1)-bhat(z,j-1))/ss
       endif
enddo
do j=1,dimj
       khat(x,j) =  bhat(z,j)
       khat(z,j) = -bhat(x,j)
       kappa(j) = khat(x,j)*kvec(x,j)+khat(z,j)*kvec(z,j)
       !normalize khat
!       khat(x,j) = khat(x,j)/kappa(j)
!       khat(z,j) = khat(z,j)/kappa(j)
enddo
! kdotgradb  - in this code, mu0 = 4pie-7*1.0e9
do j=1,dimj
              kdotgradb(j) = (khat(x,j)*(bhat(x,j)*s_dxBx(j)+bhat(z,j)*s_dxBz(j))  &
              +               khat(z,j)*(bhat(x,j)*s_dzBx(j)+bhat(z,j)*s_dzBz(j)))
              kdotgradp(j) = (khat(x,j)*s_dxp(j)+khat(z,j)*s_dzp(j)) ! units!
enddo

! dbds
do j=1,dimj
       if(j==1)then
              ss=sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
              dbds(j) = (Bmag(j+1)-Bmag(j))/ss
              db2ds(j) = (Bmag(j+1)**2-Bmag(j)**2)/ss
              dkappads(j) = (kappa(j+1)-kappa(j))/ss
              dkdotgradpds(j) = (kdotgradp(j+1) - kdotgradp(j))/ss
              dbkdotgradbds(j) = (Bmag(j+1)*kdotgradb(j+1)-Bmag(j)*kdotgradb(j))/ss
       elseif(j==dimj)then
              ss=sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2)
              dbds(j) = (Bmag(j)-Bmag(j-1))/ss
              db2ds(j) = (Bmag(j)**2-Bmag(j-1)**2)/ss
              dkappads(j) = (kappa(j)-kappa(j-1))/ss
              dkdotgradpds(j) = (kdotgradp(j) - kdotgradp(j-1))/ss
              dbkdotgradbds(j) = (Bmag(j)*kdotgradb(j)-Bmag(j-1)*kdotgradb(j-1))/ss
       else
              ss = sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
              &  + sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
              dbds(j) = (Bmag(j+1)-Bmag(j-1))/ss
              db2ds(j) = (Bmag(j+1)**2-Bmag(j-1)**2)/ss
              dkappads(j) = (kappa(j+1)-kappa(j-1))/ss
              dkdotgradpds(j) = (kdotgradp(j+1) - kdotgradp(j-1))/ss
              dbkdotgradbds(j) = (Bmag(j+1)*kdotgradb(j+1)-Bmag(j-1)*kdotgradb(j-1))/ss
      endif
enddo
! d2b2ds2 & d2bds2

do j=1,dimj
       if(j==1)then
              dsdi(j) = -1.5*s0(j)+2.0*s0(j+1)-0.5*s0(j+2)
              d2sdi2(j) = 2.0*s0(j)-5.0*s0(j+1)+4.0*s0(j+2)-s0(j+3)
              d2b2ds2(j) = (2.0*Bmag(j)**2 -5*Bmag(j+1)**2 +4*Bmag(j+2)**2-Bmag(j+3)**2)/dss(j)**2
              d2bds2(j) = (2.0*Bmag(j) -5*Bmag(j+1) +4*Bmag(j+2)-Bmag(j+3))/dss(j)**2
       elseif(j==dimj)then
              dsdi(j) =  1.5*s0(j)-2.0*s0(j-1)+0.5*s0(j-2)
              d2sdi2(j) = 2.0*s0(j)-5.0*s0(j-1)+4.0*s0(j-2)-s0(j-3)
              d2b2ds2(j) = (2.0*Bmag(j)**2 -5*Bmag(j-1)**2 +4*Bmag(j-2)**2-Bmag(j-3)**2)/dss(j)**2
              d2bds2(j) = (2.0*Bmag(j) -5*Bmag(j-1) +4*Bmag(j-2)-Bmag(j-3))/dss(j)**2
       else
!              ss = 0.5*(((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
!                      + ((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2))
               ss = sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
               &  + sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
!              d2b2ds2(j) = (Bmag(j+1)**2 -2.0*Bmag(j)**2 + Bmag(j-1)**2)/ss
!              d2bds2(j) = (Bmag(j+1) -2.0*Bmag(j) + Bmag(j-1))/ss
!              d2b2ds2(j) = (Bmag(j+1)**2 -2.0*Bmag(j)**2 + Bmag(j-1)**2)/dss(j)**2
!              d2bds2(j) = (Bmag(j+1) -2.0*Bmag(j) + Bmag(j-1))/dss(j)
!               d2b2ds2(j) = (db2ds(j+1)-db2ds(j-1))/ss
!               d2bds2(j) = (dbds(j+1)-dbds(j-1))/ss
               dsdi(j) = (s0(j+1)-s0(j-1))/2.0
               d2sdi2(j) = s0(j+1)-2.0*s0(j)+s0(j-1)
               d2b2ds2(j) = (bmag(j+1)**2*(1.0-0.5*d2sdi2(j)/dsdi(j))&
                       -2.0*bmag(j)**2+bmag(j-1)**2*(1.0+0.5*d2sdi2(j)/dsdi(j)))/dsdi(j)**2
               d2bds2(j) = (bmag(j+1)*(1.0-0.5*d2sdi2(j)/dsdi(j))&
                       -2.0*bmag(j)+bmag(j-1)*(1.0+0.5*d2sdi2(j)/dsdi(j)))/dsdi(j)**2

       endif

enddo


! dbydy using divB
do j=1,dimj
       dbydy(j) = -(s_dxBx(j) + s_dzBz(j))
       C(j) = -1./Bmag(j)*(dbds(j) +dbydy(j))
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


! smooth
!call smooth(dbds,dimj,10,0)
!call smooth(d2bds2,dimj,10,0)
!call smooth(d2b2ds2,dimj,10,0)
!call smooth(kdotgradb,dimj,10,0)
!call smooth(kdotgradp,dimj,10,0)
!call smooth(dkdotgradpds,dimj,10,0)
!call gaussSmoother1D(s1,dkdotgradpds,dimj)
!call smooth(dbydy,dimj,10,0)
!call smooth(C,dimj,10,0)
! either of these changes the normal mode frequency
jmid = int((dimj+1)/2.)
jend = int(0.9*jmid)
!write(6,*)'smoothing dkappads from 2 to ',jend,' and ',2*jmid-jend,' to ',dimj-1
!call smooth(dkappads,dimj,10,2,jend)
!call smooth(dkappads,dimj,10,2*jmid-jend,dimj-1)
!write(6,*)'smoothing dbkdotgradbds from 2 to ',jend,' and ',2*jmid-jend,' to ',dimj-1
!call smooth(dbkdotgradbds,dimj,10,2,jend)
!call smooth(dbkdotgradbds,dimj,10,2*jmid-jend,dimj-1)
!call gaussSmoother1D(dbkdotgradbds,s1,dimj)
!write(6,*)'smoothing dkdotgradpds from 2 to ',jend,' and ',2*jmid-jend,' to ',dimj-1
!call smooth(dkdotgradpds,dimj,10,2,jend)
!call smooth(dkdotgradpds,dimj,10,2*jmid-jend,dimj-1)

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
!        c3(j) =fs(j)*(-0.5*d2b2ds2(j) + 2.0*fa(j)*dbds(j)**2)
        c3(j) =fs(j)*(-0.5*d2b2ds2(j) + fa(j)/(2.0*Bmag(j)**2)*db2ds(j)**2)
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
!              c5(j) = -Bmag(j)*(kdotgradb(j+1)-kdotgradb(j))/ss+(kappa(j+1)-kappa(j))/ss*Bmag(j)**2*(fa(j)-fs(j))+&
!                      kappa(j)*dbds(j)*Bmag(j)*(1.0-4.0*fs(j)**2)+dbydy(j)*(kdotgradb(j)-kappa(j)*Bmag(j))

        if(c5method==1)then
            c5(j) = -C(j)*Bmag(j)*(kdotgradb(j)-kappa(j)*Bmag(j))+dkappads(j)*Bmag(j)**2*(fa(j)-fs(j))+&
                      kappa(j)*dbds(j)*Bmag(j)*(2.0-4.0*fs(j)**2)-dbkdotgradbds(j)
        elseif(c5method==2)then
            c5(j) = C(j)*kdotgradp(j)+dkappads(j)*Bmag(j)**2*(fa(j)-fs(j))+&
                       kappa(j)*dbds(j)*Bmag(j)*(2.0-4.0*fs(j)**2)-dbkdotgradbds(j)
        elseif(c5method==3)then
            c5(j) = C(j)*kdotgradp(j)+dkappads(j)*Bmag(j)**2*(fa(j)-fs(j)-1.) &
                    -kappa(j)*dbds(j)*Bmag(j)*(4.0*fs(j)**2)+mu0*dkdotgradpds(j)
        elseif(c5method==4)then
            c5(j) = C(j)*kdotgradp(j)-2.0*dkappads(j)*Bmag(j)**2*fs(j) -4.0*fs(j)**2*kappa(j)*dbds(j)*Bmag(j)&
                    +mu0*dkdotgradpds(j)

        else
            write(5,*)'Wrong c5 method'
            stop
        end if

       elseif(j==dimj)then
               ss=sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2)
!              c5(j) = -kappa(j)*Bc4(j)*Bmag(j)  &
!              &     + (Bc3(j)*Bmag(j)-Bc3(j-1)*Bmag(j-1))/ss          &
!             &  + Bc4(j)*(khat(x,j)*(bhat(x,j)*s_dxBx(j)+bhat(z,j)*s_dxBz(j))  &
!             &  +         khat(z,j)*(bhat(x,j)*s_dzBx(j)+bhat(z,j)*s_dzBz(j)))
!              c5(j) = -Bmag(j)*(kdotgradb(j)-kdotgradb(j-1))/ss+(kappa(j)-kappa(j-1))/ss*Bmag(j)**2*(fa(j)-fs(j))+&
!                      kappa(j)*dbds(j)*Bmag(j)*(1.0-4.0*fs(j)**2)+dbydy(j)*(kdotgradb(j)-kappa(j)*Bmag(j))
         if(c5method==1)then
              c5(j) = -C(j)*Bmag(j)*(kdotgradb(j)-kappa(j)*Bmag(j))+dkappads(j)*Bmag(j)**2*(fa(j)-fs(j))+&
                       kappa(j)*dbds(j)*Bmag(j)*(2.0-4.0*fs(j)**2)-dbkdotgradbds(j)
          elseif(c5method==2)then
              c5(j) = C(j)*mu0*kdotgradp(j)+dkappads(j)*Bmag(j)**2*(fa(j)-fs(j))+&
                     kappa(j)*dbds(j)*Bmag(j)*(2.0-4.0*fs(j)**2)-dbkdotgradbds(j)
          elseif(c5method==3)then
              c5(j) = C(j)*mu0*kdotgradp(j)+dkappads(j)*Bmag(j)**2*(fa(j)-fs(j)-1.)&
                      -kappa(j)*dbds(j)*Bmag(j)*(4.0*fs(j)**2)+mu0*dkdotgradpds(j)
          elseif(c5method==4)then
            c5(j) = C(j)*mu0*kdotgradp(j)-2.0*dkappads(j)*Bmag(j)**2*fs(j) -4.0*fs(j)**2*kappa(j)*dbds(j)*Bmag(j)&
                    +mu0*dkdotgradpds(j)


        else
              write(5,*)'Wrong c5 method'
              stop
           end if
       else
              ss = sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
              &  + sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
!             c5(j) = -kappa(j)*Bc4(j)*Bmag(j)  &
!             &     + (Bc3(j+1)*Bmag(j+1)-Bc3(j-1)*Bmag(j-1))/ss          &
!             &  + Bc4(j)*(khat(x,j)*(bhat(x,j)*s_dxBx(j)+bhat(z,j)*s_dxBz(j))  &
!             &  +         khat(z,j)*(bhat(x,j)*s_dzBx(j)+bhat(z,j)*s_dzBz(j)))
!               c5(j) = -Bmag(j)*(kdotgradb(j+1)-kdotgradb(j-1))/ss+(kappa(j+1)-kappa(j-1))/ss*Bmag(j)**2*(fa(j)-fs(j))+&
!                       kappa(j)*Bmag(j)*dbds(j)*(1.0-4.0*fs(j)**2)+dbydy(j)*(kdotgradb(j)-kappa(j)*Bmag(j))
!              ! new version jan 16,18
         if(c5method==1)then
! uses magnetic field
              c5(j) = -C(j)*Bmag(j)*(kdotgradb(j)-kappa(j)*Bmag(j))+dkappads(j)*Bmag(j)**2*(fa(j)-fs(j))+&
                       kappa(j)*dbds(j)*Bmag(j)*(2.0-4.0*fs(j)**2)-dbkdotgradbds(j)
!              c5(j) = -C(j)*Bmag(j)*(kdotgradb(j)-kappa(j)*Bmag(j))+dkappads(j)*Bmag(j)**2*(fa(j)-fs(j))+&
!                       kappa(j)*dbds(j)*Bmag(j)*(2.0-4.0*fs(j)**2)-(Bmag(j+1)*kdotgradb(j+1)-Bmag(j-1)*kdotgradb(j-1))/ss


          elseif(c5method==2)then
!             c5(j) = C(j)*(khat(x,j)*s_dxp(j)+khat(z,j)*s_dzp(j))*mu0+(kappa(j+1)-kappa(j-1))/ss*Bmag(j)**2*(fa(j)-fs(j))+&
!                      kappa(j)*dbds(j)*Bmag(j)*(2.0-4.0*fs(j)**2)-(Bmag(j+1)*kdotgradb(j+1)-Bmag(j-1)*kdotgradb(j-1))/ss
             c5(j) = C(j)*mu0*kdotgradp(j)+dkappads(j)*Bmag(j)**2*(fa(j)-fs(j))+&
                      kappa(j)*dbds(j)*Bmag(j)*(2.0-4.0*fs(j)**2)-dbkdotgradbds(j)

! use pressure gradient in last term
          elseif(c5method==3)then
              c5(j) = C(j)*mu0*kdotgradp(j)+dkappads(j)*Bmag(j)**2*(fa(j)-fs(j)-1.)&
                      -kappa(j)*dbds(j)*Bmag(j)*(4.0*fs(j)**2)+mu0*dkdotgradpds(j)

           elseif(c5method==4)then
            c5(j) = C(j)*mu0*kdotgradp(j)-2.0*dkappads(j)*Bmag(j)**2*fs(j) -4.0*fs(j)**2*kappa(j)*dbds(j)*Bmag(j)&
                    +mu0*dkdotgradpds(j)

        else
             write(5,*)'Wrong c5 method'
              stop
           endif

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
!              c10(j) = 2.0*(Bmag(j)*kappa(j))**2*(fa(j)-fs(j)) -2.0*Bmag(j)*kappa(j)*kdotgradb(j) - &
!                       dbds(j)**2+Bmag(j)*d2b2ds2(j)-2.0*dbydy(j)*dbds(j)+Bmag(j)*(dbydy(j+1)-dbydy(j))/ss  - &
!                       dbydy(j)**2  
              c10(j) = 2.0*(Bmag(j)*kappa(j))**2*(fa(j)-fs(j)) -2.0*Bmag(j)*kappa(j)*kdotgradb(j) - &
                        (C(j)*Bmag(j))**2 -Bmag(j) * (Bmag(j+1)*C(j+1)-Bmag(j)*C(j))/ss
!            c10(j) = 2.0*(Bmag(j)*kappa(j))**2*(fa(j)-fs(j)-1.) +2.0*kappa(j)*kdotgradp(j) - &
!                      (C(j)*Bmag(j))**2 -Bmag(j) * (Bmag(j+1)*C(j+1)-Bmag(j)*C(1))/ss
       elseif(j==dimj)then
              ss=sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2)
!              c10(j) = 2.0*Bmag(j)*Bc3(j)*kappa(j)-Bc4(j)**2 &
!              &      + Bmag(j)*(Bc4(j)-Bc4(j-1))/ss
!              c10(j) = 2.0*(Bmag(j)*kappa(j))**2*(fa(j)-fs(j)) -2.0*Bmag(j)*kappa(j)*kdotgradb(j) - &
!                       dbds(j)**2+Bmag(j)*d2b2ds2(j)-2.0*dbydy(j)*dbds(j)+Bmag(j)*(dbydy(j)-dbydy(j-1))/ss  - &
!                       dbydy(j)**2  
              c10(j) = 2.0*(Bmag(j)*kappa(j))**2*(fa(j)-fs(j)) -2.0*Bmag(j)*kappa(j)*kdotgradb(j) - &
                       (C(j)*Bmag(j))**2 -Bmag(j) * (Bmag(j)*C(j)-Bmag(j-1)*C(j-1))/ss
!            c10(j) = 2.0*(Bmag(j)*kappa(j))**2*(fa(j)-fs(j)-1.) +2.0*kappa(j)*kdotgradp(j) - &
!                      (C(j)*Bmag(j))**2 -Bmag(j) * (Bmag(j)*C(j)-Bmag(j-1)*C(j-1))/ss
       else
              ss = sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
              &  + sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
!              c10(j) = 2.0*Bmag(j)*Bc3(j)*kappa(j)-Bc4(j)**2 &
!              &      + Bmag(j)*(Bc4(j+1)-Bc4(j-1))/ss

!              c10(j) = 2.0*(Bmag(j)*kappa(j))**2*(fa(j)-fs(j)) -2.0*Bmag(j)*kappa(j)*kdotgradb(j) - &
!                       dbds(j)**2+Bmag(j)*d2bds2(j)-2.0*dbydy(j)*dbds(j)+Bmag(j)*(dbydy(j+1)-dbydy(j-1))/ss  - &
!                       dbydy(j)**2                      
              c10(j) = 2.0*(Bmag(j)*kappa(j))**2*(fa(j)-fs(j)) -2.0*Bmag(j)*kappa(j)*kdotgradb(j) - &
                        (C(j)*Bmag(j))**2 -Bmag(j) * (Bmag(j+1)*C(j+1)-Bmag(j-1)*C(j-1))/ss
!            c10(j) = 2.0*(Bmag(j)*kappa(j))**2*(fa(j)-fs(j)-1.) +2.0*kappa(j)*kdotgradp(j) - &
!                      (C(j)*Bmag(j))**2 -Bmag(j) * (Bmag(j+1)*C(j+1)-Bmag(j-1)*C(j-1))/ss
       endif
enddo

! smooth
! frt now smooth c3
!       call smooth(c1,dimj,10,2,dimj-1)
!       call smooth(c2,dimj,10,2,dimj-1)
!       call smooth(c3,dimj,10,2,dimj-1)
!       call smooth(c4,dimj,10,2,dimj-1)
!       call smooth(c5,dimj,10,2,dimj-1)
!      call smooth(c6,dimj,10,2,dimj-1)
!       call smooth(c7,dimj,10,2,dimj-1)
!       call smooth(c8,dimj,10,2,dimj-1)
!       call smooth(c9,dimj,10,2,dimj-1)
!       call smooth(c10,dimj,10,2,dimj-1)
! does not make much difference
!        call symmetry(c1,dimj)
!        call symmetry(c2,dimj)
!        call symmetry(c3,dimj)
!        call symmetry(c4,dimj)
!        call symmetry(c5,dimj)
!        call symmetry(c6,dimj)
!        call symmetry(c7,dimj)
!        call symmetry(c8,dimj)
!        call symmetry(c9,dimj)
!        call symmetry(c10,dimj)



if(method ==3) then
        write(6,*)'****** computing coefficients using the 3rd approach..****'
        jmid = (dimj+1)/2
        do j=1,dimj
         xequator = state(x,jmid)
         eps0 =5.0e-2
         g3s(j) = (pBack(xequator+eps,0.0,-1)-pBack(xequator-eps,0.0,-1))/(2.0*eps)
         g3s(j) = g3s(j) *mu0/Bz(xequator,0.0,-1)
         ! renormalize 
         eps = eps0*(sqrt(state(x,j)**2+state(z,j)**2)/abs(xequator))**3
         ! now compute fs
         if(j==1)then
                dx = state(x,2)-state(x,1)
                dz = state(z,2)-state(z,1)
         elseif(j==dimj)then
                dx = state(x,dimj)-state(x,dimj-1)
                dz = state(z,dimj)-state(z,dimj-1)
         else
                dx = state(x,j+1)-state(x,j-1)
                dz = state(z,j+1)-state(z,j-1)
         end if
         ds = sqrt(dx**2+dz**2)
         dx = dx/ds
         dz = dz/ds

         x0 = state(x,j) - eps*sign(dz,state(z,j))
         y0 = 0.0
         z0 = state(z,j) + eps*sign(dx,state(z,j))
!         x01 = state(x,j) + eps*sign(dz,state(z,j))
!         y01 = 0.0
!         z01 = state(z,j) - eps*sign(dx,state(z,j))

         ! reset if inside the boundary, including eps
         if(wallBound)then
                 if(x0 < -boundary)then
                  x0 = -boundary
                  eps = sqrt((x0-state(x,j))**2+(z0-state(z,j))**2)
                 endif
         else
                 if(abs(sqrt(x0**2+z0**2)-boundary)<1.0e-8)then ! inside boundary
                         write(6,*)' inside boundary, resetting location..'
                         x0 = x0 *boundary/sqrt(x0**2+z0**2)
                         z0 = z0 *boundary/sqrt(x0**2+z0**2)
                         eps = sqrt((x0-state(x,j))**2+(z0-state(z,j))**2)
                 end if
         endif


         if(z0 > 0.0)then
                idir = -1
         else
                idir = 1
         endif
         call find_equator_tracerRK4a(x0,y0,z0,idir,isopen,xe)
!         call find_equator_tracerRK4a(x01,y01,z01,idir,isopen,xe1)
         write(6,'(a,i4,a,f8.4,a,f8.4,a,f8.4)')' j=',j,' x=',state(x,j),' z=',state(z,j),' xequator=',xe
         if(xe >= 100)then
                 write(6,*)'equatorial tracing not found..'
                 write(6,*)' j=',j,' for point on line:',state(x,j),state(z,j)
                 write(6,*)' isopen =',isopen
                 write(6,*)' starting point =',x0,z0
         end if

         f3s(j) =s_Bz(jmid)/sqrt(s_Bx(j)**2+s_Bz(j)**2)*abs(xe-xequator)/eps

        end do

        ! now compute the new c5 and c10
        do j=1,dimj
         if(j==1)then
              ss=sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
              c5(j) = -2.0*fs(j)*dkappads(j) -4.0*kappa(j)*bmag(j)*fs(j)**2*dbds(j) &
                      +bmag(j) *g3s(j)*(f3s(j+1)-f3s(j))/ss-f3s(j)*g3s(j)*dbydy(j)
        elseif(j==dimj)then
              ss=sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2)
              c5(j) = -2.0*fs(j)*dkappads(j) -4.0*kappa(j)*bmag(j)*fs(j)**2*dbds(j) &
                      +bmag(j) *g3s(j)*(f3s(j)-f3s(j-1))/ss-f3s(j)*g3s(j)*dbydy(j)
        else
              ss = sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
              &  + sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
              c5(j) = -2.0*fs(j)*dkappads(j)-4.0*kappa(j)*bmag(j)*fs(j)**2*dbds(j) &
                      +bmag(j) *g3s(j)*(f3s(j+1)-f3s(j-1))/ss-f3s(j)*g3s(j)*dbydy(j)
        endif

        c10(j) = -4.0*fs(j) +bmag(j)**2*kappa(j)**2 + 2.0 *bmag(j)*g3s(j) *kappa(j) * f3s(j) 
        if(j==1)then
              ss=sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
             c10(j) = c10(j)   &
                      -dbds(j)**2+Bmag(j)*d2b2ds2(j)-2.0*dbydy(j)*dbds(j)+Bmag(j)*(dbydy(j+1)-dbydy(j))/ss  - &
                      dbydy(j)**2  
        elseif(j==dimj)then
              ss=sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2)
             c10(j) = c10(j)   & 
                      - dbds(j)**2+Bmag(j)*d2b2ds2(j)-2.0*dbydy(j)*dbds(j)+Bmag(j)*(dbydy(j)-dbydy(j-1))/ss  - &
                      dbydy(j)**2  
        else
              ss = sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
              &  + sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
              c10(j) = c10(j)   &
                     -dbds(j)**2+Bmag(j)*d2bds2(j)-2.0*dbydy(j)*dbds(j)+Bmag(j)*(dbydy(j+1)-dbydy(j-1))/ss  - &
                      dbydy(j)**2                     

         end if
        end do 

end if


! ------setup fop matix--------
fop = 0.0

!upper left quadrant
do j=2,dimj-1
!       ss = 0.5*sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
!       &  + 0.5*sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
!       dsdi = (s0(j+1)-s0(j-1))/2.0
!       d2sdi2 = s0(j+1)-2.0*s0(j)+s0(j-1)
       ss = dsdi(j)

       delta = d2sdi2(j)/dsdi(j)
       if(iopen ==0)then
          fop(j,j) = c3(j)-2.0*c1(j)/ss**2
          if(j.ne.dimj-1)fop(j,j+1)= c1(j)*(1.0-delta/2.)/ss**2+c2(j)/(2.0*ss)
          if(j.ne.2)fop(j,j-1) = c1(j)*(1.0+delta/2.)/ss**2-c2(j)/(2.0*ss)
       else if(iopen==1)then
          if(j.ne.dimj-1)fop(j,j+1)= c1(j)*(1.0-delta/2.)/ss**2+c2(j)/(2.0*ss)
          fop(j,j) = c3(j)-2.0*c1(j)/dsdi(j)**2
          if(j.ne.2)fop(j,j-1) = c1(j)*(1.0+delta/2.)/ss**2-c2(j)/(2.0*ss)

          if(j==2)then
             fop(j,j) = -2.0*c1(j)/ss**2+c3(j)
             fop(j,j+1)= c1(j)*(1.0-delta/2.)/ss**2+c2(j)/2.0/ss
          end if
                  
          if(j==dimj-1)then
             fop(j,j) = -2.0*c1(j)/ss**2+c3(j)
             fop(j,j-1) = c1(j)*(1.0+delta/2.)/ss**2-c2(j)/2.0/ss
          end if
       else
           stop 'wrong iopen 1' 
       end if

!       if((j.ne.dimj-1.and.iopen==0).or.(j<dimj-2.and.iopen==1))fop(j,j+1)= c1(j)/(ss**2)+c2(j)/(2.0*ss)
!       fop(j,j) = c3(j)-2.0*c1(j)/(ss**2)
!       if((j.ne.2.and.iopen==0).or.(j > 3.and.iopen==1))fop(j,j-1) = c1(j)/(ss**2)-c2(j)/(2.0*ss)

!       if((j.ne.dimj-1.and.iopen==0).or.(j<dimj-2.and.iopen==1))fop(j,j+1)= c1(j)*(1.0-0.5*d2sdi2/dsdi)/dsdi**2+c2(j)/(2.0*ss)
!       fop(j,j) = c3(j)-2.0*c1(j)/dsdi**2
!       if((j.ne.2.and.iopen==0).or.(j > 3.and.iopen==1))fop(j,j-1) = c1(j)*(1.0+0.5*d2sdi2/dsdi)/dsdi**2-c2(j)/(2.0*ss)

!       if(iopen==1.and.j==2)then
!               fop(j,j) = -2.0*c1(j)/dsdi**2
!!               fop(j,j-1) = 0.0
!!               fop(j,j+1) = c2(j)/ss
!       end if
!       if(iopen==1.and.j==dimj-1)then
!                fop(j,j) = -2.0*c1(j)/dsdi**2
!!                fop(j,j-1) = -c2(j)/ss
!!                fop(j,j+1) = 0.0
!        end if

enddo


!upper right quadrant
do j=2,dimj-1
!       ss = 0.5*sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
!       &  + 0.5*sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
       ss = dsdi(j)
       delta = d2sdi2(j)/dsdi(j)
       if(iopen==0)then
         if(j.ne.dimj-1)fop(j,j+1+dimj) = c4(j)/(2.*ss)
         if(j.ne.2)fop(j,j-1+dimj) = -c4(j)/(2.*ss)
         fop(j,j+dimj) = c5(j)
       else if(iopen==1)then
         if(j.ne.dimj-1)fop(j,j+1+dimj) = c4(j)/(2.*ss)
         if(j.ne.2)fop(j,j-1+dimj) = -c4(j)/(2.*ss)
         fop(j,j+dimj) = c5(j)

         if(j==2)then
           geom = state(x,1)/(2.0*state(z,1))
           alpha = 1.5 +  dsdi(1)*C(1)
           fop(j,j+dimj) = 2.0*c1(j)*(1.0+delta/2.0)*geom/alpha/ss**2 &
                   -c2(j)*geom/alpha/ss-c4(j)/alpha/ss + c5(j)
           fop(j,j+1+dimj) = -c1(j)*(1.0+delta/2.0)*geom/2.0/ss**2/alpha&
                   +c2(j)*geom/4.0/alpha/ss+c4(j)*(1.0+1.0/(2.0*alpha))/2./ss 
         end if

         if(j==dimj-1)then
           geom = state(x,dimj)/(2.0*state(z,dimj))
           alpha = 1.5 - dsdi(dimj)*C(dimj)
           fop(j,j+dimj) = 2.0*c1(j)*(1.0-delta/2.0)*geom/alpha/ss**2 &
                   +c2(j)*geom/alpha/ss+c4(j)/alpha/ss + c5(j)
           fop(j,j-1+dimj) = -c1(j)*(1.0-delta/2.0)*geom/2.0/ss**2/alpha&
                   -c2(j)*geom/4.0/alpha/ss-c4(j)*(1.0+1.0/(2.0*alpha))/2.0/ss 
         end if
       else
         stop 'wrong iopen 2'
       end if

!       if(j.ne.dimj-1)fop(j,j+1+dimj) = c4(j)/(2.0*ss)
!       fop(j,j+dimj) = c5(j)
!       if(j.ne.2)fop(j,j-1+dimj) = -c4(j)/(2.0*ss)
  
!       if(iopen==1.and.j==2)then
!        fop(j,j+dimj) = c5(j) + C(j)*c4(j)
!        fop(j,j+1+dimj) = 0.0
!       end if
!       if(iopen==1.and.j==dimj-1)then
!       fop(j,j+dimj) = c5(j) + C(j)*c4(j)
!       fop(j,j-1+dimj) = 0.0
!       end if
enddo

!lower left quadrant
do j=2,dimj-1
!       ss = 0.5*sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
!       &  + 0.5*sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
       ss = dsdi(j)
       if(iopen==0)then
         if(j.ne.dimj-1)fop(j+dimj,j+1) =  c6(j)/(2.0*ss)
         fop(j+dimj,j) = c7(j)
         if(j.ne.2)fop(j+dimj,j-1) = -c6(j)/(2.0*ss)
       else if(iopen==1)then
         if(j.ne.dimj-1)fop(j+dimj,j+1) =  c6(j)/(2.0*ss)
         fop(j+dimj,j) = c7(j)
         if(j.ne.2)fop(j+dimj,j-1) = -c6(j)/(2.0*ss)

         if(j==2)then
           fop(j+dimj,j) = c7(j) 
           fop(j+dimj,j+1) = c6(j)/(2.0*ss) 
         end if
         if(j==dimj-1)then
           fop(j+dimj,j) = c7(j)
           fop(j+dimj,j-1) = -c6(j)/(2.0*ss)
         end if
       else
          stop 'wrong iopen 3'
       end if
!       if((j.ne.dimj-1.and.iopen==0).or.(j<dimj-2.and.iopen==1))fop(j+dimj,j+1) =  c6(j)/(2.0*ss)
!       fop(j+dimj,j) = c7(j)
!       if((j.ne.2.and.iopen==0).or.(j > 3.and.iopen==1))fop(j+dimj,j-1) = -c6(j)/(2.0*ss)
!       if(iopen==1.and.j==2)then
!!             fop(j+dimj,j-1) = 0.0
!             fop(j+dimj,j) = 0.0
!!             fop(j+dimj,j+1) =  c6(j)/(ss)
!       end if
!       if(iopen==1.and.j==dimj-1)then
!!             fop(j+dimj,j-1) = - c6(j)/(ss) 
!             fop(j+dimj,j) = 0.0
!!            fop(j+dimj,j+1) = 0.0
!       end if

enddo

!lower right quadrant
do j=2,dimj-1
!       ss = 0.5*sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
!       &  + 0.5*sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
       ss = dsdi(j)
       delta = d2sdi2(j)/dsdi(j)

       if(iopen==0)then
         if(j.ne.dimj-1) fop(j+dimj,j+1+dimj) = c8(j)*(1.0-delta/2.0)/ss**2+c9(j)/(2.0*ss)
         fop(j+dimj,j+dimj) = c10(j)-2.0*c8(j)/ss**2
         if(j.ne.2)fop(j+dimj,j-1+dimj) = c8(j)*(1.0+delta/2.0)/ss**2-c9(j)/(2.0*ss)
       elseif(iopen==1)then
         if(j.ne.dimj-1) fop(j+dimj,j+1+dimj) = c8(j)*(1.0-delta/2.0)/ss**2+c9(j)/(2.0*ss)
         fop(j+dimj,j+dimj) = c10(j)-2.0*c8(j)/ss**2
         if(j.ne.2)fop(j+dimj,j-1+dimj) = c8(j)*(1.0+delta/2.0)/ss**2-c9(j)/(2.0*ss)

         if(j==2)then
          geom = state(x,1)/(2.0*state(z,1))
          alpha = 1.5 +  dsdi(1)*C(1)
          write(6,*)'1-alpha =',alpha
          fop(j+dimj,j+dimj) = -c6(j)/(alpha*ss)+c8(j)*(-2.0+2.0/alpha+delta/alpha)/ss**2 &
                  -c9(j)/(alpha*ss) +c10(j)
          fop(j+dimj,j+1+dimj) = c6(j)*geom/(4.0*alpha*ss)&
                  +c8(j)*(1.0-1.0/(2.0*alpha)-delta/2.0-delta/(4.0*alpha))/ss**2&
                  +c9(j)*(1.0+1.0/(2.0*alpha))/(2.0*ss)
         end if

         if(j==dimj-1)then
          geom = state(x,dimj)/(2.0*state(z,dimj))
          alpha = 1.5 - dsdi(dimj)*C(dimj)
          write(6,*)'dimj-alpha =',alpha
          fop(j+dimj,j+dimj) = c6(j)/(alpha*ss)+c8(j)*(-2.0+2.0/alpha+delta/alpha)/ss**2 &
                  +c9(j)/(alpha*ss) +c10(j)
          fop(j+dimj,j-1+dimj) = -c6(j)*geom/(4.0*alpha*ss)&
                  +c8(j)*(1.0-1.0/(2.0*alpha)+delta/2.0+delta/(4.0*alpha))/ss**2&
                  -c9(j)*(1.0+1.0/(2.0*alpha))/(2.0*ss)
         end if
       else
               stop 'wrong iopen 4'
       end if


!      if(j.ne.dimj-1) fop(j+dimj,j+1+dimj) = c8(j)/(ss**2)+c9(j)/(2.0*ss)
!       fop(j+dimj,j+dimj) = c10(j)-2.0*c8(j)/(ss**2)
!       if(j.ne.2)fop(j+dimj,j-1+dimj) = c8(j)/(ss**2)-c9(j)/(2.0*ss)

!       if(j.ne.dimj-1) fop(j+dimj,j+1+dimj) = c8(j)*(1.0-0.5*d2sdi2/dsdi)/dsdi**2+c9(j)/(2.0*ss)
!       fop(j+dimj,j+dimj) = c10(j)-2.0*c8(j)/dsdi**2
!       if(j.ne.2)fop(j+dimj,j-1+dimj) = c8(j)*(1.0+0.5*d2sdi2/dsdi)/dsdi**2-c9(j)/(2.0*ss)

!       if(iopen==1.and.j==2)then
!!        fop(j+dimj,j+dimj) = c10(j)-1.0*c8(j)/(ss**2)*(1.0+C(j)*ss) +c9(j)*C(j)
!!        fop(j+dimj,j+1+dimj) = 1.0*c8(j)/(ss**2)
!        fop(j+dimj,j+dimj) = c10(j)-1.0*c8(j)/(dsdi**2)*(1.0+C(j)*ss) +c9(j)*C(j)
!        fop(j+dimj,j+1+dimj) = 1.0*c8(j)/(dsdi**2)*(1.0-0.5*d2sdi2/dsdi)
!       end if
!       if(iopen==1.and.j==dimj-1)then
!!        fop(j+dimj,j+dimj) = c10(j)-1.0*c8(j)/(ss**2)*(1.0-C(j)*ss) +c9(j)*C(j)
!!        fop(j+dimj,j-1+dimj) = 1.0*c8(j)/(ss**2)
!        fop(j+dimj,j+dimj) = c10(j)-1.0*c8(j)/(dsdi**2)*(1.0-C(j)*ss) +c9(j)*C(j)
!        fop(j+dimj,j-1+dimj) = 1.0*c8(j)/(dsdi**2)*(1.0+0.5*d2sdi2/dsdi)

!       end if

enddo

open(100,file=trim(datapath)//'.debug.dat',status='unknown',recl=500)

open(101,file=trim(datapath)//'.debug1.dat',status='unknown',recl=500)

open(102,file=trim(datapath)//'.debug2.dat',status='unknown',recl=500)

open(103,file=trim(datapath)//'.debug3.dat',status='unknown',recl=500)

open(104,file=trim(datapath)//'.mass.dat',status='unknown',recl=500)

write(6,*)'dimj =',dimj
write(100,*)'j, state(x,j), state(z,j), state(p,j), state(d,j), s_Bx(j), s_Bz(j), s_ptotal(j), Bmag(j) s_dxP sdzp'
do j=2,dimj-1
write(100,*)j,state(x,j),state(z,j), state(p,j), state(d,j), s_Bx(j), s_Bz(j), s_ptotal(j), Bmag(j), s_dxp(j),s_dzp(j)
enddo
write(101,*)'j,khat(x,j), khat(z,j), cs(j), ca(j), Bc1(j), Bc2(j), Bc3(j), Bc4(j) '
do j=2,dimj-1
write(101,*) j,khat(x,j),khat(z,j),cs(j),ca(j),Bc1(j),Bc2(j),Bc3(j),Bc4(j)
enddo
write(102,*)'j, c1(j), c2(j), c3(j), c4(j), c5(j), c6(j), c7(j), c8(j), c9(j), c10(j) C'
do j=2,dimj-1
write(102,*)j,c1(j), c2(j), c3(j), c4(j), c5(j), c6(j), c7(j), c8(j), c9(j), c10(j), C(j)
enddo


write(103,*)'j, s_dxBx(j), s_dxBz(j), s_dzBx(j), s_dzBz(j), s_dxPtot(j), s_dzPtot(j) dbydy(j) dbds(j)', &
'd2bds2(j) d2b2ds2(j) kappa(j) kdotgradb kdotgradp'
do j=2,dimj-1
write(103,*)j,s_dxBx(j), s_dxBz(j), s_dzBx(j), s_dzBz(j), s_dxPtot(j),&
s_dzPtot(j),dbydy(j),dbds(j),d2bds2(j),d2b2ds2(j),kappa(j),kdotgradb(j),kdotgradp(j)
enddo

write(104,*)' j x z mass rho b'
do j=2,dimj-1
 write(104,*)j,state(x,j), state(z,j),state(mf,j),state(d,j),Bmag(j)
enddo

close(100)
close(101)
close(102)
close(103)
close(104)

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
deallocate(kvec)
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
deallocate(dsdi)
deallocate(d2sdi2)

deallocate(C)
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
deallocate(s_dxP)
deallocate(s_dzP)
deallocate(dbydy)
deallocate(dbds)
deallocate(dkappads)
deallocate(d2b2ds2)
deallocate(d2bds2)
deallocate(kdotgradb)
deallocate(kdotgradp)
deallocate(dkdotgradpds)
deallocate(dbkdotgradbds)
deallocate(dss)
deallocate(f3s)
deallocate(g3s)
deallocate(s0)

end subroutine fileigen

subroutine smooth(a,idim,ipass,starti,iend)
        implicit none
        integer :: idim,ipass,istart,iend,starti
        real :: a(idim),atemp(idim)
        integer :: i,it
        real, parameter :: weight =2.0

        atemp = a

        if(starti==0)then
                istart =2
                iend = idim-1
        else
                istart = starti
!                iend = idim - istart + 1
        end if

        do it =1,ipass

        do i=istart,iend
          atemp(i) = (a(i-1) + weight*a(i) + a(i+1))/(2.+weight)
        end do

       a = atemp
       end do
       
end subroutine smooth


subroutine gaussSmoother1D(a,x,gi)
! For a symmetrical kernel of two arguments $K(x_1,x_2) = K(x_2,x_1)$ the smoothing routine calculates
! a smoothed interpolant $\hat y$ of a data set $y$ at a point $x_o$
! $\hat y(x_o) = \frac{\sum\limits_{i=1}^N K(x_o,x_i)y(x_i) }{\sum\limits_{i=1}^N K(x_o,x_i)}$
! This procedure may be commutatively applied to each dimension of an N-dimensional data set
!
! a=original grid, x=xaxis grid, z= z axis grid, gi=x grid size, gj = z grid size, c=smoothed output
! 1D version frt 2/18
!
     implicit none
     integer,intent(in) :: gi
     real,dimension(gi),intent(inout) :: a
     real,dimension(gi),intent(in) :: x
     real,dimension(gi) :: b
     integer :: i,l
     real :: c0,c1
     real :: kern1d !kernel function
     real :: w !smoothing weight
     real :: kf ! kernal function result
     real, parameter :: sigma = 3.0 ! how many sigmas to do smoothing
     
     !w = 0.25 !0.25 smoothed out grid scale oscillations and preserved force balance
     w = (x(2)-x(1)) !scales with different grid resolutions, set to max grid spacing
     
     !smooth over x-direction
     do i = 1,gi
         c0 = 0.0; c1 = 0.0;
         !perform numerator sum
         do l = 1,gi
          if(abs(x(i)-x(l))<sigma*w)then ! for speed, limit the range to sigma
          kf = kern1d(x(i),x(l),w)
          c0 = c0 + a(l)*kf
       !perform denominator sum
          c1 = c1 + kf
         end if
        enddo
        b(i) = c0/c1
     enddo
! finally return the new smoothed value
     a = b
     return
     
end subroutine gaussSmoother1D

real function kern1d(x0,x1,w)     
     real,intent(in) :: x0,x1,w
     real,save :: norm = 1/(sqrt(8.0*atan(1.0))) !normalization const.
     !gaussian kernel for smoothing routine
     kern1d = (norm/w)*exp(-(x0-x1)**2/(2*w**2))
end function kern1d

subroutine symmetry(a,ndim)
        ! enforces symmetry/antisymmetry
        implicit none
        integer :: ndim
        real :: a(ndim), atemp(ndim),average
        integer :: n,nmid

        nmid = (ndim+1)/2

        do n=1,nmid-1
          average = (abs(a(n))+abs(a(ndim-n+1)))/2.0
          a(n) = sign(average , a(n))
          a(ndim-n+1) = sign(average ,a(ndim-n+1))
         end do
         ! middle point

!         if (a(nmid-1)*a(nmid+1) < 0.0)a(nmid) = 0.0
         a(nmid) = 0.5*(a(nmid+1)+a(nmid-1))

         return
         end subroutine symmetry
 ! 1 D Interpolation program
! 2/18 frt
real function interp1d(x,xi,yi,idim)
        implicit none
        integer, intent(in):: idim
        integer :: i,ig
        real, intent(in):: xi(idim),yi(idim),x
        real :: interpg

        interp1d = 0.0

        ! first check to see if data is montonic
        if(xi(2) < xi(1))then
                write(6,*)' Interp1d: x seems not montonic'
                return
        end if

        ! check to see if the data is in range
        if(x < xi(1) .or. x > xi(idim))then
                write(6,*)' Interp1d: a x out of range'
                return
        endif

        ! now find index
        do i=1,idim-1
         if(x >= xi(i) .and. x <= xi(i+1))then
                 ig = i
                 exit
         end if
        end do

         ! now interpolate

!         interp1d = (-yi(i) * (x-xi(i+1)) + yi(i+1)*(x-xi(i)))/(xi(i+1)-xi(i))
        if(ig>2.and.ig<idim-2)then
              interp1d = interpg(x,xi(ig-1:ig+2),yi(ig-1:ig+2),4)
!              interp1d = interpg(x,xi(ig:ig+1),yi(ig:ig+1),2)
        else
              interp1d = interpg(x,xi(ig:ig+1),yi(ig:ig+1),2)
        end if


         return
end function interp1d

real function interpg(x,xi,yi,idim)
      implicit none
      integer,intent(in) :: idim
      real,intent(in) :: xi(idim),yi(idim),x
      real :: v
      integer :: i,j
!% General function to interpolate between data points
!% using Lagrange polynomial
!% Inputs
!%   x    Vector of x coordinates of data points (imax values)
!%   y    Vector of y coordinates of data points (imax values)
!%   xi   The x value where interpolation is computed
!% Output
!%   yi   The interpolation polynomial evaluated at xi

!%* Calculate yi = p(xi) using Lagrange polynomial

   interpg =0.0;
   do i=1,idim 
     v=1;
     do j=1,idim   
        if(i.ne.j)then
          v=v*(x-xi(j))/(xi(i)-xi(j));
        end if
     enddo

        interpg=interpg+v*yi(i);
   end do
return 
end function interpg
