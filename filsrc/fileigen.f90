subroutine fileigen
!eigenvalue and eigenvector solver for the linearized EOM
use stateMod
use stateIndex
use constants
use dnams
implicit none
real :: Bx,Bz,pBack,dxBx,dxBz,dzBx,dzBz,ptotal,dxPtot,dzPtot
real,allocatable,dimension(:) :: s_Bx,s_Bz,s_pBack,s_dxBx,s_dxBz,s_dzBx
real,allocatable,dimension(:) :: s_dzBz,s_ptotal,s_dxPtot,s_dzPtot
real,allocatable,dimension(:,:) :: khat,bhat
real,allocatable,dimension(:) :: kappa,ck,Bmag,Bc1,Bc2,Bc3,Bc4,cs,ca
real,allocatable,dimension(:) :: C1,C2,C3,C4,C5,C6,C7,C8,C9,C10
real,allocatable,dimension(:,:) :: fop
integer :: i,j
real :: ss,xx,zz

write(*,*) 'Calculating force operator'

allocate(khat(2,dimj))
allocate(bhat(2,dimj))

allocate(kappa(dimj))
allocate(Bmag(dimj))
allocate(ck(dimj))
allocate(cs(dimj))
allocate(ca(dimj))
allocate(Bc1(dimj))
allocate(Bc2(dimj))
allocate(Bc3(dimj))
allocate(Bc4(dimj))

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


do j=1,dimj
	Bmag(j) = sqrt(s_Bx(j)**2+s_Bz(j)**2)
enddo

do j=1,dimj
	bhat(x,j) = s_Bx(j)/Bmag(j)
	bhat(z,j) = s_Bz(j)/Bmag(j)
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

! do j=1,dimj
! 	ck(j)=khat(x,j)*(khat(x,j)*(s_dxBx(j)/Bmag(j)-s_Bx(j)/(Bmag(j)**2)*(bhat(x,j)*s_dxBx(j)+bhat(z,j)*s_dxBz(j))) &
! 				  &+ khat(z,j)*(s_dzBx(j)/Bmag(j)-s_Bx(j)/(Bmag(j)**2)*(bhat(x,j)*s_dzBx(j)+bhat(z,j)*s_dzBz(j))))&
! 	   & +khat(z,j)*(khat(x,j)*(s_dxBz(j)/Bmag(j)-s_Bz(j)/(Bmag(j)**2)*(bhat(x,j)*s_dxBx(j)+bhat(z,j)*s_dxBz(j))) &
! 			      &+ khat(z,j)*(s_dzBz(j)/Bmag(j)-s_Bz(j)/(Bmag(j)**2)*(bhat(x,j)*s_dzBx(j)+bhat(z,j)*s_dzBz(j))))
! enddo

do j=1,dimj
	ck(j)= (1/Bmag(j))*(khat(x,j)*(khat(x,j)*s_dxBx(j)+khat(z,j)*s_dzBx(j)) &
		&			   +khat(z,j)*(khat(x,j)*s_dxBz(j)+khat(z,j)*s_dzBz(j)))
enddo

do j=1,dimj
	cs(j) = (cm/Re)*sqrt(gamma*state(p,j)/(nano*state(d,j)*mi*cmC))
enddo

do j=1,dimj
	ca(j) = (cm/Re)*Bmag(j)/(nano*sqrt(state(d,j)*mi*mu0*cmC/nano))
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
	c1(j)=Bmag(j)*Bc1(j)
enddo

do j=1,dimj
	if(j==1)then
		ss=sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
		c2(j) = Bmag(j)*Bc2(j)+(Bc1(j+1)*Bmag(j+1)-Bc1(j)*Bmag(j))/ss
	elseif(j==dimj)then
		ss=sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2)
		c2(j) = Bmag(j)*Bc2(j)+(Bc1(j)*Bmag(j)-Bc1(j-1)*Bmag(j-1))/ss
	else
		ss = sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
		&  + sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
		c2(j) = Bmag(j)*Bc2(j)+(Bc1(j+1)*Bmag(j+1)-Bc1(j-1)*Bmag(j-1))/ss
	endif
enddo


do j=1,dimj
	if(j==1)then
		ss=sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
		c3(j) = (Bc2(j+1)*Bmag(j+1)-Bc1(j)*Bmag(j))/ss
	elseif(j==dimj)then
		ss=sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2)
		c3(j) = (Bc2(j)*Bmag(j)-Bc1(j-1)*Bmag(j-1))/ss
	else
		ss = sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
		&  + sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
		c3(j) = (Bc2(j+1)*Bmag(j+1)-Bc1(j-1)*Bmag(j-1))/ss
	endif
enddo

do j=1,dimj
	c4(j)= -(Bmag(j)**2)*kappa(j)+Bc3(j)*Bmag(j)+Bmag(j) &
		&  *(khat(x,j)*(bhat(x,j)*s_dxBx(j)+bhat(z,j)*s_dxBz(j)) &
		&  + khat(z,j)*(bhat(x,j)*s_dzBx(j)+bhat(z,j)*s_dzBz(j)))
enddo

do j=1,dimj
	if(j==1)then
		ss=sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
		c5(j) = -kappa(j)*Bc4(j)*Bmag(j)  &
		&     + (Bc3(j+1)*Bmag(j+1)-Bc3(j)*Bmag(j))/ss          &
		&  + Bc4(j)*(khat(x,j)*(bhat(x,j)*s_dxBx(j)+bhat(z,j)*s_dxBz(j))  &
		&  +         khat(z,j)*(bhat(x,j)*s_dzBx(j)+bhat(z,j)*s_dzBz(j)))
	elseif(j==dimj)then
		ss=sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2)
		c5(j) = -kappa(j)*Bc4(j)*Bmag(j)  &
		&     + (Bc3(j)*Bmag(j)-Bc3(j-1)*Bmag(j-1))/ss          &
		&  + Bc4(j)*(khat(x,j)*(bhat(x,j)*s_dxBx(j)+bhat(z,j)*s_dxBz(j))  &
		&  +         khat(z,j)*(bhat(x,j)*s_dzBx(j)+bhat(z,j)*s_dzBz(j)))
	else
		ss = sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
		&  + sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
		c5(j) = -kappa(j)*Bc4(j)*Bmag(j)  &
		&     + (Bc3(j+1)*Bmag(j+1)-Bc3(j-1)*Bmag(j-1))/ss          &
		&  + Bc4(j)*(khat(x,j)*(bhat(x,j)*s_dxBx(j)+bhat(z,j)*s_dxBz(j))  &
		&  +         khat(z,j)*(bhat(x,j)*s_dzBx(j)+bhat(z,j)*s_dzBz(j)))
	endif
enddo

do j=1,dimj
	c6(j) = 2.0*Bc1(j)*Bmag(j)*kappa(j)
enddo

do j=1,dimj
	c7(j) = 2.0*Bc2(j)*Bmag(j)*kappa(j)
enddo

do j=1,dimj
	c8(j) = Bmag(j)**2
enddo

do j=1,dimj
	if(j==1)then
		ss=sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
		c9(j) = Bmag(j)*(Bmag(j+1)-Bmag(j))/ss
	elseif(j==dimj)then
		ss=sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2)
		c9(j) = Bmag(j)*(Bmag(j)-Bmag(j-1))/ss
	else
		ss = sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
		&  + sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
		c9(j) = Bmag(j)*(Bmag(j+1)-Bmag(j-1))/ss
	endif
enddo

do j=1,dimj
	if(j==1)then
		ss=sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
		c10(j) = 2.0*Bmag(j)*Bc3(j)*kappa(j)-Bc4(j)**2 &
		&      + Bmag(j)*(Bc4(j+1)-Bc4(j))/ss
	elseif(j==dimj)then
		ss=sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2)
		c10(j) = 2.0*Bmag(j)*Bc3(j)*kappa(j)-Bc4(j)**2 &
		&      + Bmag(j)*(Bc4(j)-Bc4(j-1))/ss
	else
		ss = sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
		&  + sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
		c10(j) = 2.0*Bmag(j)*Bc3(j)*kappa(j)-Bc4(j)**2 &
		&      + Bmag(j)*(Bc4(j+1)-Bc4(j-1))/ss
	endif
enddo

fop = 0.0

!upper left quadrant
do j=2,dimj-1
	ss = 0.5*sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
	&  + 0.5*sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
	fop(j,j+1) = c1(j)/(ss**2)+c2(j)/(2.0*ss)
	fop(j,j) = c3(j)-2.0*c1(j)/(ss**2)
	fop(j,j-1) = c1(j)/(ss**2)-c2(j)/(2.0*ss)
enddo

!upper right quadrant
do j=2,dimj-1
	ss = 0.5*sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
	&  + 0.5*sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
	fop(j,j+1+dimj) = c4(j)/(2.0*ss)
	fop(j,j+dimj) = c5(j)
	fop(j,j-1+dimj) = -c4(j)/(2.0*ss)
enddo

!lower left quadrant
do j=2,dimj-1
	ss = 0.5*sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
	&  + 0.5*sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
	fop(j+dimj,j+1) =  c6(j)/(2.0*ss)
	fop(j+dimj,j) = c7(j)
	fop(j+dimj,j-1) = -c6(j)/(2.0*ss)
enddo

!lower right quadrant
do j=2,dimj-1
	ss = 0.5*sqrt((state(x,j)-state(x,j-1))**2+(state(z,j)-state(z,j-1))**2) &
	&  + 0.5*sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2)
	fop(j+dimj,j+1+dimj) = c8(j)/(ss**2)+c9(j)/(2.0*ss)
	fop(j+dimj,j+dimj) = c10(j)-2.0*c8(j)/(ss**2)
	fop(j+dimj,j-1+dimj) = c8(j)/(ss**2)-c9(j)/(2.0*ss)
enddo

open(unit=111, file=trim(datapath)//'.fop', status="replace")
do j = 2,dimj-1
	write(111,*) fop(2:dimj-1,j),fop(dimj+2:2*dimj-1,j)
enddo
do j = dimj+2,2*dimj-1
	write(111,*) fop(2:dimj-1,j),fop(dimj+2:2*dimj-1,j)
enddo
close(111)




deallocate(khat)
deallocate(bhat)

deallocate(kappa)
deallocate(Bmag)
deallocate(ck)
deallocate(cs)
deallocate(ca)
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

end subroutine fileigen