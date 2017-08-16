!******************************************************
! Diagnostic routines
!******************************************************

subroutine loadState
	use stateMod
	implicit none
	integer :: i
	open(unit=11,  file="state.temp",  status="old")
	do i = 1,dimj
		read(11,*) state(:,i)
	end do
	close(11)
	open(unit=11,  file="newState.temp",  status="old")
	do i = 1,dimj
		read(11,*) newState(:,i)
	end do
	close(11)
	open(unit=11,  file="oldState.temp",  status="old")
	do i = 1,dimj
		read(11,*) oldState(:,i)
	end do
	close(11)
	write(6,*) 'temp files loaded'
end subroutine loadState

subroutine dumpState
	use stateMod
	implicit none
	integer :: i
	open(unit=11,file="state.temp",  status="replace")
	do i = 1,dimj
		write(11,*) state(:,i)
	end do
	close(11)
	open(unit=11,file="newState.temp",  status="replace")
	do i = 1,dimj
		write(11,*) newState(:,i)
	end do
	close(11)
	open(unit=11,file="oldState.temp",  status="replace")
	do i = 1,dimj
		write(11,*) oldState(:,i)
	end do
	close(11)
	write(6,*) 'temp files written'
end subroutine dumpState

subroutine testState
	use stateMod
	use stateIndex
	use initCond
	use parametersBackground
	use constants
	use filConstants
	use potential
	use cont
	use errorMod
	implicit none
	real :: pBack,Bx,Bz
	real :: l,ptotal,ptot
	real :: x0,x1,z0,z1
	real :: dxdz
	real :: Vbg1,Vbg2,Kbg1,Kbg2,filK
	integer :: i,j,k
	real :: Aback
	real :: c0,c1,c2,u,a,b0,c00,c11,c22,c33,c44
	character(len=100) :: fileName
	logical :: notDone
	integer :: bufsize
	real,dimension(x:z,dimj) :: pos
	logical :: isopen
	integer :: ilast,ifirst
	real :: dipx,dipy,dipz
	real :: bilinearInterp
	real,allocatable,dimension(:) :: apexi,apex_z,springx,springz,filBetai,filMassi,midPres,midBz,midV
	real :: dk
	real :: filBeta
	real :: timeSlice
	real :: upscalingFactor
	real :: dxBx,dxBz,dzBx,dzBz,dxPtot,dzPtot
	real :: massFil,filVol
	real,allocatable,dimension(:) :: accxx,accxz,acczz,acczx,springxpm,springzpm,springMass
	real,dimension(x:z,dimj) :: acc1,acc2
	real :: springH
	integer :: potGrid,mid1,mid2,mid3,mid4
	real,allocatable,dimension(:) :: potposx,potposz
	real,allocatable,dimension(:,:) :: accxGrid,acczGrid
	real :: apexPotx0,apexPotz0,PotAmp

	j = dimj/2+1
	
	l = 0.5*sqrt((state(x,j+1)-state(x,j))**2+(state(z,j+1)-state(z,j))**2) &
	& + 0.5*sqrt((state(x,j-1)-state(x,j))**2+(state(z,j-1)-state(z,j))**2)
	
	write(*,*) Ptotal(state(x,j),state(z,j),j),&
	& state(d,j)**gamma*state(pdg,j)+1/(2.0*mu0)*(state(d,j)*l*Re/state(mf,j))**2,&
	& state(p,j) + state(b,j)**2/(2.0*mu0)
	
	write(*,*) state(d,j)*l*Re/state(mf,j),state(b,j)
	write(*,*) state(d,j)**gamma*state(pdg,j),state(p,j)
	write(*,*) state(d,j)*l*Re/state(b,j), state(mf,j)
	write(*,*) state(p,j)/state(d,j)**gamma, state(pdg,j)
	write(*,*) state(t,j),state(p,j)/(state(d,j)*kb)
	
end subroutine
