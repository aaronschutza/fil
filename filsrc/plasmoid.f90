!************************ Plasmoid functions ************************
! These funcitons operate separately on the reconnected plasmoid
!********************************************************************

subroutine reconnect
	use stateMod
	use stateIndex
	use errorMod
	implicit none
	integer :: i,j
	real :: z1,z2,z3,z4
	real :: x1,x2,x3,x4
	logical :: noApex
	logical :: reconl
	logical :: firstPlas
	integer,dimension(dimj) :: reconpoints
	integer :: added
	reconpoints = 0
	j = 1
	noApex = .true.
	reconl = .false.
	do while(noApex .and. j < dimj)
		z1 = state(z,j)
		z2 = state(z,j+1)
		z3 = state(z,dimj-j)
		z4 = state(z,dimj+1-j)
		
		if (z1 > z4 .and. z2 > z3) then
			!no recon
		else if (z4 > z1 .and. z3 > z2) then
			!no recon
		else if (z1 > z4 .and. z3 > z2) then
			!recon
			reconl = .true.
			reconpoints(j) = 1
			write(*,*) 'recon point found'
		else if (z4 > z1 .and. z2 > z3) then
			!recon
			reconl = .true.
			reconpoints(j) = 1
			write(*,*) 'recon point found'
		else if (z1 > z4 .and. z2 == z3) then
			!apex
			noApex = .false.
			if (reconl) then
				write(*,*) 'apex found'
			end if
		else if (z4 > z1 .and. z2 == z3) then
			!apex
			noApex = .false.
			if (reconl) then
				write(*,*) 'apex found'
			end if
		else
			!no recon
		end if
		
		j = j+1
		
	end do
	
	if (reconl) then
		
		do j = 1,dimj
			if (reconpoints(j)==1) then
				i = j
				write(*,*) 'cutting at point ',i
				exit
			end if
		end do
		x1 = state(x,i)
		x2 = state(x,i+1)
		x3 = state(x,dimj-i)
		x4 = state(x,dimj+1-i)
		
		z1 = state(z,i)
		z2 = state(z,i+1)
		z3 = state(z,dimj-i)
		z4 = state(z,dimj+1-i)
		
		added = dimj-2*i - 1
		dimjp = dimjp + added
		if(.not. allocated(statePlas))then
			allocate(statePlas(dimi,dimjp))
			allocate(statePlasTemp(dimi,dimjp))
			allocate(newStatePlas(dimi,dimj))
			allocate(oldStatePlas(dimi,dimj))
			firstPlas = .true.
		else
			statePlasTemp = statePlas
			deallocate(statePlas)
			allocate(statePlas(dimi,dimjp))
			statePlas(:,(added/2)+1:(dimjp-added/2)) = statePlasTemp
			deallocate(statePlasTemp)
			allocate(statePlasTemp(dimi,dimjp))
			deallocate(newStatePlas)
			deallocate(oldStatePlas)
			allocate(newStatePlas(dimi,dimj))
			allocate(oldStatePlas(dimi,dimj))
			firstPlas = .false.
		end if
		
		state(x,i+1) = (-x3 * z1 * x2 + x4 * z1 * x2 - x4 * z3 * x2 + x4 * z3 * x1 &
					&- x4 * z2 * x1 + z4 * x3 * x2 - z4 * x3 * x1 + x3 * z2 * x1) &
						& / (x3* z2 - x3 * z1 - x4 * z2 + x4 * z1 - z3 * x2 + z3 * x1 + z4 * x2 - z4 * x1)
		state(z,i+1) = (-x4 * z3 * z2 + x4 * z3 * z1 + z4 * x3 * z2 - z4 * x3 * z1 &
					&+ z3 * x1 * z2 - z3 * z1 * x2 - z4 * x1 * z2 + z4 * z1 * x2) &
						& / (x3 * z2 - x3 * z1 - x4 * z2 + x4 * z1 - z3 * x2 + z3 * x1 + z4 * x2 - z4 * x1)
		
		state(x,dimj-i) = (-x3 * z1 * x2 + x4 * z1 * x2 - x4 * z3 * x2 + x4 * z3 * x1 &
					& - x4 * z2 * x1 + z4 * x3 * x2 - z4 * x3 * x1 + x3 * z2 * x1) &
						& / (x3* z2 - x3 * z1 - x4 * z2 + x4 * z1 - z3 * x2 + z3 * x1 + z4 * x2 - z4 * x1)
		state(z,dimj-i) = (-x4 * z3 * z2 + x4 * z3 * z1 + z4 * x3 * z2 - z4 * x3 * z1&
					&+ z3 * x1 * z2 - z3 * z1 * x2 - z4 * x1 * z2 + z4 * z1 * x2) &
						& / (x3 * z2 - x3 * z1 - x4 * z2 + x4 * z1 - z3 * x2 + z3 * x1 + z4 * x2 - z4 * x1)
		
		state(vx,i+1) = (state(vx,i+1)+state(vx,dimj-i))/2.0
		state(vz,i+1) = (state(vz,i+1)+state(vz,dimj-i))/2.0
		
		state(vx,dimj-i) = (state(vx,dimj-i)+state(vx,i+1))/2.0
		state(vz,dimj-i) = (state(vz,dimj-i)+state(vz,i+1))/2.0
		
		stateTemp = state
		
		if(firstPlas)then
			statePlas = state(:,(i+2):(dimj-i))
		else
			statePlas(:,1:(added/2)) = state(:,(i+2):(dimj/2+1))
			statePlas(:,(dimjp-added/2+1):dimjp) = state(:,(dimj/2+2):(dimj-i))
		endif
		
		deallocate(state)
		allocate(state(dimi,dimj-added))
		
		state(:,1:(i+1)) = stateTemp(:,1:(i+1))
		state(:,(dimj-added-i+1):(dimj-added)) = stateTemp(:,(dimj-i+1):(dimj))
		
		dimj = dimj - added
		deallocate(stateTemp)
		allocate(stateTemp(dimi,dimj))
		deallocate(newState)
		allocate(newState(dimi,dimj))
		deallocate(oldState)
		allocate(oldState(dimi,dimj))
		write(*,*) 'added',added
		write(*,*) 'dimj',dimj
		write(*,*) 'dimjp',dimjp
		write(*,*) 'dim',dimj+dimjp
		write(*,*) ' '
		
		
		!if(firstPlas)then
!			open(unit=112, file="reconTest", status="replace")
!			close(112)
!		endif
!		
!		open(unit=112, file="reconTest",position="append", status="old")
!			do i = 1,dimj
!				write(112,*) state(x:vz,i)
!			end do
!			write(112,*) ' '
!			write(112,*) ' '
!			if (dimjp /= 0) then
!				do i = 1,dimjp
!					write(112,*) statePlas(x:vz,i)
!				end do
!				write(112,*) ' '
!				write(112,*) ' '
!			end if
!		close(112)
		
		
		errorFlag2 = .true.
	end if
end subroutine reconnect

subroutine reconnect2
	use stateMod
	use stateIndex
	use errorMod
	implicit none
	integer :: i,j
	real :: z1,z2,z3,z4
	real :: x1,x2,x3,x4
	logical :: noApex
	logical :: reconl
	logical :: firstPlas
	integer,dimension(dimj) :: reconpoints
	integer :: added
	reconpoints = 0
	j = 1
	noApex = .true.
	reconl = .false.
	do while(noApex .and. j < dimj)
		z1 = state(z,j)
		z2 = state(z,j+1)
		z3 = state(z,dimj-j)
		z4 = state(z,dimj+1-j)
		
		if (z1 > z4 .and. z2 > z3) then
			!no recon
		else if (z4 > z1 .and. z3 > z2) then
			!no recon
		else if (z1 > z4 .and. z3 > z2) then
			!recon
			reconl = .true.
			reconpoints(j) = 1
			write(*,*) 'recon point found'
		else if (z4 > z1 .and. z2 > z3) then
			!recon
			reconl = .true.
			reconpoints(j) = 1
			write(*,*) 'recon point found'
		else if (z1 > z4 .and. z2 == z3) then
			!apex
			noApex = .false.
			if (reconl) then
				write(*,*) 'apex found'
			end if
		else if (z4 > z1 .and. z2 == z3) then
			!apex
			noApex = .false.
			if (reconl) then
				write(*,*) 'apex found'
			end if
		else
			!no recon
		end if
		
		j = j+1
		
	end do
	
	if (reconl) then
		
		do j = 1,dimj
			if (reconpoints(j)==1) then
				i = j
				write(*,*) 'cutting at point ',i
				exit
			end if
		end do
		
		
		added = dimj - 2*i
		dimjp = dimjp + added
		
		write(*,*) 'allocating'
		
		allocate(statePlas(dimi,dimjp))
		allocate(statePlasTemp(dimi,dimjp))
		allocate(newStatePlas(dimi,dimjp))
		allocate(oldStatePlas(dimi,dimjp))
		
		write(*,*) 'configuring'
		
		stateTemp = state
		

		statePlas = state(:,(i+1):(dimj-i))
		newStatePlas = statePlas
		statePlasTemp = statePlas
		oldStatePlas = statePlas
		
		deallocate(state)
		allocate(state(dimi,2*i))
		
		state(:,1:i) = stateTemp(:,1:i)
		state(:,(i+1):(2*i)) = stateTemp(:,(dimj-i+1):(dimj))
		
		dimj = 2*i
		deallocate(stateTemp)
		allocate(stateTemp(dimi,dimj))
		deallocate(newState)
		allocate(newState(dimi,dimj))
		deallocate(oldState)
		allocate(oldState(dimi,dimj))
		write(*,*) 'added',added
		write(*,*) 'dimj',dimj
		write(*,*) 'dimjp',dimjp
		write(*,*) 'dim',dimj+dimjp
		write(*,*) ' '
		
		stateTemp = state
		newState = state
		oldState = state
		
		
		!~ write(*,*) 'state ='
		!~ do j = 1,dimj
			!~ write(*,*) state(x:vz,j)
		!~ enddo
		!~ write(*,*) 'statep ='
		!~ do j = 1,dimjp
			!~ write(*,*) statePlas(x:vz,j)
		!~ enddo
		
		!if(firstPlas)then
!			open(unit=112, file="reconTest", status="replace")
!			close(112)
!		endif
!		
!		open(unit=112, file="reconTest",position="append", status="old")
!			do i = 1,dimj
!				write(112,*) state(x:vz,i)
!			end do
!			write(112,*) ' '
!			write(112,*) ' '
!			if (dimjp /= 0) then
!				do i = 1,dimjp
!					write(112,*) statePlas(x:vz,i)
!				end do
!				write(112,*) ' '
!				write(112,*) ' '
!			end if
!		close(112)
		
		errorFlag2 = .true.
		write(*,*) 'recon turned off'
	end if
end subroutine reconnect2

subroutine derivs1p(s,time,tau,f) ! updates derivatives of attributes given a state S
	use dimM
	use stateIndex
	use boundMod
	use errorMod
	implicit none
	real,dimension(dimi,dimjp) :: s,f
	real,dimension(vx:vz,dimjp) :: acc
	real :: time,tau
	f(x,:)     = s(vx,:)
	f(z,:)     = s(vz,:)
	call accelp(s,acc)
	f(vx:vz,:) = acc
	f(d:t,:)   = 0.0
	f(b,:)     = 0.0
	f(pdg,:)   = 0.0
	f(mf,:)     = 0.0
	
	return
end subroutine derivs1p

subroutine algebra1p(s,time,tau,sn) ! updates attributes via algebraic equations
	use dimM
	use stateIndex
	use constants
	use initCond
	use errorMod
	use variables1,only: n
	implicit none
	real,dimension(dimi,dimjp) :: s
	real,dimension(dimi,dimjp) :: sn
	real :: time,tau,Bz,Bx,pBack,k,pb
	real,dimension(dimjp) :: Rho
	integer :: i
	real :: l
	sn = s
	! find the number density Rho and everything else is updated
	call findRho1p(s,Rho)
	
	sn(d,:) = Rho ! number density
	
	if(errorFlag)then
		write(6,*) 'n=',n,'t=',time
	end if
	
	do i = 1,dimjp
		if(Rho(i)<0.0)then
			write(6,*) 'Rho is negative at index=',i
		end if
		sn(p,i) = s(pdg,i)*(Rho(i)**gamma) ! P = P_init*(Rho_init)^-gamma * (Rho)^gamma
	end do
	
	sn(t,:) = s(p,:)/( kb*s(d,:) ) ! P = n*kb*t
	
	
	do i = 1,dimjp
		!new way of finding Bfil with simpler mass conservation formula (analytically the same as other method)
		if (i == dimjp) then
			l = 0.5*sqrt( (s(x,1)-s(x,i))**2 + (s(z,1)-s(z,i))**2 ) &
			& + 0.5*sqrt( (s(x,i-1)-s(x,i))**2 + (s(z,i-1)-s(z,i))**2 )
			s(b,i) = s(d,i)*Re*l/s(mf,i)
		else if (i == 1) then
			l = 0.5*sqrt( (s(x,i+1)-s(x,i))**2 + (s(z,i+1)-s(z,i))**2 ) &
			& + 0.5*sqrt( (s(x,dimjp)-s(x,i))**2 + (s(z,dimjp)-s(z,i))**2 )
			s(b,i) = s(d,i)*Re*l/s(mf,i)
		else
			l = 0.5*sqrt( (s(x,i+1)-s(x,i))**2 + (s(z,i+1)-s(z,i))**2 ) &
			& + 0.5*sqrt( (s(x,i-1)-s(x,i))**2 + (s(z,i-1)-s(z,i))**2 )
			s(b,i) = s(d,i)*Re*l/s(mf,i)
		end if
	end do
	
end subroutine algebra1p

subroutine smoothingp(s,so,i,sn) 
! smooting term added to central difference equations from
! Hoffmann & Chiang with smoothed boundary points
	use stateIndex
	use dimM
	implicit none
	real :: wFactor = 0.1
	real :: smoothFactor = 4.0 !4.0
	real,dimension(dimi,dimjp) :: s,sn,so
	integer :: i,j
	sn = s
	
	do j = 1,dimj
		if(j==1)then
			sn(i,j) = s(i,j) - wFactor*( so(i,dimjp-1) - 4.0*so(i,dimjp) + 6.0*so(i,j) - 4.0*so(i,j+1) + so(i,j+2) )
		elseif(j==2)then
			sn(i,j) = s(i,j) - wFactor*( so(i,dimjp) - 4.0*so(i,j-1) + 6.0*so(i,j) - 4.0*so(i,j+1) + so(i,j+2) )
		elseif(j==dimjp-1)then
			sn(i,j) = s(i,j) - wFactor*( so(i,j-2) - 4.0*so(i,j-1) + 6.0*so(i,j) - 4.0*so(i,j+1) + so(i,1) )
		elseif(j==dimjp)then
			sn(i,j) = s(i,j) - wFactor*( so(i,j-2) - 4.0*so(i,j-1) + 6.0*so(i,j) - 4.0*so(i,1) + so(i,2) )
		else
			sn(i,j) = s(i,j) - wFactor*( so(i,j-2) - 4.0*so(i,j-1) + 6.0*so(i,j) - 4.0*so(i,j+1) + so(i,j+2) )
		endif
	end do
	
end subroutine smoothingp
	
subroutine findRho1p(s,Rho) ! 1st attept to update number density using algebraic apporach
	use dimM
	use stateIndex
	use constants
	use errorMod
	implicit none
	real,dimension(dimi,dimjp),intent(in) :: s
	real,dimension(dimjp) :: Rho
	integer :: i,j,k
	integer :: point
	real :: coefA,coefB,l,coef,Bx,Bz,pBack
	do i = 1,dimjp
		if (i == 1) then
			l = 0.5*sqrt((s(x,i+1)-s(x,i))**2+(s(z,i+1)-s(z,i))**2) &
			& + 0.5*sqrt((s(x,dimjp)-s(x,i))**2+(s(z,dimjp)-s(z,i))**2)
			coef = 2.0*mu0*s(mf,i)*s(mf,i)/(l*l*Re*Re)
			coefA = coef*s(pdg,i)
			coefB = coef*(pBack(s(x,i),s(z,i),i)+((Bx(s(x,i),s(z,i),i))**2+(Bz(s(x,i),s(z,i),i))**2)/(2.0*mu0))
			call newtonRho2(coefA,coefB,s(d,i),Rho(i),i)
		else if (i == dimjp) then
			l = 0.5*sqrt((s(x,1)-s(x,i))**2+(s(z,1)-s(z,i))**2) &
			& + 0.5*sqrt((s(x,i-1)-s(x,i))**2+(s(z,i-1)-s(z,i))**2)
			coef = 2.0*mu0*s(mf,i)*s(mf,i)/(l*l*Re*Re)
			coefA = coef*s(pdg,i)
			coefB = coef*(pBack(s(x,i),s(z,i),i)+((Bx(s(x,i),s(z,i),i))**2+(Bz(s(x,i),s(z,i),i))**2)/(2.0*mu0))
			call newtonRho2(coefA,coefB,s(d,i),Rho(i),i)
		else
			l = 0.5*sqrt((s(x,i+1)-s(x,i))**2+(s(z,i+1)-s(z,i))**2) &
			& + 0.5*sqrt((s(x,i-1)-s(x,i))**2+(s(z,i-1)-s(z,i))**2)
			coef = 2.0*mu0*s(mf,i)*s(mf,i)/(l*l*Re*Re)
			coefA = coef*s(pdg,i)
			coefB = coef*(pBack(s(x,i),s(z,i),i)+((Bx(s(x,i),s(z,i),i))**2+(Bz(s(x,i),s(z,i),i))**2)/(2.0*mu0))
			call newtonRho2(coefA,coefB,s(d,i),Rho(i),i)
		end if
	end do
end subroutine findRho1p

subroutine eulerp(s,so,time,tau,sn,derivs,algebra)
	use dimM
	implicit none
	external :: derivs,algebra
	real,dimension(dimi,dimjp) :: s,so,sn,f
	real :: time,tau
	call derivs(s,time,tau,f)
	sn = s + tau*f
	call algebra(sn,time,tau,sn)
end subroutine eulerp

subroutine backEulerp(s,so,time,tau,sn,derivs,algebra)
	use dimM
	implicit none
	external :: derivs,algebra
	real,dimension(dimi,dimjp) :: s,so,sn,f
	real :: time,tau
	call derivs(s,time,tau,f)
	sn = so + tau*f
	call algebra(sn,time+tau,tau,sn)
end subroutine backEulerp

subroutine trapezoidalp(s,so,time,tau,sn,derivs,algebra)
	use dimM
	implicit none
	external :: derivs,algebra
	real,dimension(dimi,dimjp) :: s,so,sn,f,fo
	real :: time,tau
	call derivs(s,time,tau,f)
	call derivs(so,time,tau,fo)
	sn = so + tau*0.5*(f + fo)
	call algebra(sn,time+tau,tau,sn)
end subroutine trapezoidalp

subroutine rk2p(s,so,time,tau,sn,derivs,algebra)
	use dimM
	implicit none
	external :: derivs,algebra
	real,dimension(dimi,dimjp) :: s,so,sn,sa,f,k1,k2
	real :: time,tau
	call derivs(s,time,tau,f)
	k1 = tau*f
	call algebra( s+k1, time+tau,tau, sa)
	call derivs(    sa, time+tau,tau,  f)
	k2 = tau*f
	sn = s + 0.5*(k1 + k2)
	call algebra(sn,time+tau,tau,sn)
end subroutine rk2p

subroutine rk4p(s,so,time,tau,sn,derivs,algebra)
	use dimM
	use errorMod
	implicit none
	external :: derivs,algebra
	real,dimension(dimi,dimjp) :: s,so,sn,sa,f,k1,k2,k3,k4
	real :: time,tau
	
	call derivs(s,time,tau,f)
	k1 = tau*f
	call algebra(s+0.5*k1, time+0.5*tau ,tau, sa)
	call derivs(       sa, time+0.5*tau ,tau,  f)
	k2 = tau*f
	call algebra(s+0.5*k2, time+0.5*tau ,tau, sa)
	call derivs(       sa, time+0.5*tau ,tau,  f)
	k3 = tau*f
	call algebra(    s+k3, time+tau     ,tau, sa)
	call derivs(       sa, time+tau     ,tau,  f)
	k4 = tau*f
	sn = s + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0
	call algebra(sn,time+tau,tau,sn)
end subroutine rk4p

subroutine accelp(s,acc) ! accel of mass elements in Re/s^2. Central differnce methods used for spatial derivatives
	use dimM
	use stateIndex
	use constants
	implicit none
	real,dimension(dimi,dimjp) :: s
	real,dimension(vx:vz,dimjp) :: acc,a1,a2,a3
	real,dimension(dimjp) :: Bfil,Bfil_x,Bfil_z,Bmed_x,Bmed_z,dxdl,dzdl
	integer :: i
	real :: l,dx1,dz1,dx2,dz2,dx3,dz3,Bx,Bz,pBack
	
	! Bmed is the background medium field
	do i = 1,dimjp
		Bmed_x(i) = Bx(s(x,i),s(z,i),i)
		Bmed_z(i) = Bz(s(x,i),s(z,i),i)
	end do
	! magnitude of the filament magnetic field using pressure balance
	Bfil(:) = s(b,:)
	! one sided difference for boundaries ******************
	i = 1
	l = sqrt( (s(x,i+1) - s(x,i)  )**2 + (s(z,i+1) - s(z,i)  )**2 ) &
		& + sqrt( (s(x,i)   - s(x,dimjp))**2 + (s(z,i)   - s(z,dimjp))**2 )
		dxdl(i) = (s(x,i+1) - s(x,dimjp))/l
		dzdl(i) = (s(z,i+1) - s(z,dimjp))/l
	
	! second boundary **************************************
	i = dimjp
	l = sqrt( (s(x,1) - s(x,i)  )**2 + (s(z,1) - s(z,i)  )**2 ) &
		& + sqrt( (s(x,i)   - s(x,i-1))**2 + (s(z,i)   - s(z,i-1))**2 )
		dxdl(i) = (s(x,1) - s(x,i-1))/l
		dzdl(i) = (s(z,1) - s(z,i-1))/l
	
	! central difference for dX/dl ****************************************
	do i = 2,dimjp-1
		!l = sqrt( (s(x,i+1) - s(x,i-1))**2 + (s(z,i+1) - s(z,i-1))**2 )
		l = sqrt( (s(x,i+1) - s(x,i)  )**2 + (s(z,i+1) - s(z,i)  )**2 ) &
		& + sqrt( (s(x,i)   - s(x,i-1))**2 + (s(z,i)   - s(z,i-1))**2 )
		dxdl(i) = (s(x,i+1) - s(x,i-1))/l
		dzdl(i) = (s(z,i+1) - s(z,i-1))/l
	end do
	
	Bfil_x = Bfil*dxdl
	Bfil_z = Bfil*dzdl
	
	! first term of accel equation ***************************************************************************
	a1(vx:vz,1) = 0.0; a1(vx:vz,dimjp) = 0.0
	do i = 1,dimjp
		if (i == 1) then
			l = sqrt( (s(x,i+1) - s(x,i)  )**2 + (s(z,i+1) - s(z,i)  )**2 ) &
			& + sqrt( (s(x,i)   - s(x,dimjp))**2 + (s(z,i)   - s(z,dimjp))**2 )
			a1(vx,i) = Bfil(i)/(mu0*nano*mi*cm*s(d,i)*Re*Re) * ( (Bfil_x(i+1) - Bmed_x(i+1)) - (Bfil_x(dimjp) - Bmed_x(dimjp)) )/l
			a1(vz,i) = Bfil(i)/(mu0*nano*mi*cm*s(d,i)*Re*Re) * ( (Bfil_z(i+1) - Bmed_z(i+1)) - (Bfil_z(dimjp) - Bmed_z(dimjp)) )/l
		else if (i == dimjp) then
			l = sqrt( (s(x,1) - s(x,i)  )**2 + (s(z,1) - s(z,i)  )**2 ) &
			& + sqrt( (s(x,i)   - s(x,i-1))**2 + (s(z,i)   - s(z,i-1))**2 )
			a1(vx,i) = Bfil(i)/(mu0*nano*mi*cm*s(d,i)*Re*Re) * ( (Bfil_x(1) - Bmed_x(1)) - (Bfil_x(i-1) - Bmed_x(i-1)) )/l
			a1(vz,i) = Bfil(i)/(mu0*nano*mi*cm*s(d,i)*Re*Re) * ( (Bfil_z(1) - Bmed_z(1)) - (Bfil_z(i-1) - Bmed_z(i-1)) )/l
		else
			l = sqrt( (s(x,i+1) - s(x,i)  )**2 + (s(z,i+1) - s(z,i)  )**2 ) &
			& + sqrt( (s(x,i)   - s(x,i-1))**2 + (s(z,i)   - s(z,i-1))**2 )
			a1(vx,i) = Bfil(i)/(mu0*nano*mi*cm*s(d,i)*Re*Re) * ( (Bfil_x(i+1) - Bmed_x(i+1)) - (Bfil_x(i-1) - Bmed_x(i-1)) )/l
			a1(vz,i) = Bfil(i)/(mu0*nano*mi*cm*s(d,i)*Re*Re) * ( (Bfil_z(i+1) - Bmed_z(i+1)) - (Bfil_z(i-1) - Bmed_z(i-1)) )/l
		end if
	end do
	! second term of accel equation *******************************************************************************************
	
	do i = 1,dimjp
		a2(vx,i) = (Bfil_x(i)-Bmed_x(i))*(Bx(s(x,i)+h,s(z,i),i)-Bx(s(x,i)-h,s(z,i),i))/(2.0*h*mu0*nano*mi*cm*s(d,i)*Re*Re) &
			&+ (Bfil_z(i)-Bmed_z(i))*(Bx(s(x,i),s(z,i)+h,i)-Bx(s(x,i),s(z,i)-h,i))/(2.0*h*mu0*nano*mi*cm*s(d,i)*Re*Re)
		a2(vz,i) = (Bfil_x(i)-Bmed_x(i))*(Bz(s(x,i)+h,s(z,i),i)-Bz(s(x,i)-h,s(z,i),i))/(2.0*h*mu0*nano*mi*cm*s(d,i)*Re*Re) &
			&+(Bfil_z(i)-Bmed_z(i))*(Bz(s(x,i),s(z,i)+h,i)-Bz(s(x,i),s(z,i)-h,i))/(2.0*h*mu0*nano*mi*cm*s(d,i)*Re*Re)
	end do

	acc = a1 + a2 ! in units of Re/s^2
	
end subroutine accelp

subroutine deltaDPTp(s,ddpt) ! currently the derivatives of d t and p are zero
	use dimM
	use stateIndex
	implicit none
	real,dimension(dimi,dimjp) :: s
	real,dimension(d:t,dimjp) :: ddpt
	ddpt = 0.0
	
end subroutine deltaDPTp
