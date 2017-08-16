
!********************************* numerical methods *****************************************************

subroutine adaptive(numMethod,derivs,algebra,numMethodp,derivsp,algebrap,id) ! dynamic time step
!checks the CLF condition for the configuration of the mass elements by finding the fast mode velocity (v alfven)
! from the internal parameters of the filament, and adjusts the time step accordingly
	use stateIndex
	use stateMod
	use constants
	use errorMod
	implicit none
	external :: numMethod,derivs,algebra,numMethodp,derivsp,algebrap
	real,dimension(dimj) :: vA
	real,dimension(dimj) :: l
	real :: c0,c1
	integer :: i,id
	logical::plasStopL
	
	!$OMP master
	do i = 1,dimj
		c0 = state(d,i)*mi*mu0*cmC/nano
		vA(i) = cm*state(b,i)/(nano*sqrt(c0))
	end do

	i = 1
	l(i) = Re*sqrt((state(x,i)-state(x,i+1))**2+(state(z,i)-state(z,i+1))**2)
	do i = 2,dimj-1
		c0 = Re*sqrt((state(x,i)-state(x,i+1))**2+(state(z,i)-state(z,i+1))**2)
		c1 = Re*sqrt((state(x,i)-state(x,i-1))**2+(state(z,i)-state(z,i-1))**2)
		l(i) = min(c0,c1)
	end do
	i = dimj
	l(i) =  Re*sqrt((state(x,i)-state(x,i-1))**2+(state(z,i)-state(z,i-1))**2)
	
	
	if(.false. .and. dimjp /= 0)then
		do i = 1,dimjp
			if(i == 1)then
				c0 = Re*sqrt((statePlas(x,i)-statePlas(x,i+1))**2+(statePlas(z,i)-statePlas(z,i+1))**2)
				c1 = Re*sqrt((statePlas(x,i)-statePlas(x,dimjp))**2+(statePlas(z,i)-statePlas(z,dimjp))**2)
				l(i) = min(c0,c1)
			elseif(i==dimjp)then
				c0 = Re*sqrt((statePlas(x,i)-statePlas(x,1))**2+(statePlas(z,i)-statePlas(z,1))**2)
				c1 = Re*sqrt((statePlas(x,i)-statePlas(x,i-1))**2+(statePlas(z,i)-statePlas(z,i-1))**2)
				l(i) = min(c0,c1)
			else
				c0 = Re*sqrt((statePlas(x,i)-statePlas(x,i+1))**2+(statePlas(z,i)-statePlas(z,i+1))**2)
				c1 = Re*sqrt((statePlas(x,i)-statePlas(x,i-1))**2+(statePlas(z,i)-statePlas(z,i-1))**2)
				l(i) = min(c0,c1)
			endif
		end do	
	endif
	
	c0 = l(1)/vA(1)
	do i = 2,dimj
		c0 = min(l(i)/vA(i),c0)
	end do
	
	tau = c0
	
	tau = max(tau,tau_min_limit)
	!$OMP end master
		!$OMP barrier

	call numMethod(state,oldState,time,tau,newState,derivs,algebra,id)
		!$OMP barrier

	!$OMP master
	call smoothing(newState,state,vx,newState)
	call smoothing(newState,state,vz,newState)
	
	if(.false. .and. dimjp/=0)then
		call numMethodp(statePlas,oldStatePlas,time,tau,newStatePlas,derivsp,algebrap,id)
		call smoothingp(newStatePlas,statePlas,vx,newStatePlas)
		call smoothingp(newStatePlas,statePlas,vz,newStatePlas)
	endif
	!$OMP end master
	
end subroutine adaptive

subroutine fixedtau(numMethod,derivs,algebra,numMethodp,derivsp,algebrap,id) 
! calls the smoothing subroutines, fixed time step
	use stateIndex
	use stateMod
	implicit none
	external :: numMethod,derivs,algebra,numMethodp,derivsp,algebrap
	integer :: id
	
	call numMethod(state,oldState,time,tau,newState,derivs,algebra,id)
		!$OMP barrier
	
	!$OMP master
		call smoothing(newState,state,vx,newState)
		call smoothing(newState,state,vz,newState)
	!$OMP end master

	!$OMP master
	if(dimjp/=0)then
		call numMethodp(statePlas,oldStatePlas,time,tau,newStatePlas,derivsp,algebrap)
		call smoothingp(newStatePlas,statePlas,vx,newStatePlas)
		call smoothingp(newStatePlas,statePlas,vz,newStatePlas)
	endif
	!$OMP end master
	
end subroutine fixedtau

subroutine smoothing(s,so,i,sn) 
! smooting term added to central difference equations from
! Hoffmann & Chiang with smoothed boundary points
	use stateIndex
	use dimM
	implicit none
	real :: wFactor = 0.1
	real :: smoothFactor = 4.0 !4.0
	real,dimension(dimi,dimj) :: s,sn,so
	integer :: i,j
	sn = s
	do j = 3,dimj-2
		sn(i,j) = s(i,j)-wFactor*(so(i,j-2)-4.0*so(i,j-1)+6.0*so(i,j)-4.0*so(i,j+1)+so(i,j+2))
	end do
	j = 2
	sn(i,j) = (s(i,j-1)+smoothFactor*s(i,j)+s(i,j+1))/(2.0+smoothFactor)
	j = dimj-1
	sn(i,j) = (s(i,j-1)+smoothFactor*s(i,j)+s(i,j+1))/(2.0+smoothFactor)
end subroutine smoothing

subroutine euler(s,so,time,tau,sn,derivs,algebra,id)
	use dimM
	use stateIndex
	implicit none
	integer :: id
	external :: derivs,algebra
	real,dimension(dimi,dimj) :: s,so,sn,f
	real :: time,tau
	integer :: jmin,jmax
	jmin = thr_domain(id+1,1)
	jmax = thr_domain(id+1,2)
	call derivs(s,time,tau,f,id)
	!$OMP barrier
	sn(:,jmin:jmax) = s(:,jmin:jmax) + tau*f(:,jmin:jmax)
	!$OMP barrier
	call algebra(sn,time+tau,tau,sn,id)
	!$OMP barrier
end subroutine euler

subroutine backEuler(s,so,time,tau,sn,derivs,algebra,id)
	use dimM
	implicit none
	integer :: id
	external :: derivs,algebra
	real,dimension(dimi,dimj) :: s,so,sn,f
	real :: time,tau
	integer :: jmin,jmax
	jmin = thr_domain(id+1,1)
	jmax = thr_domain(id+1,2)
	call derivs(s,time,tau,f,id)
	!$OMP barrier
	sn(:,jmin:jmax) = so(:,jmin:jmax) + tau*f(:,jmin:jmax)
	!$OMP barrier
	call algebra(sn,time+tau,tau,sn,id)
	!$OMP barrier
end subroutine backEuler

subroutine trapezoidal(s,so,time,tau,sn,derivs,algebra,id)
	use dimM
	implicit none
	integer :: id
	external :: derivs,algebra
	real,dimension(dimi,dimj) :: s,so,sn,f,fo
	real :: time,tau
	integer :: jmin,jmax
	jmin = thr_domain(id+1,1)
	jmax = thr_domain(id+1,2)
	call derivs(s,time,tau,f,id)
	!$OMP barrier
	call derivs(so,time,tau,fo,id)
	!$OMP barrier
	sn(:,jmin:jmax) = so(:,jmin:jmax) + tau*0.5*(f(:,jmin:jmax) + fo(:,jmin:jmax))
	!$OMP barrier
	call algebra(sn,time+tau,tau,sn,id)
	!$OMP barrier
end subroutine trapezoidal

subroutine rk2(s,so,time,tau,sn,derivs,algebra,id)
	use dimM
	use rkVectors
	implicit none
	integer :: id
	external :: derivs,algebra
	real,dimension(dimi,dimj) :: s,so,sn,f
	real :: time,tau
	integer :: jmin,jmax
	jmin = thr_domain(id+1,1)
	jmax = thr_domain(id+1,2)
	call derivs(s,time,tau,f,id)
	!$OMP barrier
	k1(:,jmin:jmax) = tau*f(:,jmin:jmax)
	s1(:,jmin:jmax) = s(:,jmin:jmax)+k1(:,jmin:jmax)
	!$OMP barrier
	call algebra(s1,time+tau,tau,sa,id)
	!$OMP barrier
	call derivs(sa,time+tau,tau,f,id)
	!$OMP barrier
	k2(:,jmin:jmax) = tau*f(:,jmin:jmax)
	!$OMP barrier
	sn(:,jmin:jmax) = s(:,jmin:jmax) + 0.5*(k1(:,jmin:jmax)+k2(:,jmin:jmax))
	!$OMP barrier
	call algebra(sn,time+tau,tau,sn,id)
	!$OMP barrier
end subroutine rk2

subroutine rk4(s,so,time,tau,sn,derivs,algebra,id)
	use dimM
	use rkVectors
	implicit none
	integer :: id
	external :: derivs,algebra
	real,dimension(dimi,dimj) :: s,so,sn,f
	real :: time,tau
	integer :: jmin,jmax
	jmin = thr_domain(id+1,1)
	jmax = thr_domain(id+1,2)
	call derivs(s,time,tau,f,id)
	!$OMP barrier
	k1(:,jmin:jmax) = tau*f(:,jmin:jmax)
	s1(:,jmin:jmax) = s(:,jmin:jmax) + 0.5*k1(:,jmin:jmax)
	!$OMP barrier
	call algebra(s1,time+0.5*tau,tau,sa,id)
	!$OMP barrier

	call derivs(sa,time+0.5*tau,tau,f,id)
	!$OMP barrier
	k2(:,jmin:jmax) = tau*f(:,jmin:jmax)
	s2(:,jmin:jmax) = s(:,jmin:jmax) + 0.5*k2(:,jmin:jmax)
	!$OMP barrier
	call algebra(s2,time+0.5*tau ,tau,sa,id)
	!$OMP barrier

	call derivs(sa,time+0.5*tau,tau,f,id)
	!$OMP barrier
	k3(:,jmin:jmax) = tau*f(:,jmin:jmax)
	s3(:,jmin:jmax) = s(:,jmin:jmax) + k3(:,jmin:jmax)
	!$OMP barrier
	call algebra(s3,time+tau,tau,sa,id)
	!$OMP barrier
	
	call derivs(sa,time+tau,tau,f,id)
	!$OMP barrier
	k4(:,jmin:jmax) = tau*f(:,jmin:jmax)
	!$OMP barrier
	sn(:,jmin:jmax) = s(:,jmin:jmax) &
	& + (k1(:,jmin:jmax) + 2.0*k2(:,jmin:jmax) + 2.0*k3(:,jmin:jmax) + k4(:,jmin:jmax))/6.0
	!$OMP barrier
	call algebra(sn,time+tau,tau,sn,id)
	!$OMP barrier
end subroutine rk4
