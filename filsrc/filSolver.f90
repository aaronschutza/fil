
!********************************* deriv and algebra routines *****************************************

subroutine derivs1(s,time,tau,f,id) ! updates derivatives of attributes given a state S
	use dimM
	use stateIndex
	use boundMod
	use parametersBackground, only: wallBound
	implicit none
	integer :: id,jmin,jmax
	real,dimension(dimi,dimj) :: s,f
	real,dimension(vx:vz,dimj) :: acc
	real,dimension(d:t,dimj) :: ddpt
	real :: boundVz1,boundVz2,inertVz1,inertVz2
	real :: boundVz1_test,boundVz2_test
	real :: roundVz1,roundVz2,roundVx1,roundVx2
	real :: time,tau
	
	jmin = thr_domain(id+1,1)
	jmax = thr_domain(id+1,2)

	if (jmin == 1) then
		if(radialDamping)then
			s(vx,1)    = roundVx1(s,time)
			s(vz,1)    = roundVz1(s,time)
		elseif(wallBound)then
			s(vx,1)    = 0.0
			s(vz,1)    = boundVz1(s,time)
		else
			s(vx,1)    = 0.0
			s(vz,1)    = 0.0
		endif
	endif
	if (jmax == dimj) then
		if(radialDamping)then
			s(vx,dimj) = roundVx2(s,time)
			s(vz,dimj) = roundVz2(s,time)
		elseif(wallBound)then
			s(vx,dimj) = 0.0
			s(vz,dimj) = boundVz2(s,time)
		else
			s(vx,dimj) = 0.0
			s(vz,dimj) = 0.0 
		endif
	endif
	f(x,jmin:jmax)     = s(vx,jmin:jmax)
	f(z,jmin:jmax)     = s(vz,jmin:jmax)
	call accel4(s,acc,jmin,jmax)
	f(vx:vz,jmin:jmax) = acc(:,jmin:jmax)
	if(use_drag_force)then
		if(time>drag_tmin.and.time<drag_tmax)then 
			call accel_drag(s,acc,jmin,jmax)
			f(vx:vz,jmin:jmax) = f(vx:vz,jmin:jmax) + acc(:,jmin:jmax)
		endif
	endif
	if(use_driver_force)then
		if(time>drive_tmin.and.time<drive_tmax)then
			call accel_drive(s,acc,jmin,jmax,time)
			f(vx:vz,jmin:jmax) = f(vx:vz,jmin:jmax) + acc(:,jmin:jmax)
		endif
	endif
	call deltaDPT(s,ddpt,jmin,jmax)
	f(d:t,jmin:jmax)   = ddpt(:,jmin:jmax)
	f(b,jmin:jmax)     = 0.0
	f(pdg,jmin:jmax)   = 0.0
	f(mf,jmin:jmax)    = 0.0
end subroutine derivs1

subroutine algebra1(s,time,tau,sn,id) ! updates attributes via algebraic equations
	use dimM
	use stateIndex
	use constants
	use initCond
	use errorMod
	use variables1,only: n
	use boundMod
	implicit none
	real,dimension(dimi,dimj) :: s
	real,dimension(dimi,dimj) :: sn
	real :: time,tau,Bz,Bx,pBack,k,pb
	real,dimension(dimj) :: Rho
	integer :: i,id,jmin,jmax
	real :: l
	
	if(id == -1)then
		jmin = 1
		jmax = dimj
	else
		jmin = thr_domain(id+1,1)
		jmax = thr_domain(id+1,2)
	endif
	!impose radial position boundary condition
	if(radialDamping)then
		if(jmin==1)then
			l = sqrt(s(x,1)**2+s(z,1)**2)
			s(x,1) = boundary*s(x,1)/l
			s(z,1) = boundary*s(z,1)/l
		endif
		if(jmax==dimj)then
			l = sqrt(s(x,dimj)**2+s(z,dimj)**2)
			s(x,dimj) = boundary*s(x,dimj)/l
			s(z,dimj) = boundary*s(z,dimj)/l
		endif
	endif
	
	sn(:,jmin:jmax) = s(:,jmin:jmax)
	! find the number density Rho and everything else is updated
	call findRho2(s,Rho,id)
	
	sn(d,jmin:jmax) = Rho(jmin:jmax) ! number density
	
	if(errorFlag)then
		write(6,*) 'n=',n,'t=',time
	end if
	
	do i = jmin,jmax
		if(Rho(i)<0.0)then
			write(6,*) 'Rho is negative at index=',i
		end if
		sn(p,i) = s(pdg,i)*(Rho(i)**gamma) ! P = P_init*(Rho_init)^-gamma * (Rho)^gamma
	end do
	
	sn(t,jmin:jmax) = s(p,jmin:jmax)/(kb*s(d,jmin:jmax)) ! P = n*kb*t
	
	do i = jmin,jmax
		!new way of finding Bfil with simpler mass conservation formula (analytically the same as other method)
		if (i == dimj) then
			l = 0.5*sqrt((s(x,i)-s(x,i-1))**2+(s(z,i)-s(z,i-1))**2)
			s(b,i) = s(d,i)*Re*l/s(mf,i)
		else if (i == 1) then
			l = 0.5*sqrt((s(x,i)-s(x,i+1))**2+(s(z,i)-s(z,i+1))**2)
			s(b,i) = s(d,i)*Re*l/s(mf,i)
		else
			l = 0.5*sqrt((s(x,i+1)-s(x,i))**2+(s(z,i+1)-s(z,i))**2) &
			& + 0.5*sqrt((s(x,i-1)-s(x,i))**2+(s(z,i-1)-s(z,i))**2)
			s(b,i) = s(d,i)*Re*l/s(mf,i)
		end if
	end do
	
end subroutine algebra1

subroutine findRho1(s,Rho,id) ! 1st attept to update number density using algebraic apporach
	use dimM
	use stateIndex
	use constants
	use errorMod
	implicit none
	real,dimension(dimi,dimj),intent(in) :: s
	real,dimension(dimj) :: Rho
	integer :: i,j,k,id,jmin,jmax
	integer :: point
	real :: coefA,coefB,l,coef,Bx,Bz,pBack

	jmin = thr_domain(id+1,1)
	jmax = thr_domain(id+1,2)
	do i = 1,dimj
		if (i == 1) then
			l = 0.5*sqrt((s(x,i)-s(x,i+1))**2+(s(z,i)-s(z,i+1))**2)
			coef = 2.0*mu0*s(mf,i)*s(mf,i)/(l*l*Re*Re)
			coefA = coef*s(pdg,i)
			coefB = coef*(pBack(s(x,i),s(z,i),i)+((Bx(s(x,i),s(z,i),i))**2+(Bz(s(x,i),s(z,i),i))**2)/(2.0*mu0))
			call newtonRho2(coefA,coefB,s(d,i),Rho(i),i)
		else if (i == dimj) then
			l = 0.5*sqrt((s(x,i)-s(x,i-1))**2+(s(z,i)-s(z,i-1))**2)
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

end subroutine findRho1

subroutine findRho2(s,Rho,id) ! 1st attept to update number density using algebraic apporach
	use dimM
	use stateIndex
	use constants
	use errorMod
	implicit none
	real,dimension(dimi,dimj),intent(in) :: s
	real,dimension(dimj) :: Rho
	integer :: i,j,k,id,jmin,jmax
	integer :: point
	real :: coefA,coefB,l,coef,Bx,Bz,pBack,ptotal
	jmin = thr_domain(id+1,1)
	jmax = thr_domain(id+1,2)

	do i = jmin,jmax
		if (i == 1) then
			l = 0.5*sqrt((s(x,i)-s(x,i+1))**2+(s(z,i)-s(z,i+1))**2)
			coef = 2.0*mu0*s(mf,i)*s(mf,i)/(l*l*Re*Re)
			coefA = coef*s(pdg,i)
			coefB = coef*ptotal(s(x,i),s(z,i),i)
			call newtonRho2(coefA,coefB,s(d,i),Rho(i),i)
		else if (i == dimj) then
			l = 0.5*sqrt((s(x,i)-s(x,i-1))**2+(s(z,i)-s(z,i-1))**2)
			coef = 2.0*mu0*s(mf,i)*s(mf,i)/(l*l*Re*Re)
			coefA = coef*s(pdg,i)
			coefB = coef*ptotal(s(x,i),s(z,i),i)
			call newtonRho2(coefA,coefB,s(d,i),Rho(i),i)
		else
			l = 0.5*sqrt((s(x,i+1)-s(x,i))**2+(s(z,i+1)-s(z,i))**2) &
			& + 0.5*sqrt((s(x,i-1)-s(x,i))**2+(s(z,i-1)-s(z,i))**2)
			coef = 2.0*mu0*s(mf,i)*s(mf,i)/(l*l*Re*Re)
			coefA = coef*s(pdg,i)
			coefB = coef*ptotal(s(x,i),s(z,i),i)
			call newtonRho2(coefA,coefB,s(d,i),Rho(i),i)
		end if
	end do

end subroutine findRho2

subroutine newtonRho2(A,B,ro,rn,n) ! better newton method for density equaiton
	use constants,only: gamma
	use errorMod
	implicit none
	real,intent(in) :: A,B,ro
	real,intent(out) :: rn
	real :: x,xo,f,fp
	integer :: i,n
	integer :: maxTry = 5000
	real :: tolerance
	tolerance = 1.0e-8
	x = ro
	do i = 1,maxTry
		f = x**2 + A*(x)**(gamma) - B
		fp = 2.0*x + gamma*A*(x)**(gamma-1)
		xo = x
		x = x - f/fp
		if ( abs(f/fp) < tolerance ) then
			rn = x
			return
		end if
	end do
	rn = x
	!write(6,*) 'Error, method did not converge!'
	!write(6,*) 'Rho=',x,'Rho0=',ro,'A=',A,'B=',B,'index=',n,'abs(f)=',abs(f)
	!errorFlag = .true.
	!call dumpState
	!stop
end subroutine newtonRho2
