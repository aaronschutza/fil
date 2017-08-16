!********************* Special functions and routines *****************************
! functions for calculating specific quantities from filament configuration
! auxiliary functions for running special tests
!**********************************************************************************

function Kbg1(x)
! x_e dependent entropy function [nPa*(Re/nT)^gamma]
	use constants
	implicit none
	real :: Kbg1
	real :: x,pBack,Vbg1
	Kbg1 = pBack(x,0.0,-1)*(Vbg1(x))**(gamma)
end function Kbg1

function Kbg2(z)
! z_bound dependent entropy function [nPa*(Re/nT)^gamma]
	use constants
	use filConstants
	implicit none
	real :: Kbg2
	real :: z,pBack,Vbg2
	Kbg2 = pBack(wall,z,-1)*(Vbg2(z))**(gamma)
end function Kbg2

function filVol(s)
! filament volume [Re/nT]
	use dimM
	use stateIndex
	use constants
	use filConstants
	implicit none
	real :: filVol
	real,dimension(dimi,dimj) :: s
	real :: l
	integer :: i
	filVol = 0.0
	do i = 1,dimj
		if (i == dimj) then
			l = 0.5*sqrt( (s(x,i)-s(x,i-1))**2 + (s(z,i)-s(z,i-1))**2 )
			filVol = filVol + l/s(b,i)
		else if (i == 1) then
			l = 0.5*sqrt( (s(x,i)-s(x,i+1))**2 + (s(z,i)-s(z,i+1))**2 )
			filVol = filVol + l/s(b,i)
		else
			l = 0.5*sqrt( (s(x,i+1)-s(x,i))**2 + (s(z,i+1)-s(z,i))**2 ) &
			& + 0.5*sqrt( (s(x,i-1)-s(x,i))**2 + (s(z,i-1)-s(z,i))**2 )
			filVol = filVol + l/s(b,i)
		end if
	end do
end function filVol

real function filBeta(s) ! field line averaged beta
	use dimM
	use stateIndex
	use constants
	use filConstants
	implicit none
	real :: filVol
	real,dimension(dimi,dimj) :: s
	real :: l
	integer :: i
	real :: c0,pBack
	c0 = 0.0
	do i = 1,dimj
		if (i == dimj) then
			l = 0.5*sqrt( (s(x,i)-s(x,i-1))**2 + (s(z,i)-s(z,i-1))**2 )
			c0 = c0 + s(p,i)*l/s(b,i)**3
		else if (i == 1) then
			l = 0.5*sqrt( (s(x,i)-s(x,i+1))**2 + (s(z,i)-s(z,i+1))**2 )
			c0 = c0 + s(p,i)*l/s(b,i)**3
		else
			l = 0.5*sqrt( (s(x,i+1)-s(x,i))**2 + (s(z,i+1)-s(z,i))**2 ) &
			& + 0.5*sqrt( (s(x,i-1)-s(x,i))**2 + (s(z,i-1)-s(z,i))**2 )
			c0 = c0 + s(p,i)*l/s(b,i)**3
		end if
	end do
	filBeta = 2*mu0*c0/filVol(s)
end function

function Vbg1(x)
! x_e dependent background volume [Re/nT]
	use constants
	use parametersBackground
	use filConstants
	implicit none
	real :: Vbg1
	real :: x
	real :: u,c0,c1,c2
	u = alpha*(x-wall)
	c0 = 4*ht*exp(u)
	c1 = acos(exp(-u))
	c2 = pi*alpha*A0*exp(-wall*alpha)
	Vbg1 = c0*c1/c2
end function Vbg1

function Vbg2(z)
! z_bound dependent background volume [Re/nT]
	use constants
	use parametersBackground
	use filConstants
	implicit none
	real :: Vbg2
	real :: z
	real :: zeta,c0,c1
	zeta = 0.5*pi*z/ht
	c0 = 4*ht*zeta
	c1 = pi*alpha*A0*exp(-wall*alpha)*cos(zeta)
	Vbg2 = c0/c1
end function Vbg2

subroutine findMf2(s,mFil) 
! find the mass of each element by integrating dm = Rho/B ds to adjacent elements [imposed to be constant] ( 1/(nT*cm^2) ) inverse flux
	use dimM
	use stateIndex
	use constants
	use filConstants
	implicit none
	real,dimension(dimi,dimj) :: s
	real,dimension(dimj) :: mFil
	real :: l
	integer :: i
	do i = 1,dimj
		if (i == dimj) then
			l = 0.5*sqrt( (s(x,i)-s(x,i-1))**2 + (s(z,i)-s(z,i-1))**2 )
			mFil(i) = Re*l*s(d,i)/s(b,i)
		else if (i == 1) then
			l = 0.5*sqrt( (s(x,i)-s(x,i+1))**2 + (s(z,i)-s(z,i+1))**2 )
			mFil(i) = Re*l*s(d,i)/s(b,i)
		else
			l = 0.5*sqrt( (s(x,i+1)-s(x,i))**2 + (s(z,i+1)-s(z,i))**2 ) &
			& + 0.5*sqrt( (s(x,i-1)-s(x,i))**2 + (s(z,i-1)-s(z,i))**2 )
			mFil(i) = Re*l*s(d,i)/s(b,i)
		end if
	end do
	totalM = 0.0
	do i = 1,dimj
		totalM = totalM + mFil(i)
	end do
end subroutine findMf2

real function massFil(s) 
! find the mass by integrating dm = Rho/B ds to adjacent elements [imposed to be constant] ( 1/(nT*cm^2) ) inverse flux
	use dimM
	use stateIndex
	use constants
	use filConstants
	implicit none
	real,dimension(dimi,dimj) :: s
	real :: l,c0
	integer :: i
	c0 = 0.0
	do i = 1,dimj
		if (i == dimj) then
			l = 0.5*sqrt( (s(x,i)-s(x,i-1))**2 + (s(z,i)-s(z,i-1))**2 )
			c0 = c0 + Re*l*s(d,i)/s(b,i)
		else if (i == 1) then
			l = 0.5*sqrt( (s(x,i)-s(x,i+1))**2 + (s(z,i)-s(z,i+1))**2 )
			c0 = c0 + Re*l*s(d,i)/s(b,i)
		else
			l = 0.5*sqrt( (s(x,i+1)-s(x,i))**2 + (s(z,i+1)-s(z,i))**2 ) &
			& + 0.5*sqrt( (s(x,i-1)-s(x,i))**2 + (s(z,i-1)-s(z,i))**2 )
			c0 = c0 + Re*l*s(d,i)/s(b,i)
		end if
	end do
	
	massFil = c0
	
end function

function filK(s)
	use dimM
	use stateIndex
	use constants
	use filConstants
	implicit none
	real :: filK
	real,dimension(dimi,dimj) :: s
	real :: l
	integer :: i
	filK = 0.0
	do i = 1,dimj
		if (i == dimj) then
			l = 0.5*sqrt((s(x,i)-s(x,i-1))**2+(s(z,i)-s(z,i-1))**2)
			filK = filK + (l/s(b,i))*(s(p,i)**(1/gamma))
		else if (i == 1) then
			l = 0.5*sqrt((s(x,i)-s(x,i+1))**2+(s(z,i)-s(z,i+1))**2)
			filK = filK + (l/s(b,i))*(s(p,i)**(1/gamma))
		else
			l = 0.5*sqrt((s(x,i+1)-s(x,i))**2+(s(z,i+1)-s(z,i))**2) &
			& + 0.5*sqrt((s(x,i-1)-s(x,i))**2+(s(z,i-1)-s(z,i))**2)
			filK = filK + (l/s(b,i))*(s(p,i)**(1/gamma))
		end if
	end do
	
	filK = filK**gamma
	
end function filK

function fil_KinNRG(s)
	!kinetic energy per unit magnetic flux [J/Weber]
	use dimM
	use stateIndex
	use constants
	use filConstants
	implicit none
	real :: fil_KinNRG
	real,dimension(dimi,dimj) :: s
	real :: l
	integer :: i
	fil_KinNRG = 0.0

	do i = 1,dimj
		if (i == 1) then
			l = sqrt( (s(x,i+1)-s(x,i))**2 + (s(z,i+1)-s(z,i))**2)
		elseif (i == dimj)then
			l = sqrt( (s(x,i)-s(x,i-1))**2 + (s(z,i)-s(z,i-1))**2)
		else
			l = 0.5*sqrt( (s(x,i+1)-s(x,i-1))**2 + (s(z,i+1)-s(z,i-1))**2)
		endif
		fil_KinNRG = fil_KinNRG &
		& + 0.5*s_0(d,i)*l*sqrt((s(vx,i))**2+(s(vz,i))**2)/s(b,i)
	enddo

	fil_KinNRG = fil_KinNRG*mi*(Re**3)*nano

end function fil_KinNRG

subroutine findH(x,z,h) ! adaptive step size for background derivatives
! set to constant for now (works with position dependence)
	implicit none
	real :: x,z,h
	h = 0.04
end subroutine

function pulse(t)
	use constants
	implicit none
	real :: t
	real :: pulse
	real :: dt
	dt = 60
	!if (t<dt)then
!		
!		pulse = 0.01*sin(pi*t/dt)
!	else
!		pulse = 0.0
!	endif
	pulse = 0.0
end function pulse

function alphaT(t)
	use parametersBackground
	implicit none
	real :: alphaT
	real :: t
	
	alphaT = alpha0/(1-(0.001/13.0)*t)
	
end function alphaT

