!************************* Boundary Condition Functions (filament endpoints) *****************************

function sigPT(t,x,z,J) ! Time dependent Pedersen conductivity
	use boundMod
	use constants
	implicit none
	real :: sigPT
	real :: t
	real :: x,z
	real :: J
	logical :: jlogic
	
	!if(mod(n,1000) == 0)then
	!	write(*,*) z,J
	!end if
	
	! if(abs(J)>Jthreshold .and. .not. jswitch)then
	! 	jswitch = .true.
	! end if
	
	! if(jSticky)then
	! 	jlogic = jswitch
	! else
	! 	jlogic = abs(J)>Jthreshold
	! end if
	
	! if(Jzjump)then
	! 	jlogic = abs(z) < Jthreshold
	! end if
	
	! if(jlogic)then
	! 	sigPT = sigmaP*condScale
	! else
	! 	sigPT = sigmaP*condScale
	! end if

	sigPT = sigmaP*condScale
	
end function sigPT

function boundVz1(s,time) ! boundary velocity for the first element [Re/s] (conducting wall in z direction)
	use dimM
	use stateIndex
	use constants
	implicit none
	real,dimension(dimi,dimj) :: s
	real :: boundVz1,l,dxdl1,dxdl2,dzdl1,dzdl2,dx3,dz3,Bfil,Bfilx,Bfilz,dBz,dBx,dxdl,dzdl,Bz,time
	real :: sigPT,pulse
	l = sqrt( (s(x,1)-s(x,1+2) )**2 + (s(z,1)-s(z,1+2))**2)
	dxdl1  = -(s(x,1)-s(x,1+2))/l
    dzdl1  = -(s(z,1)-s(z,1+2))/l
    l = sqrt( (s(x,1)-s(x,1+1) )**2 + (s(z,1)-s(z,1+1))**2)
    dxdl2  = -(s(x,1)-s(x,1+1))/l
    dzdl2  = -(s(z,1)-s(z,1+1))/l
    dx3  = 2.0*dxdl2 - dxdl1
    dz3  = 2.0*dzdl2 - dzdl1
    l = sqrt(dx3*dx3+dz3*dz3)
    dxdl  = dx3/l
    dzdl  = dz3/l
    Bfil   = s(b,1)
    Bfilx  = dxdl*Bfil
    Bfilz  = dzdl*Bfil
    !dBx = Bfilx - Bx(s(x,dimj),s(z,dimj))
    dBz = Bfilz - Bz(s(x,1),s(z,1),1)
    l = 1.0 !(1.0+39.0*time/3600.0) 
	boundVz1 = pulse(time) - cm*dBz/(Re*mu0*sigPT(time,s(x,1),s(z,1),dBz)*l*Bfilx)
end function boundVz1

function boundVz2(s,time) ! boundary velocity for the last element [Re/s] (conducting wall in z direction)
	use dimM
	use stateIndex
	use constants
	implicit none
	real,dimension(dimi,dimj) :: s
	real :: boundVz2,l,dxdl1,dxdl2,dzdl1,dzdl2,dx3,dz3,Bfil,Bfilx,Bfilz,dBz,dBx,dxdl,dzdl,Bz,time
	real :: sigPT,pulse
	l = sqrt( (s(x,dimj)-s(x,dimj-2) )**2 + (s(z,dimj)-s(z,dimj-2))**2)
	dxdl1  = (s(x,dimj)-s(x,dimj-2))/l
    dzdl1  = (s(z,dimj)-s(z,dimj-2))/l
    l = sqrt( (s(x,dimj)-s(x,dimj-1) )**2 + (s(z,dimj)-s(z,dimj-1))**2)
    dxdl2  = (s(x,dimj)-s(x,dimj-1))/l
    dzdl2  = (s(z,dimj)-s(z,dimj-1))/l
    dx3  = 2.0*dxdl2 - dxdl1
    dz3  = 2.0*dzdl2 - dzdl1
    l = sqrt(dx3*dx3+dz3*dz3)
    dxdl  = dx3/l
    dzdl  = dz3/l
    Bfil   = s(b,dimj)
    Bfilx  = dxdl*Bfil
    Bfilz  = dzdl*Bfil
    !dBx = Bfilx - Bx(s(x,dimj),s(z,dimj))
    dBz = Bfilz - Bz(s(x,dimj),s(z,dimj),dimj)
    l = 1.0 !(1.0+39.0*time/3600.0)
	boundVz2 = -pulse(time) - cm*dBz/(Re*mu0*sigPT(time,s(x,dimj),s(z,dimj),dBz)*l*Bfilx)
end function boundVz2

function roundVz1(s,time) ! boundary velocity for the first element [Re/s] (conducting shell)
	use dimM
	use stateIndex
	use constants
	implicit none
	real,dimension(dimi,dimj) :: s
	real :: roundVz1,l,dxdl1,dxdl2,dzdl1,dzdl2,dx3,dz3,Bfil,Bfilx,Bfilz,dBz,dBx,dxdl,dzdl,Bz,Bx,time
	real :: sigPT,pulse
	real :: BfilR,dBtheta
	l = sqrt( (s(x,1)-s(x,1+2) )**2 + (s(z,1)-s(z,1+2))**2)
	dxdl1  = -(s(x,1)-s(x,1+2))/l
    dzdl1  = -(s(z,1)-s(z,1+2))/l
    l = sqrt( (s(x,1)-s(x,1+1) )**2 + (s(z,1)-s(z,1+1))**2)
    dxdl2  = -(s(x,1)-s(x,1+1))/l
    dzdl2  = -(s(z,1)-s(z,1+1))/l
    dx3  = 2.0*dxdl2 - dxdl1
    dz3  = 2.0*dzdl2 - dzdl1
    l = sqrt(dx3*dx3+dz3*dz3)
    dxdl  = dx3/l
    dzdl  = dz3/l
    Bfil   = s(b,1)
    Bfilx  = dxdl*Bfil
    Bfilz  = dzdl*Bfil
    dBx = Bfilx - Bx(s(x,1),s(z,1),1)
    dBz = Bfilz - Bz(s(x,1),s(z,1),1)
    l = sqrt(s(x,1)**2+s(z,1)**2)
    BfilR = (Bfilx*s(x,1)+Bfilz*s(z,1))/l
    dBtheta = (-dBx*s(z,1)+dBz*s(x,1))/l
	roundVz1 = cm*dBtheta*s(x,1)/(Re*mu0*sigPT(time,s(x,1),s(z,1),dBz)*BfilR*l)
end function roundVz1

function roundVz2(s,time) ! boundary velocity for the last element [Re/s] (conducting shell)
	use dimM
	use stateIndex
	use constants
	implicit none
	real,dimension(dimi,dimj) :: s
	real :: roundVz2,l,dxdl1,dxdl2,dzdl1,dzdl2,dx3,dz3,Bfil,Bfilx,Bfilz,dBz,dBx,dxdl,dzdl,Bz,Bx,time
	real :: sigPT,pulse
	real :: BfilR,dBtheta
	l = sqrt( (s(x,dimj)-s(x,dimj-2) )**2 + (s(z,dimj)-s(z,dimj-2))**2)
	dxdl1  = (s(x,dimj)-s(x,dimj-2))/l
    dzdl1  = (s(z,dimj)-s(z,dimj-2))/l
    l = sqrt( (s(x,dimj)-s(x,dimj-1) )**2 + (s(z,dimj)-s(z,dimj-1))**2)
    dxdl2  = (s(x,dimj)-s(x,dimj-1))/l
    dzdl2  = (s(z,dimj)-s(z,dimj-1))/l
    dx3  = 2.0*dxdl2 - dxdl1
    dz3  = 2.0*dzdl2 - dzdl1
    l = sqrt(dx3*dx3+dz3*dz3)
    dxdl  = dx3/l
    dzdl  = dz3/l
    Bfil   = s(b,dimj)
    Bfilx  = dxdl*Bfil
    Bfilz  = dzdl*Bfil
    dBx = Bfilx - Bx(s(x,dimj),s(z,dimj),dimj)
    dBz = Bfilz - Bz(s(x,dimj),s(z,dimj),dimj)
    l = sqrt(s(x,dimj)**2+s(z,dimj)**2)
    BfilR = (Bfilx*s(x,dimj)+Bfilz*s(z,dimj))/l
    dBtheta = (-dBx*s(z,dimj)+dBz*s(x,dimj))/l
	roundVz2 = cm*dBtheta*s(x,dimj)/(Re*mu0*sigPT(time,s(x,1),s(z,dimj),dBz)*BfilR*l)
end function roundVz2

function roundVx1(s,time) ! boundary velocity for the first element [Re/s] (conducting shell)
	use dimM
	use stateIndex
	use constants
	implicit none
	real,dimension(dimi,dimj) :: s
	real :: roundVx1,l,dxdl1,dxdl2,dzdl1,dzdl2,dx3,dz3,Bfil,Bfilx,Bfilz,dBz,dBx,dxdl,dzdl,Bz,Bx,time
	real :: sigPT,pulse
	real :: BfilR,dBtheta
	l = sqrt( (s(x,1)-s(x,1+2) )**2 + (s(z,1)-s(z,1+2))**2)
	dxdl1  = -(s(x,1)-s(x,1+2))/l
    dzdl1  = -(s(z,1)-s(z,1+2))/l
    l = sqrt( (s(x,1)-s(x,1+1) )**2 + (s(z,1)-s(z,1+1))**2)
    dxdl2  = -(s(x,1)-s(x,1+1))/l
    dzdl2  = -(s(z,1)-s(z,1+1))/l
    dx3  = 2.0*dxdl2 - dxdl1
    dz3  = 2.0*dzdl2 - dzdl1
    l = sqrt(dx3*dx3+dz3*dz3)
    dxdl  = dx3/l
    dzdl  = dz3/l
    Bfil   = s(b,1)
    Bfilx  = dxdl*Bfil
    Bfilz  = dzdl*Bfil
    dBx = Bfilx - Bx(s(x,1),s(z,1),1)
    dBz = Bfilz - Bz(s(x,1),s(z,1),1)
    l = sqrt(s(x,1)**2+s(z,1)**2)
    BfilR = (Bfilx*s(x,1)+Bfilz*s(z,1))/l
    dBtheta = (-dBx*s(z,1)+dBz*s(x,1))/l
	roundVx1 = -cm*dBtheta*s(z,1)/(Re*mu0*sigPT(time,s(x,1),s(z,1),dBz)*BfilR*l)
end function roundVx1

function roundVx2(s,time) ! boundary velocity for the last element [Re/s] (conducting shell)
	use dimM
	use stateIndex
	use constants
	implicit none
	real,dimension(dimi,dimj) :: s
	real :: roundVx2,l,dxdl1,dxdl2,dzdl1,dzdl2,dx3,dz3,Bfil,Bfilx,Bfilz,dBz,dBx,dxdl,dzdl,Bz,Bx,time
	real :: sigPT,pulse
	real :: BfilR,dBtheta
	l = sqrt( (s(x,dimj)-s(x,dimj-2) )**2 + (s(z,dimj)-s(z,dimj-2))**2)
	dxdl1  = (s(x,dimj)-s(x,dimj-2))/l
    dzdl1  = (s(z,dimj)-s(z,dimj-2))/l
    l = sqrt( (s(x,dimj)-s(x,dimj-1) )**2 + (s(z,dimj)-s(z,dimj-1))**2)
    dxdl2  = (s(x,dimj)-s(x,dimj-1))/l
    dzdl2  = (s(z,dimj)-s(z,dimj-1))/l
    dx3  = 2.0*dxdl2 - dxdl1
    dz3  = 2.0*dzdl2 - dzdl1
    l = sqrt(dx3*dx3+dz3*dz3)
    dxdl  = dx3/l
    dzdl  = dz3/l
    Bfil   = s(b,dimj)
    Bfilx  = dxdl*Bfil
    Bfilz  = dzdl*Bfil
    dBx = Bfilx - Bx(s(x,dimj),s(z,dimj),dimj)
    dBz = Bfilz - Bz(s(x,dimj),s(z,dimj),dimj)
    l = sqrt(s(x,dimj)**2+s(z,dimj)**2)
    BfilR = (Bfilx*s(x,dimj)+Bfilz*s(z,dimj))/l
    dBtheta = (-dBx*s(z,dimj)+dBz*s(x,dimj))/l
	roundVx2 = -cm*dBtheta*s(z,dimj)/(Re*mu0*sigPT(time,s(x,1),s(z,dimj),dBz)*BfilR*l)
end function roundVx2

function boundVz1_test(time) ! boundary velocity for the first element [Re/s] (wall with parameterized velocity)
	use constants
	implicit none
	real :: boundVz1_test,time
	if(time>js)then
		boundVz1_test = 0.0
	else
		boundVz1_test = -jolt*sin(time*2.0*Pi/jt)
	end if
end function boundVz1_test

function boundVz2_test(time) ! boundary velocity for the last element [Re/s] (wall with parameterized velocity)
	use constants
	implicit none
	real :: boundVz2_test,time
	if(time>js)then
		boundVz2_test = 0.0
	else
		boundVz2_test = jolt*sin(time*2.0*Pi/jt)
	end if
end function boundVz2_test

