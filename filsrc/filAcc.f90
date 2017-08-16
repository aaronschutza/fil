! subroutines for calculating the acceleration of mass elements

subroutine accel4(s,acc,jmin,jmax) ! accel of mass elements in Re/s^2. Interpolated Gradients on Ptot used for spatial derivatives
	use dimM
	use stateIndex
	use constants
	implicit none
	real,dimension(dimi,dimj) :: s
	real,dimension(vx:vz,dimj) :: acc,a1,a2,a3
	real,dimension(dimj) :: Bfil_x,Bfil_z,Bmed_x,Bmed_z,dxdl,dzdl,ds
	integer :: i,jmin,jmax,rmin,rmax
	real :: l,dx1,dz1,dx2,dz2,dx3,dz3,Bx,Bz,pBack,dxPtot,dzPtot,ptotal,cfactor_1,cfactor_2
	
	if(jmin==1)then
		rmin = jmin
	else
		rmin = jmin-1
	endif

	if(jmax==dimj)then
		rmax = jmax
	else
		rmax = jmax+1
	endif

	do i = rmin,rmax
		Bmed_x(i) = Bx(s(x,i),s(z,i),i)
		Bmed_z(i) = Bz(s(x,i),s(z,i),i)
		if (i==1) then
			! one sided difference for boundaries ******************
			l = sqrt( (s(x,i)-s(x,i+2))**2 + (s(z,i)-s(z,i+2))**2 )
			dx1  = -( s(x,i) - s(x,i+2) )/l
			dz1  = -( s(z,i) - s(z,i+2) )/l
			
			l = sqrt( (s(x,i)-s(x,i+1))**2 + (s(z,i)-s(z,i+1))**2 )
			dx2  = -( s(x,i) - s(x,i+1) )/l
			dz2  = -( s(z,i) - s(z,i+1) )/l
			dx3  = 2.0*dx2 - dx1
			dz3  = 2.0*dz2 - dz1
			
			l = sqrt(dx3*dx3 + dz3*dz3)
			dxdl(i) = dx3/l
			dzdl(i) = dz3/l
		
		elseif (i==dimj) then
			! second boundary **************************************
			l = sqrt( (s(x,i)-s(x,i-2))**2 + (s(z,i)-s(z,i-2))**2 )
			dx1  = ( s(x,i) - s(x,i-2) )/l
			dz1  = ( s(z,i) - s(z,i-2) )/l
					
			l = sqrt( (s(x,i)-s(x,i-1))**2 + (s(z,i)-s(z,i-1))**2 )
			dx2  = ( s(x,i) - s(x,i-1) )/l
			dz2  = ( s(z,i) - s(z,i-1) )/l
			dx3  = 2.0*dx2 - dx1
			dz3  = 2.0*dz2 - dz1
			
			l = sqrt(dx3*dx3 + dz3*dz3)
			dxdl(i) = dx3/l
			dzdl(i) = dz3/l
		else
		! central difference for dX/dl ****************************************
			ds(i) = sqrt((s(x,i+1)-s(x,i))**2+(s(z,i+1)-s(z,i))**2) &
				& + sqrt((s(x,i)-s(x,i-1))**2+(s(z,i)-s(z,i-1))**2)
			dxdl(i) = (s(x,i+1) - s(x,i-1))/ds(i)
			dzdl(i) = (s(z,i+1) - s(z,i-1))/ds(i)
		endif
	end do

	Bfil_x(rmin:rmax) = s(b,rmin:rmax)*dxdl(rmin:rmax)
	Bfil_z(rmin:rmax) = s(b,rmin:rmax)*dzdl(rmin:rmax)
	
	cfactor_1 = mu0*nano*mi*cm*Re*Re
	cfactor_2 = mi*cm*Re*Re*nano
	do i = jmin,jmax
		if(i==1)then
			acc(vx:vz,i) = 0.0
		elseif(i==dimj)then
			acc(vx:vz,i) = 0.0
		else
        	acc(vx,i) = s(b,i)*(Bfil_x(i+1)-Bfil_x(i-1))/(cfactor_1*s(d,i)*ds(i)) &
        			& - dxPtot(s(x,i),s(z,i),i)/(cfactor_2*s(d,i))
        	acc(vz,i) = s(b,i)*(Bfil_z(i+1)-Bfil_z(i-1))/(cfactor_1*s(d,i)*ds(i)) &
        			& - dzPtot(s(x,i),s(z,i),i)/(cfactor_2*s(d,i))
        endif
	end do
	
end subroutine accel4

subroutine deltaDPT(s,ddpt,jmin,jmax) ! the derivatives of d t and p are zero (updated algebraically)
	use dimM
	use stateIndex
	implicit none
	real,dimension(dimi,dimj) :: s
	real,dimension(d:t,dimj) :: ddpt
	integer :: jmin,jmax
	ddpt(:,jmin:jmax) = 0.0
	
end subroutine deltaDPT

subroutine accel_drag(s,acc,jmin,jmax) ! accel of mass elements in Re/s^2. Interpolated Gradients on Ptot used for spatial derivatives
	use dimM
	use stateIndex
	use constants
	use boundMod
	implicit none
	real,dimension(dimi,dimj) :: s
	real,dimension(vx:vz,dimj) :: acc
	integer :: i,jmin,jmax
	real :: ui
	do i=jmin,jmax
		ui = sqrt(s(vx,i)**2+s(vz,i)**2)
		acc(vx,i) = -abs(drag_coeff)*ui*s(vx,i)
		acc(vz,i) = -abs(drag_coeff)*ui*s(vz,i)
	enddo
end subroutine accel_drag

subroutine accel_drive(s,acc,jmin,jmax,time) ! sinusoidal driving force
	use dimM
	use stateIndex
	use constants
	use boundMod
	implicit none
	real,dimension(dimi,dimj) :: s
	real,dimension(vx:vz,dimj) :: acc
	integer :: i,jmin,jmax
	real :: time
	real :: drive_omega
	
	drive_omega = drive_omega_0 + (drive_omega_f-drive_omega_0)*(time-drive_tmin)/(drive_tmax-drive_tmin)
	do i = jmin,jmax
		acc(vx,i) = drive_coeff*sin(drive_omega*(time-drive_tmin))
		acc(vz,i) = 0.0
	enddo
end subroutine accel_drive

