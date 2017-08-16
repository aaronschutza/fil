!************** Linear displacement solver routines ***********************
! Used to analyze restoring forces and model motion in terms of low
! amplitude displacements
!**************************************************************************

subroutine littleSolver2d(xg,zg,accxg,acczg,dimij,cx0,cz0)
	use dnams
	implicit none
	integer :: dimij
	real,dimension(dimij) :: xg,zg
	real,dimension(dimij,dimij) :: accxg,acczg
	integer :: i,j,k,kmax
	real :: x0,z0,x,z,coef1,ax,az,vx,vz,cx0,cz0
	real :: t,tau,tmax
	
	t = 0.0
	tau = 0.06
	tmax = 120.0
	kmax = ceiling(tmax/tau)
	
	x0 = 0.0
	z0 = 0.5
	
	x = x0
	z = z0
	
	vx = 0.0
	vz = 0.0
	
	open(unit=11, file=trim(dataPath)//".ls", status="replace")
	write(11,*) cx0,cz0
	write(11,*) t,x,z,vx,vz
	do k = 1,kmax
		
		call littleBilinear(xg,zg,accxg,acczg,dimij,x,z,ax,az)
		
		vx = vx + tau*ax
		vz = vz + tau*az
		x = x + tau*vx
		z = z + tau*vz
		t = t+tau
		
		write(11,*) t,x,z,vx,vz
	enddo
	
	close(11)
	
end subroutine

subroutine littleBilinear(xg,zg,accxg,acczg,dimij,x,z,ax,az)
	implicit none
	integer :: dimij
	real,dimension(dimij) :: xg,zg
	real,dimension(dimij,dimij) :: accxg,acczg
	integer :: i,j
	real :: x,z,coef1,ax,az
	
	do i = 1,dimij-1
		if (xg(i)<=x .and. x<=xg(i+1))then
			exit
		end if
	end do
	
	do j = 1,dimij-1
		if (zg(j)<=z .and. z<=zg(j+1))then
			exit
		end if
	end do
	
	coef1 = (xg(i+1)-xg(i))*(zg(j+1)-zg(j))
	ax = accxg(i,j)*(xg(i+1)-x)*(zg(j+1)-z)/coef1 &
	& +  accxg(i+1,j)*(x-xg(i))*(zg(j+1)-z)/coef1 &
	& +  accxg(i,j+1)*(xg(i+1)-x)*(z-zg(j))/coef1 &
	& +  accxg(i+1,j+1)*(x-xg(i))*(z-zg(j))/coef1
	az = acczg(i,j)*(xg(i+1)-x)*(zg(j+1)-z)/coef1 &
	& +  acczg(i+1,j)*(x-xg(i))*(zg(j+1)-z)/coef1 &
	& +  acczg(i,j+1)*(xg(i+1)-x)*(z-zg(j))/coef1 &
	& +  acczg(i+1,j+1)*(x-xg(i))*(z-zg(j))/coef1

end subroutine