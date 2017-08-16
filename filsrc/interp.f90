!*********************** Background interpolation routines ******************************

real function interpolate(mesh,x,z,elnum)
	use parametersBackground, only: interpType
	implicit none
	real,intent(in) :: x
	real,intent(in) :: z
    integer :: elnum
    character(len=*) :: mesh
    real :: bilinearInterp,bicubicInterp
    real :: dipx,dipy,dipz
    select case(interpType)
    	case(0)
    		interpolate = bilinearInterp(mesh,x,z,elnum)
    	case(1)
    		interpolate = bicubicInterp(mesh,x,z,elnum)
    	case default
    		write(*,*) 'invalid case in interpolate'
    		stop
    end select
    		! if (useAnalyticDipole)then
			! 	!SUBROUTINE Dipole_modified (dipole_moment, tilt_angle, x, y, z, bx, by, bz)
			! 	call Dipole_modified(dipole_moment, 0.0, x, 0.0, z, dipx, dipy, dipz)
			! 	bilinearInterp = bilinearInterp + dipx
			! endif
			! if (useAnalyticDipole) then
			! 	!SUBROUTINE Dipole_modified (dipole_moment, tilt_angle, x, y, z, bx, by, bz)
			! 	call Dipole_modified(dipole_moment, 0.0, x, 0.0, z, dipx, dipy, dipz)
			! 	bilinearInterp = bilinearInterp + dipz
			! endif

end function

real function bilinearInterp(mesh,x,z,elnum)
	use potential
	use dimM
	implicit none
	real,intent(in) :: x
	real,intent(in) :: z
	real :: term1,coef1
    integer :: i,j,elnum
    character(len=*) :: mesh
	
	if(x>xmaxl .or. z>zmaxl .or. x<xminl .or. z<zminl)then
		write(*,*) 'outofbounds',elnum
		write(*,*) xminl,x,xmaxl
		write(*,*) zminl,z,zmaxl
		stop
	else
!********************************************************************************
		if (elnum < 1 .or. elnum > dimj) then
			do i = 1,gridi-1
				if (xg(i)<=x .and. x<=xg(i+1))then
					exit
				end if
			end do
			do j = 1,gridj-1
				if (zg(j)<=z .and. z<=zg(j+1))then
					exit
				end if
			end do
		elseif (recallPos(elnum,1) == -1 .or. recallPos(elnum,2) == -1) then
			do i = 1,gridi-1
				if (xg(i)<=x .and. x<=xg(i+1))then
					exit
				end if
			end do
			do j = 1,gridj-1
				if (zg(j)<=z .and. z<=zg(j+1))then
					exit
				end if
			end do
			recallPos(elnum,1) = i
			recallPos(elnum,2) = j
		else
			if(xg(recallPos(elnum,1))<=x .and. x<=xg(recallPos(elnum,1)+1))then
				i = recallPos(elnum,1)
			elseif(recallPos(elnum,1)<gridi-1 .and. recallPos(elnum,1) > 1)then
				if(xg(recallPos(elnum,1)+1)<=x .and. x<=xg(recallPos(elnum,1)+2))then
					recallPos(elnum,1) = recallPos(elnum,1)+1
					i = recallPos(elnum,1)
				elseif(xg(recallPos(elnum,1)-1)<=x .and. x<=xg(recallPos(elnum,1)))then
					recallPos(elnum,1) = recallPos(elnum,1)-1
					i = recallPos(elnum,1)
				else
					do i = 1,gridi-1
						if (xg(i)<=x .and. x<=xg(i+1))then
							exit
						endif
					end do
					recallPos(elnum,1) = i
				endif
			elseif(recallPos(elnum,1)==gridi-1)then
				if(xg(recallPos(elnum,1)-1)<=x .and. x<=xg(recallPos(elnum,1)))then
					recallPos(elnum,1) = recallPos(elnum,1)-1
					i = recallPos(elnum,1)
				else
					do i = 1,gridi-1
						if (xg(i)<=x .and. x<=xg(i+1))then
							exit
						end if
					end do
					recallPos(elnum,1) = i
				endif
			elseif(recallPos(elnum,1)==1)then
				if(xg(recallPos(elnum,1)+1)<=x .and. x<=xg(recallPos(elnum,1)+2))then
					recallPos(elnum,1) = recallPos(elnum,1)+1
					i = recallPos(elnum,1)
				else
					do i = 1,gridi-1
						if (xg(i)<=x .and. x<=xg(i+1))then
							exit
						endif
					end do
					recallPos(elnum,1) = i
				endif
			else
				do i = 1,gridi-1
					if (xg(i)<=x .and. x<=xg(i+1))then
						exit
					end if
				end do
			endif
!***************************************************************************************
			if(zg(recallPos(elnum,2))<=z .and. z<=zg(recallPos(elnum,2)+1))then
				j = recallPos(elnum,2)
			elseif(recallPos(elnum,2)<gridj-1 .and. recallPos(elnum,2) > 1)then
				if(zg(recallPos(elnum,2)+1)<=z .and. z<=zg(recallPos(elnum,2)+2))then
					recallPos(elnum,2) = recallPos(elnum,2)+1
					j = recallPos(elnum,2)
				elseif(zg(recallPos(elnum,2)-1)<=z .and. z<=zg(recallPos(elnum,2)))then
					recallPos(elnum,2) = recallPos(elnum,2)-1
					j = recallPos(elnum,2)
				else
					do j = 1,gridj-1
						if (zg(j)<=z .and. z<=zg(j+1))then
							exit
						endif
					end do
					recallPos(elnum,2) = j
				endif
			elseif(recallPos(elnum,2)==gridj-1)then
				if(zg(recallPos(elnum,2)-1)<=z .and. z<=zg(recallPos(elnum,2)))then
					recallPos(elnum,2) = recallPos(elnum,2)-1
					j = recallPos(elnum,2)
				else
					do j = 1,gridj-1
						if (zg(j)<=z .and. z<=zg(j+1))then
							exit
						end if
					end do
					recallPos(elnum,2) = j
				endif
			elseif(recallPos(elnum,2)==1)then
				if(zg(recallPos(elnum,2)+1)<=z .and. z<=zg(recallPos(elnum,2)+2))then
					recallPos(elnum,2) = recallPos(elnum,2)+1
					j = recallPos(elnum,2)
				else
					do j = 1,gridj-1
						if (zg(j)<=z .and. z<=zg(j+1))then
							exit
						endif
					end do
					recallPos(elnum,2) = j
				endif
			else
				do j = 1,gridj-1
					if (zg(j)<=z .and. z<=zg(j+1))then
						exit
					end if
				end do
			endif
		endif
	endif
	
	select case(mesh)
		case('Ag')
			coef1 = (xg(i+1)-xg(i))*(zg(j+1)-zg(j))
			bilinearInterp = Ag(i,j)*(xg(i+1)-x)*(zg(j+1)-z)/coef1 &
						& +  Ag(i+1,j)*(x-xg(i))*(zg(j+1)-z)/coef1 &
						& +  Ag(i,j+1)*(xg(i+1)-x)*(z-zg(j))/coef1 &
						& +  Ag(i+1,j+1)*(x-xg(i))*(z-zg(j))/coef1
		case('Bxg')
			coef1 = (xg(i+1)-xg(i))*(zg(j+1)-zg(j))
			bilinearInterp = Bxg(i,j)*(xg(i+1)-x)*(zg(j+1)-z)/coef1 &
						& +  Bxg(i+1,j)*(x-xg(i))*(zg(j+1)-z)/coef1 &
						& +  Bxg(i,j+1)*(xg(i+1)-x)*(z-zg(j))/coef1 &
						& +  Bxg(i+1,j+1)*(x-xg(i))*(z-zg(j))/coef1
		case('Bzg')
			coef1 = (xg(i+1)-xg(i))*(zg(j+1)-zg(j))
			bilinearInterp = Bzg(i,j)*(xg(i+1)-x)*(zg(j+1)-z)/coef1 &
						& +  Bzg(i+1,j)*(x-xg(i))*(zg(j+1)-z)/coef1 &
						& +  Bzg(i,j+1)*(xg(i+1)-x)*(z-zg(j))/coef1 &
						& +  Bzg(i+1,j+1)*(x-xg(i))*(z-zg(j))/coef1
		case('Pg')
			coef1 = (xg(i+1)-xg(i))*(zg(j+1)-zg(j))
			bilinearInterp = Pg(i,j)*(xg(i+1)-x)*(zg(j+1)-z)/coef1 &
						& +  Pg(i+1,j)*(x-xg(i))*(zg(j+1)-z)/coef1 &
						& +  Pg(i,j+1)*(xg(i+1)-x)*(z-zg(j))/coef1 &
						& +  Pg(i+1,j+1)*(x-xg(i))*(z-zg(j))/coef1
		case('Ng')
			coef1 = (xg(i+1)-xg(i))*(zg(j+1)-zg(j))
			bilinearInterp = Ng(i,j)*(xg(i+1)-x)*(zg(j+1)-z)/coef1 &
						& +  Ng(i+1,j)*(x-xg(i))*(zg(j+1)-z)/coef1 &
						& +  Ng(i,j+1)*(xg(i+1)-x)*(z-zg(j))/coef1 &
						& +  Ng(i+1,j+1)*(x-xg(i))*(z-zg(j))/coef1
		case('Ptotg')
			coef1 = (xg(i+1)-xg(i))*(zg(j+1)-zg(j))
			bilinearInterp = Ptotg(i,j)*(xg(i+1)-x)*(zg(j+1)-z)/coef1 &
						& +  Ptotg(i+1,j)*(x-xg(i))*(zg(j+1)-z)/coef1 &
						& +  Ptotg(i,j+1)*(xg(i+1)-x)*(z-zg(j))/coef1 &
						& +  Ptotg(i+1,j+1)*(x-xg(i))*(z-zg(j))/coef1
		case('dx_Bxg')
			coef1 = (xg(i+1)-xg(i))*(zg(j+1)-zg(j))
			bilinearInterp = dx_Bxg(i,j)*(xg(i+1)-x)*(zg(j+1)-z)/coef1 &
						& +  dx_Bxg(i+1,j)*(x-xg(i))*(zg(j+1)-z)/coef1 &
						& +  dx_Bxg(i,j+1)*(xg(i+1)-x)*(z-zg(j))/coef1 &
						& +  dx_Bxg(i+1,j+1)*(x-xg(i))*(z-zg(j))/coef1
		case('dz_Bxg')
			coef1 = (xg(i+1)-xg(i))*(zg(j+1)-zg(j))
			bilinearInterp = dz_Bxg(i,j)*(xg(i+1)-x)*(zg(j+1)-z)/coef1 &
						& +  dz_Bxg(i+1,j)*(x-xg(i))*(zg(j+1)-z)/coef1 &
						& +  dz_Bxg(i,j+1)*(xg(i+1)-x)*(z-zg(j))/coef1 &
						& +  dz_Bxg(i+1,j+1)*(x-xg(i))*(z-zg(j))/coef1
		case('dx_Bzg')
			coef1 = (xg(i+1)-xg(i))*(zg(j+1)-zg(j))
			bilinearInterp = dx_Bzg(i,j)*(xg(i+1)-x)*(zg(j+1)-z)/coef1 &
						& +  dx_Bzg(i+1,j)*(x-xg(i))*(zg(j+1)-z)/coef1 &
						& +  dx_Bzg(i,j+1)*(xg(i+1)-x)*(z-zg(j))/coef1 &
						& +  dx_Bzg(i+1,j+1)*(x-xg(i))*(z-zg(j))/coef1
		case('dz_Bzg')
			coef1 = (xg(i+1)-xg(i))*(zg(j+1)-zg(j))
			bilinearInterp = dz_Bzg(i,j)*(xg(i+1)-x)*(zg(j+1)-z)/coef1 &
						& +  dz_Bzg(i+1,j)*(x-xg(i))*(zg(j+1)-z)/coef1 &
						& +  dz_Bzg(i,j+1)*(xg(i+1)-x)*(z-zg(j))/coef1 &
						& +  dz_Bzg(i+1,j+1)*(x-xg(i))*(z-zg(j))/coef1
		case('dx_Ptotg')
			coef1 = (xg(i+1)-xg(i))*(zg(j+1)-zg(j))
			bilinearInterp = dx_Ptotg(i,j)*(xg(i+1)-x)*(zg(j+1)-z)/coef1 &
						& +  dx_Ptotg(i+1,j)*(x-xg(i))*(zg(j+1)-z)/coef1 &
						& +  dx_Ptotg(i,j+1)*(xg(i+1)-x)*(z-zg(j))/coef1 &
						& +  dx_Ptotg(i+1,j+1)*(x-xg(i))*(z-zg(j))/coef1
		case('dz_Ptotg')
			coef1 = (xg(i+1)-xg(i))*(zg(j+1)-zg(j))
			bilinearInterp = dz_Ptotg(i,j)*(xg(i+1)-x)*(zg(j+1)-z)/coef1 &
						& +  dz_Ptotg(i+1,j)*(x-xg(i))*(zg(j+1)-z)/coef1 &
						& +  dz_Ptotg(i,j+1)*(xg(i+1)-x)*(z-zg(j))/coef1 &
						& +  dz_Ptotg(i+1,j+1)*(x-xg(i))*(z-zg(j))/coef1
		case('testMesh')
			coef1 = (xg(i+1)-xg(i))*(zg(j+1)-zg(j))
			bilinearInterp = testMesh(i,j)*(xg(i+1)-x)*(zg(j+1)-z)/coef1 &
						& +  testMesh(i+1,j)*(x-xg(i))*(zg(j+1)-z)/coef1 &
						& +  testMesh(i,j+1)*(xg(i+1)-x)*(z-zg(j))/coef1 &
						& +  testMesh(i+1,j+1)*(x-xg(i))*(z-zg(j))/coef1
		case default
			write(*,*) 'error: bilinearInterp case not found'
			stop
	end select
end function

real function bicubicInterp(mesh,x,z,elnum)
	use potential
	use dimM
	implicit none
	real,intent(in) :: x
	real,intent(in) :: z
	real :: term1,coef1
    integer :: i,j,elnum
    character(len=*) :: mesh
	real,dimension(4) :: gc,gc_x,gc_z,gc_xz
	real :: ansy,ansy_x,ansy_z
	
	
	if(x>xmaxl .or. z>zmaxl .or. x<xminl .or. z<zminl)then
		write(*,*) 'outofbounds',elnum
		stop
	else
!********************************************************************************
		if (elnum < 1 .or. elnum > dimj) then
			do i = 1,gridi-1
				if (xg(i)<=x .and. x<=xg(i+1))then
					exit
				end if
			end do
			do j = 1,gridj-1
				if (zg(j)<=z .and. z<=zg(j+1))then
					exit
				end if
			end do
		elseif (recallPos(elnum,1) == -1 .or. recallPos(elnum,2) == -1) then
			do i = 1,gridi-1
				if (xg(i)<=x .and. x<=xg(i+1))then
					exit
				end if
			end do
			do j = 1,gridj-1
				if (zg(j)<=z .and. z<=zg(j+1))then
					exit
				end if
			end do
			recallPos(elnum,1) = i
			recallPos(elnum,2) = j
		else
			if(xg(recallPos(elnum,1))<=x .and. x<=xg(recallPos(elnum,1)+1))then
				i = recallPos(elnum,1)
			elseif(recallPos(elnum,1)<gridi-1 .and. recallPos(elnum,1) > 1)then
				if(xg(recallPos(elnum,1)+1)<=x .and. x<=xg(recallPos(elnum,1)+2))then
					recallPos(elnum,1) = recallPos(elnum,1)+1
					i = recallPos(elnum,1)
				elseif(xg(recallPos(elnum,1)-1)<=x .and. x<=xg(recallPos(elnum,1)))then
					recallPos(elnum,1) = recallPos(elnum,1)-1
					i = recallPos(elnum,1)
				else
					do i = 1,gridi-1
						if (xg(i)<=x .and. x<=xg(i+1))then
							exit
						endif
					end do
					recallPos(elnum,1) = i
				endif
			elseif(recallPos(elnum,1)==gridi-1)then
				if(xg(recallPos(elnum,1)-1)<=x .and. x<=xg(recallPos(elnum,1)))then
					recallPos(elnum,1) = recallPos(elnum,1)-1
					i = recallPos(elnum,1)
				else
					do i = 1,gridi-1
						if (xg(i)<=x .and. x<=xg(i+1))then
							exit
						end if
					end do
					recallPos(elnum,1) = i
				endif
			elseif(recallPos(elnum,1)==1)then
				if(xg(recallPos(elnum,1)+1)<=x .and. x<=xg(recallPos(elnum,1)+2))then
					recallPos(elnum,1) = recallPos(elnum,1)+1
					i = recallPos(elnum,1)
				else
					do i = 1,gridi-1
						if (xg(i)<=x .and. x<=xg(i+1))then
							exit
						endif
					end do
					recallPos(elnum,1) = i
				endif
			else
				do i = 1,gridi-1
					if (xg(i)<=x .and. x<=xg(i+1))then
						exit
					end if
				end do
			endif
!***************************************************************************************
			if(zg(recallPos(elnum,2))<=z .and. z<=zg(recallPos(elnum,2)+1))then
				j = recallPos(elnum,2)
			elseif(recallPos(elnum,2)<gridj-1 .and. recallPos(elnum,2) > 1)then
				if(zg(recallPos(elnum,2)+1)<=z .and. z<=zg(recallPos(elnum,2)+2))then
					recallPos(elnum,2) = recallPos(elnum,2)+1
					j = recallPos(elnum,2)
				elseif(zg(recallPos(elnum,2)-1)<=z .and. z<=zg(recallPos(elnum,2)))then
					recallPos(elnum,2) = recallPos(elnum,2)-1
					j = recallPos(elnum,2)
				else
					do j = 1,gridj-1
						if (zg(j)<=z .and. z<=zg(j+1))then
							exit
						endif
					end do
					recallPos(elnum,2) = j
				endif
			elseif(recallPos(elnum,2)==gridj-1)then
				if(zg(recallPos(elnum,2)-1)<=z .and. z<=zg(recallPos(elnum,2)))then
					recallPos(elnum,2) = recallPos(elnum,2)-1
					j = recallPos(elnum,2)
				else
					do j = 1,gridj-1
						if (zg(j)<=z .and. z<=zg(j+1))then
							exit
						end if
					end do
					recallPos(elnum,2) = j
				endif
			elseif(recallPos(elnum,2)==1)then
				if(zg(recallPos(elnum,2)+1)<=z .and. z<=zg(recallPos(elnum,2)+2))then
					recallPos(elnum,2) = recallPos(elnum,2)+1
					j = recallPos(elnum,2)
				else
					do j = 1,gridj-1
						if (zg(j)<=z .and. z<=zg(j+1))then
							exit
						endif
					end do
					recallPos(elnum,2) = j
				endif
			else
				do j = 1,gridj-1
					if (zg(j)<=z .and. z<=zg(j+1))then
						exit
					end if
				end do
			endif
		endif
	endif
	
	select case(mesh)
		case('Ag')
			call getCorners(Ag(i:i+1,j:j+1),Ag_d(i:i+1,j:j+1,:),gc,gc_x,gc_z,gc_xz)
			call bcuint(gc,gc_x,gc_z,gc_xz,xg(i),xg(i+1),zg(j),zg(j+1),x,z,ansy,ansy_x,ansy_z)
			bicubicInterp = ansy
		case('Bxg')
			call getCorners(Bxg(i:i+1,j:j+1),Bxg_d(i:i+1,j:j+1,:),gc,gc_x,gc_z,gc_xz)
			call bcuint(gc,gc_x,gc_z,gc_xz,xg(i),xg(i+1),zg(j),zg(j+1),x,z,ansy,ansy_x,ansy_z)
			bicubicInterp = ansy	
		case('Bzg')
			call getCorners(Bzg(i:i+1,j:j+1),Bzg_d(i:i+1,j:j+1,:),gc,gc_x,gc_z,gc_xz)
			call bcuint(gc,gc_x,gc_z,gc_xz,xg(i),xg(i+1),zg(j),zg(j+1),x,z,ansy,ansy_x,ansy_z)
			bicubicInterp = ansy
		case('Pg')
			call getCorners(Pg(i:i+1,j:j+1),Pg_d(i:i+1,j:j+1,:),gc,gc_x,gc_z,gc_xz)
			call bcuint(gc,gc_x,gc_z,gc_xz,xg(i),xg(i+1),zg(j),zg(j+1),x,z,ansy,ansy_x,ansy_z)
			bicubicInterp = ansy
		case('Ng')
			call getCorners(Ng(i:i+1,j:j+1),Ng_d(i:i+1,j:j+1,:),gc,gc_x,gc_z,gc_xz)
			call bcuint(gc,gc_x,gc_z,gc_xz,xg(i),xg(i+1),zg(j),zg(j+1),x,z,ansy,ansy_x,ansy_z)
			bicubicInterp = ansy
		case('Ptotg')
			call getCorners(Ptotg(i:i+1,j:j+1),Ptotg_d(i:i+1,j:j+1,:),gc,gc_x,gc_z,gc_xz)
			call bcuint(gc,gc_x,gc_z,gc_xz,xg(i),xg(i+1),zg(j),zg(j+1),x,z,ansy,ansy_x,ansy_z)
			bicubicInterp = ansy
		case('dx_Bxg')
			call getCorners(Bxg(i:i+1,j:j+1),Bxg_d(i:i+1,j:j+1,:),gc,gc_x,gc_z,gc_xz)
			call bcuint(gc,gc_x,gc_z,gc_xz,xg(i),xg(i+1),zg(j),zg(j+1),x,z,ansy,ansy_x,ansy_z)
			bicubicInterp = ansy_x
		case('dz_Bxg')
			call getCorners(Bxg(i:i+1,j:j+1),Bxg_d(i:i+1,j:j+1,:),gc,gc_x,gc_z,gc_xz)
			call bcuint(gc,gc_x,gc_z,gc_xz,xg(i),xg(i+1),zg(j),zg(j+1),x,z,ansy,ansy_x,ansy_z)
			bicubicInterp = ansy_z	
		case('dx_Bzg')
			call getCorners(Bzg(i:i+1,j:j+1),Bzg_d(i:i+1,j:j+1,:),gc,gc_x,gc_z,gc_xz)
			call bcuint(gc,gc_x,gc_z,gc_xz,xg(i),xg(i+1),zg(j),zg(j+1),x,z,ansy,ansy_x,ansy_z)
			bicubicInterp = ansy_x
		case('dz_Bzg')
			call getCorners(Bzg(i:i+1,j:j+1),Bzg_d(i:i+1,j:j+1,:),gc,gc_x,gc_z,gc_xz)
			call bcuint(gc,gc_x,gc_z,gc_xz,xg(i),xg(i+1),zg(j),zg(j+1),x,z,ansy,ansy_x,ansy_z)
			bicubicInterp = ansy_z
		case('dx_Ptotg')
			call getCorners(Ptotg(i:i+1,j:j+1),Ptotg_d(i:i+1,j:j+1,:),gc,gc_x,gc_z,gc_xz)
			call bcuint(gc,gc_x,gc_z,gc_xz,xg(i),xg(i+1),zg(j),zg(j+1),x,z,ansy,ansy_x,ansy_z)
			bicubicInterp = ansy_x
		case('dz_Ptotg')
			call getCorners(Ptotg(i:i+1,j:j+1),Ptotg_d(i:i+1,j:j+1,:),gc,gc_x,gc_z,gc_xz)
			call bcuint(gc,gc_x,gc_z,gc_xz,xg(i),xg(i+1),zg(j),zg(j+1),x,z,ansy,ansy_x,ansy_z)
			bicubicInterp = ansy_z
		case('testMesh')
			call getCorners(testMesh(i:i+1,j:j+1),testMesh_d(i:i+1,j:j+1,:),gc,gc_x,gc_z,gc_xz)
			call bcuint(gc,gc_x,gc_z,gc_xz,xg(i),xg(i+1),zg(j),zg(j+1),x,z,ansy,ansy_x,ansy_z)
			bicubicInterp = ansy
		case default
			write(*,*) 'error: bicubicInterp case not found'
			stop
	end select
end function

subroutine getCorners(a,b,gc,gc_x,gc_z,gc_xz)
	implicit none
	real,dimension(4) :: gc,gc_x,gc_z,gc_xz
	real,dimension(2,2) :: a
	real,dimension(2,2,3) :: b
	
	gc(1) = a(1,1)
	gc(2) = a(2,1)
	gc(3) = a(2,2)
	gc(4) = a(1,2)
	
	gc_x(1) = b(1,1,1)
	gc_x(2) = b(2,1,1)
	gc_x(3) = b(2,2,1)
	gc_x(4) = b(1,2,1)
	
	gc_z(1) = b(1,1,2)
	gc_z(2) = b(2,1,2)
	gc_z(3) = b(2,2,2)
	gc_z(4) = b(1,2,2)
	
	gc_xz(1) = b(1,1,3)
	gc_xz(2) = b(2,1,3)
	gc_xz(3) = b(2,2,3)
	gc_xz(4) = b(1,2,3)
	
end subroutine
