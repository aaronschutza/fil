!*********************** Background Tracing Routines ***************************

subroutine tracerRK4(xo,yo,zo,idir,isopen,pos,ifirst,ilast)
        use dimM
        use stateIndex
        use potential
        use initCond
        use parametersBackground, only: useGrid,wallBound
        implicit none
        real :: xo,yo,zo
        integer :: idir,ilast,ifirst
        real,dimension(x:z,dimj) :: pos
        real :: Bx,Bz,dir
        integer :: i,j,k
        real,dimension(x:z) :: k1,k2,k3,k4
        logical :: xl,yl,zl,bl,isopen,tracePrint
        interface bn
                function bn (position)
                        use stateIndex
                        implicit none
                        real position(x:z)
                        real bn(x:z)
                end function bn
        end interface
        
        tracePrint = .false.
        
        
        
        step = 0.1
        
        do while(.true.)
                if(tracePrint) write(*,*) xo,zo
                
                dir = sign(1,idir)
                pos(x,ifirst) = xo
                pos(z,ifirst) = zo
                
                do i = ifirst,dimj
                
                        if(i==dimj)then
                                ilast = dimj+1
                                exit
                        endif
                
                        k1 = dir*step*bn(pos(:,i))
                        k2 = dir*step*bn(pos(:,i)+0.5*k1)
                        k3 = dir*step*bn(pos(:,i)+0.5*k2)
                        k4 = dir*step*bn(pos(:,i)+k3)
                        pos(:,i+1) = pos(:,i)+(k1+2.0*k2+2.0*k3+k4)/6.0
                        
                        !write(*,*) pos(:,i+1)
                        
                        xl = pos(x,i+1) > xmaxl .or. pos(x,i+1) < xminl
                        yl = .false.
                        zl = pos(z,i+1) > zmaxl .or. pos(z,i+1) < zminl
                        
                        if(wallBound)then
                                bl = abs(pos(x,i+1)) < boundary
                        else
                                bl = sqrt(pos(x,i+1)**2+pos(z,i+1)**2) < boundary
                        endif
                        
                        if(bl .or. xl .or. yl .or. zl)then
                                ilast = i+1
                                exit
                        endif
                enddo
                
                if(wallBound)then        
                        if(ilast<dimj)then
                                if(tracePrint) write(*,*) 'step too big'
                                step = (ilast*step)/dimj
                                if(tracePrint) write(*,*) sqrt(pos(x,ilast)**2+pos(z,ilast)**2)
                        endif
        
                        if(ilast == dimj)then
                                if(bl .and. abs(abs(pos(x,ilast))-boundary)<1.0e-8)then
                                        isopen = .false.
                                        !write(*,*) 'found endpoint'
                                        !write(*,*) sqrt(pos(x,dimj)**2+pos(z,dimj)**2)
                                        exit
                                else
                                        step = step - abs(abs(pos(x,ilast))-boundary)/dimj
                                endif
                                if(xl .or. yl .or. zl)then
                                        isopen = .true.
                                        exit
                                endif
                        
                        endif
        
                        if(ilast == dimj+1)then
                                if(tracePrint) write(*,*) 'step too small'
                                step = step + abs(abs(pos(x,ilast))-boundary)/dimj
                                if(tracePrint) write(*,*) sqrt(pos(x,dimj)**2+pos(z,dimj)**2)
                        endif
                
                else                
                        if(ilast<dimj)then
                                if(tracePrint) write(*,*) 'step too big'
                                step = (ilast*step)/dimj
                                if(tracePrint) write(*,*) sqrt(pos(x,ilast)**2+pos(z,ilast)**2)
                        endif
        
                        if(ilast == dimj)then
                                if(bl .and. abs(sqrt(pos(x,ilast)**2+pos(z,ilast)**2)-boundary)<1.0e-8)then
                                        isopen = .false.
                                        !write(*,*) 'found endpoint'
                                        !write(*,*) sqrt(pos(x,dimj)**2+pos(z,dimj)**2)
                                        exit
                                else
                                        step = step - abs(sqrt(pos(x,ilast)**2+pos(z,ilast)**2)-boundary)/dimj
                                endif
                                if(xl .or. yl .or. zl)then
                                        isopen = .true.
                                        exit
                                endif
                        
                        endif
        
                        if(ilast == dimj+1)then
                                if(tracePrint) write(*,*) 'step too small'
                                step = step + abs(sqrt(pos(x,ilast)**2+pos(z,ilast)**2)-boundary)/dimj
                                if(tracePrint) write(*,*) sqrt(pos(x,dimj)**2+pos(z,dimj)**2)
                        endif
                endif
        enddo
        
end subroutine tracerRK4


subroutine tracerRK4a(xo,yo,zo,idir,isopen,pos,ifirst,ilast)
        use dimM
        use stateIndex
        use potential
        use initCond
        use parametersBackground, only: useGrid,wallBound
        implicit none
        real :: xo,yo,zo
        integer :: idir,ilast,ifirst
        real,dimension(x:z,dimj) :: pos
        real :: Bx,Bz,dir,Bmid,scale,length
        integer :: i,j,k
        real,dimension(x:z) :: k1,k2,k3,k4
        logical :: xl,yl,zl,bl,isopen,tracePrint
        interface bn
                function bn (position)
                        use stateIndex
                        implicit none
                        real position(x:z)
                        real bn(x:z)
                end function bn
        end interface
        
        tracePrint = .true.
        
        
        
        step = 0.01
        
        do while(.true.)
                if(tracePrint) write(*,*) xo,zo
                
                dir = sign(1,idir)
                pos(x,ifirst) = xo
                pos(z,ifirst) = zo
                Bmid = sqrt(Bx(xo,zo,-1)**2+Bz(xo,zo,-1)**2)
                do i = ifirst,dimj
                
                        if(i==dimj)then
                                ilast = dimj+1
                                exit
                        endif
                        scale = sqrt(Bx(pos(x,i),pos(z,i),-1)**2+Bz(pos(x,i),pos(z,i),-1)**2)/Bmid
                        

                        k1 = dir*scale*step*bn(pos(:,i))
                        
                        xl = pos(x,i)+0.5*k1(x) > xmaxl .or. pos(x,i)+0.5*k1(x) < xminl
                        zl = pos(z,i)+0.5*k1(z) > zmaxl .or. pos(z,i)+0.5*k1(z) < zminl
                        if(xl .or. zl)then
                                ilast = i+1
                                if (xl .and. tracePrint) write(11,*) 'hit x wall ',pos(:,i)+k1
                                if (zl .and. tracePrint) write(11,*) 'hit z wall ',pos(:,i)+k1
                                exit
                        endif
                        
                        k2 = dir*scale*step*bn(pos(:,i)+0.5*k1)
                        xl = pos(x,i)+0.5*k2(x) > xmaxl .or. pos(x,i)+0.5*k2(x) < xminl
                        zl = pos(z,i)+0.5*k2(z) > zmaxl .or. pos(z,i)+0.5*k2(z) < zminl
                        if(xl .or. zl)then
                                ilast = i+1
                                if (xl .and. tracePrint) write(11,*) 'hit x wall ',pos(:,i)+k2
                                if (zl .and. tracePrint) write(11,*) 'hit z wall ',pos(:,i)+k2
                                exit
                        endif
                        
                        k3 = dir*scale*step*bn(pos(:,i)+0.5*k2)
                        xl = pos(x,i)+k3(x) > xmaxl .or. pos(x,i)+k3(x) < xminl
                        zl = pos(z,i)+k3(z) > zmaxl .or. pos(z,i)+k3(z) < zminl
                        if(xl .or. zl)then
                                ilast = i+1
                                if (xl .and. tracePrint) write(11,*) 'hit x wall ',pos(:,i)+k3
                                if (zl .and. tracePrint) write(11,*) 'hit z wall ',pos(:,i)+k3
                                exit
                        endif
                        
                        k4 = dir*scale*step*bn(pos(:,i)+k3)
                        pos(:,i+1) = pos(:,i)+(k1+2.0*k2+2.0*k3+k4)/6.0
                        xl = pos(x,i+1) > xmaxl .or. pos(x,i+1) < xminl
                        zl = pos(z,i+1) > zmaxl .or. pos(z,i+1) < zminl
                        if(xl .or. zl)then
                                ilast = i+1
                                if (xl .and. tracePrint) write(11,*) 'hit x wall ',pos(:,i+1)
                                if (zl .and. tracePrint) write(11,*) 'hit z wall ',pos(:,i+1)
                                exit
                        endif
                        

                        pos(:,i+1) = pos(:,i)+(k1+2.0*k2+2.0*k3+k4)/6.0
                        
                        !write(*,*) pos(:,i+1)
                        
                        xl = pos(x,i+1) > xmaxl .or. pos(x,i+1) < xminl
                        yl = .false.
                        zl = pos(z,i+1) > zmaxl .or. pos(z,i+1) < zminl
                        
                        if(wallBound)then
                                bl = abs(pos(x,i+1)) < boundary
                        else
                                bl = sqrt(pos(x,i+1)**2+pos(z,i+1)**2) < boundary
                        endif
                        
                        if(bl .or. xl .or. yl .or. zl)then
                                ilast = i+1
                                exit
                        endif
                enddo
                
                if(wallBound)then        
                        if(ilast<dimj)then
                                if(tracePrint) write(*,*) 'step too big'
                                step = (ilast*step)/dimj
                                if(tracePrint) write(*,*) sqrt(pos(x,ilast)**2+pos(z,ilast)**2)
                        endif
        
                        if(ilast == dimj)then
                                if(bl .and. abs(abs(pos(x,ilast))-boundary)<1.0e-8)then
                                        isopen = .false.
                                        !write(*,*) 'found endpoint'
                                        !write(*,*) sqrt(pos(x,dimj)**2+pos(z,dimj)**2)
                                        exit
                                else
                                        step = step - abs(abs(pos(x,ilast))-boundary)/dimj
                                endif
                                if(xl .or. yl .or. zl)then
                                        isopen = .true.
                                        exit
                                endif
                        
                        endif
        
                        if(ilast == dimj+1)then
                                if(tracePrint) write(*,*) 'step too small'
                                step = step + abs(abs(pos(x,ilast))-boundary)/dimj
                                if(tracePrint) write(*,*) sqrt(pos(x,dimj)**2+pos(z,dimj)**2)
                        endif
                
                else                
                        if(ilast<dimj)then
                                if(tracePrint) write(*,*) 'step too big',step
                                length = 0.0
                                do i = ifirst,ilast-1
                                        length = length+sqrt((pos(x,i)-pos(x,i+1))**2+(pos(z,i)-pos(z,i+1))**2)
                                enddo
                                step = step*(length/(length+abs(sqrt(pos(x,ilast)**2+pos(z,ilast)**2)-boundary)))
                                if(tracePrint) write(*,*) sqrt(pos(x,ilast)**2+pos(z,ilast)**2)
                        endif
        
                        if(ilast == dimj)then
                                if(bl .and. abs(sqrt(pos(x,ilast)**2+pos(z,ilast)**2)-boundary)<1.0e-8)then
                                        isopen = .false.
                                        !write(*,*) 'found endpoint'
                                        !write(*,*) sqrt(pos(x,dimj)**2+pos(z,dimj)**2)
                                        exit
                                else
                                        if(tracePrint) write(*,*) 'step too big 2',step
                                        step = step - 1.0e-5*abs(sqrt(pos(x,ilast)**2+pos(z,ilast)**2)-boundary)
                                        if(tracePrint) write(*,*) sqrt(pos(x,dimj)**2+pos(z,dimj)**2)
                                endif
                                if(xl .or. yl .or. zl)then
                                        isopen = .true.
                                        exit
                                endif
                        
                        endif
        
                        if(ilast == dimj+1)then
                                if(tracePrint) write(*,*) 'step too small',step
                                length = 0.0
                                do i = ifirst,ilast-1
                                        length = length+sqrt((pos(x,i)-pos(x,i+1))**2+(pos(z,i)-pos(z,i+1))**2)
                                enddo
                                step = step*(abs(sqrt(pos(x,ilast)**2+pos(z,ilast)**2)-boundary)/length+1.0)
                                if(tracePrint) write(*,*) sqrt(pos(x,dimj)**2+pos(z,dimj)**2)
                        endif
                endif
        enddo
        
end subroutine tracerRK4a

subroutine tracerRK4pres(i_in,j_in,xo,zo,pgl,pgx,pgz,pres,openPoint)
        use dimM
        use stateIndex
        use potential
        use initCond
        implicit none
        real :: xo,zo,xi,zi,l1,l2,pres
        integer :: idir,ilast,ifirst
        integer :: posDimj = 2000
        real,dimension(x:z,2000) :: pos
        real :: dir = 1.0
        integer :: i,j,k,numTries,i_in,j_in
        real,dimension(x:z) :: k1,k2,k3,k4
        logical :: xl,zl,bl
        logical :: changeDir
        logical :: openPoint
        logical :: mssg = .false.
        real :: bgpress,innerPress,cutRadius
        real,dimension(plasGridi) :: pgl,pgx,pgz
        interface bn
                function bn (position)
                        use stateIndex
                        implicit none
                        real position(x:z)
                        real bn(x:z)
                end function bn
        end interface
        
        if(mssg) write(11,*) 'starting on point',i_in,j_in,xo,zo
        
        bgpress = pgl(1)
        innerPress = pgl(plasGridi)
        cutRadius = sqrt(pgx(plasGridi)**2+pgz(plasGridi)**2)
        step = 0.1
        numTries = 0
        xl = .false.
        zl = .false.
        bl = .false.
        ifirst = 1
        changeDir = .true.
        idir = 1
        

        
        if(sqrt(xo**2+zo**2)<=cutRadius)then
                pres = innerPress
                openPoint = .true.
                return
        endif
        
        
        
        do while(.true.)
                
                numTries = numTries+1
                if (numTries>4) then
                        pres = bgpress
                        openPoint = .true.
                        if(mssg) write(11,*) 'failed to trace to plasGrid'
                        return
                endif
                
                dir = sign(1,idir)
                pos(x,ifirst) = xo
                pos(z,ifirst) = zo
                
                do i = ifirst,posDimj
                        if(i==posDimj)then
                                ilast = posDimj+1
                                exit
                        endif
                        
                        k1 = dir*step*bn(pos(:,i))
                        
                        
                        xl = pos(x,i)+0.5*k1(x) > xmaxl .or. pos(x,i)+0.5*k1(x) < xminl
                        zl = pos(z,i)+0.5*k1(z) > zmaxl .or. pos(z,i)+0.5*k1(z) < zminl
                        if(xl .or. zl)then
                                ilast = i+1
                                if (xl .and. mssg) write(11,*) 'hit x wall ',pos(:,i)+k1
                                if (zl .and. mssg) write(11,*) 'hit z wall ',pos(:,i)+k1
                                exit
                        endif
                        
                        k2 = dir*step*bn(pos(:,i)+0.5*k1)
                        xl = pos(x,i)+0.5*k2(x) > xmaxl .or. pos(x,i)+0.5*k2(x) < xminl
                        zl = pos(z,i)+0.5*k2(z) > zmaxl .or. pos(z,i)+0.5*k2(z) < zminl
                        if(xl .or. zl)then
                                ilast = i+1
                                if (xl .and. mssg) write(11,*) 'hit x wall ',pos(:,i)+k2
                                if (zl .and. mssg) write(11,*) 'hit z wall ',pos(:,i)+k2
                                exit
                        endif
                        
                        k3 = dir*step*bn(pos(:,i)+0.5*k2)
                        xl = pos(x,i)+k3(x) > xmaxl .or. pos(x,i)+k3(x) < xminl
                        zl = pos(z,i)+k3(z) > zmaxl .or. pos(z,i)+k3(z) < zminl
                        if(xl .or. zl)then
                                ilast = i+1
                                if (xl .and. mssg) write(11,*) 'hit x wall ',pos(:,i)+k3
                                if (zl .and. mssg) write(11,*) 'hit z wall ',pos(:,i)+k3
                                exit
                        endif
                        
                        k4 = dir*step*bn(pos(:,i)+k3)
                        pos(:,i+1) = pos(:,i)+(k1+2.0*k2+2.0*k3+k4)/6.0
                        xl = pos(x,i+1) > xmaxl .or. pos(x,i+1) < xminl
                        zl = pos(z,i+1) > zmaxl .or. pos(z,i+1) < zminl
                        if(xl .or. zl)then
                                ilast = i+1
                                if (xl .and. mssg) write(11,*) 'hit x wall ',pos(:,i+1)
                                if (zl .and. mssg) write(11,*) 'hit z wall ',pos(:,i+1)
                                exit
                        endif
                        
                        
                        if(mssg) write(11,*) 'posx',pos(x,i),'posz',pos(z,i),sqrt((pos(x,i-1)-pos(x,i))**2+(pos(z,i-1)-pos(z,i))**2)
                        
                        if(sqrt((pos(x,i))**2+(pos(z,i))**2)<=cutRadius)then
                                ilast = i+1
                                if(mssg) write(11,*) 'hit cut off at r=',cutRadius
                                exit
                        endif
                enddo
                
                do i = 1,ilast-2
                  do j = 1,plasGridi-1
                    call line_intersects(pos(x,i),pos(z,i),pos(x,i+1),pos(z,i+1),pgx(j),pgz(j),pgx(j+1),pgz(j+1),bl,xi,zi,i,j,mssg)
                        if(bl)then
                                if(mssg) write(11,*) 'found an intercept at',xi,zi
                                if(mssg) write(11,*) 'num tries',numTries
                                if(mssg) then
                                        write(11,*) 'pgrid point 1:',pgx(j),pgz(j)
                                        write(11,*) 'pgrid point 2:',pgx(j+1),pgz(j+1)
                                        write(11,*) 'pos grid 1:',pos(x,i),pos(z,i)
                                        write(11,*) 'pos grid 2:',pos(x,i+1),pos(z,i+1)
                                endif
                                l1 = sqrt((xi-pgx(j))**2+(zi-pgz(j))**2)
                                l2 = sqrt((xi-pgx(j+1))**2+(zi-pgz(j+1))**2)
                                pres = l1*(pgl(j+1)-pgl(j))/(l1+l2) + pgl(j)
                                openPoint = .false.
                                return
                        endif
                  enddo
                enddo
                
                if(ilast == posDimj+1)then
                        if(mssg) write(11,*) 'too short'
                        if(ChangeDir)then
                                idir = -idir
                                if(mssg) write(11,*) 'flipping'
                                ChangeDir = .false.
                        else
                                step = 2.0*step
                                openPoint = .true.
                                if(mssg) write(11,*) 'flipped and failed'
                                ChangeDir = .true.
                        endif
                endif
                
                if(ilast<=posDimj)then
                        if(mssg) write(11,*) 'hit barrier'
                        if(ChangeDir)then
                                idir = -idir
                                if(mssg) write(11,*) 'flipping'
                                ChangeDir = .false.
                        else
                                openPoint = .true.
                                if(mssg) write(11,*) 'flipped and failed'
                                ChangeDir = .true.
                        endif
                endif
        enddo
        
        
end subroutine tracerRK4pres


subroutine line_intersects(p0_x,p0_y,p1_x,p1_y,p2_x,p2_y,p3_x,p3_y,bl,x_i,y_i,index_1,index_2,mssg)
        implicit none
        real :: p0_x, p0_y, p1_x, p1_y, p2_x, p2_y, p3_x, p3_y, x_i, y_i
        logical :: bl,mssg
        real :: s1_x, s1_y, s2_x, s2_y
        real :: s, t
        integer :: index_1,index_2
        s1_x = p1_x - p0_x
        s1_y = p1_y - p0_y
        s2_x = p3_x - p2_x
        s2_y = p3_y - p2_y
 
        s = (-s1_y*(p0_x-p2_x)+s1_x*(p0_y-p2_y))/(-s2_x*s1_y+s1_x*s2_y)
        t = ( s2_x*(p0_y-p2_y)-s2_y*(p0_x-p2_x))/(-s2_x*s1_y+s1_x*s2_y)
 
        if (s >= 0 .and. s <= 1 .and. t >= 0 .and. t <= 1) then
                !if(mssg) write(*,*) 'Collision detected at:'
                !if(mssg) write(*,*) 'index line 1:',index_1
                !if(mssg) write(*,*) 'index line 2:',index_2
                bl = .true.
                x_i = p0_x+(t*s1_x)
                y_i = p0_y+(t*s1_y)
        else
                !No collision
                bl = .false.
                x_i = 0.0
                y_i = 0.0
        endif
        
end subroutine


function bn(pos)
        use stateIndex
        implicit none
        real,dimension(x:z) :: bn,pos
        real :: Bx,Bz
        real :: c0,c1,c2

        c0 = Bx(pos(x),pos(z),-1)
        c1 = Bz(pos(x),pos(z),-1)
        c2 = sqrt(c0**2+c1**2)
        
        if (c2<=0.0) then
                bn=0.0
        else
                bn =  [c0,c1]/c2
        endif
!        interface bn
!                function bn (position)
!                        use stateIndex
!                        implicit none
!                        real position(x:z)
!                        real bn(x:z)
!                end function bn
!        end interface
        
end function bn
subroutine find_equator_tracerRK4a(xo,yo,zo,idir,isopen,xe)
        use dimM
        use stateIndex
        use potential
        use initCond
        use parametersBackground, only: useGrid,wallBound
        implicit none
        real :: xo,yo,zo
        integer :: idir,ilast,ifirst
        real,dimension(x:z,dimj) :: pos
        real :: Bx,Bz,dir,dz,xe
        integer :: i,j,k
        real,dimension(x:z) :: k1,k2,k3,k4
        logical :: xl,yl,zl,bl,isopen,tracePrint
        logical, parameter :: debug = .true.
        interface bn
                function bn (position)
                        use stateIndex
                        implicit none
                        real position(x:z)
                        real bn(x:z)
                end function bn
        end interface
        
        tracePrint = .false.
        
        ifirst = 1
        ilast = dimj
        xe=100.
         
        
        step = 0.01
        
        do while(.true.)
                if(tracePrint) write(*,*) xo,zo
                
                dir = sign(1,idir)
                pos(x,ifirst) = xo
                pos(z,ifirst) = zo
                
                do i = ifirst,dimj
                
                        if(i==dimj)then
                                ilast = dimj+1
                                exit
                        endif
                
                        k1 = dir*step*bn(pos(:,i))
                        k2 = dir*step*bn(pos(:,i)+0.5*k1)
                        k3 = dir*step*bn(pos(:,i)+0.5*k2)
                        k4 = dir*step*bn(pos(:,i)+k3)
                        pos(:,i+1) = pos(:,i)+(k1+2.0*k2+2.0*k3+k4)/6.0
                        
                        !write(*,*) pos(:,i+1)
                        
                        xl = pos(x,i+1) > xmaxl .or. pos(x,i+1) < xminl
                        yl = .false.
                        zl = pos(z,i+1) > zmaxl .or. pos(z,i+1) < zminl
                        
                        if(wallBound)then
                                bl = abs(pos(x,i+1)) < boundary
                        else
                                bl = sqrt(pos(x,i+1)**2+pos(z,i+1)**2) < boundary
                        endif
                        
                        if(bl .or. xl .or. yl .or. zl)then
                                ilast = i+1
                                exit
                        endif
                                       
                        if(pos(z,i+1)*pos(z,i) <=0.0)then
                          dz = pos(z,i+1)-pos(z,i)
                          if(abs(dz)>0.0) xe = (pos(x,i) *  pos(z,i+1) -  pos(x,i+1) *  pos(z,i))/dz
                          if(debug)then
                                  write(6,'(3(a,f10.4))')' xo=',xo,' zo=',zo,' xequator = ',xe
                          endif
                          return
                        endif


                enddo
                
                if(wallBound)then        
                        if(ilast<dimj)then
                                if(tracePrint) write(*,*) 'step too big'
                                step = (ilast*step)/dimj
                                if(tracePrint) write(*,*) sqrt(pos(x,ilast)**2+pos(z,ilast)**2)
                        endif
        
                        if(ilast == dimj)then
                                if(bl .and. abs(abs(pos(x,ilast))-boundary)<1.0e-8)then
                                        isopen = .false.
                                        if(debug)then
                                          write(*,*) 'found endpoint'
                                          write(*,*) 'r =',sqrt(pos(x,i)**2+pos(z,i)**2)
                                          write(6,*)' xo =',xo,'zo =',zo
                                          write(6,*)' x =',pos(x,i),' z =',pos(z,i)
                                          write(6,*)' i=',i
                                        end if
                                        exit
                                else
                                        step = step - abs(abs(pos(x,ilast))-boundary)/dimj
                                endif
                                if(xl .or. yl .or. zl)then
                                        isopen = .true.
                                        exit
                                endif
                        
                        endif
        
                        if(ilast == dimj+1)then
                                if(tracePrint) write(*,*) 'step too small'
                                step = step + abs(abs(pos(x,ilast))-boundary)/dimj
                                if(tracePrint) write(*,*) sqrt(pos(x,dimj)**2+pos(z,dimj)**2)
                        endif
                
                else                
                        if(ilast<dimj)then
                                if(tracePrint) write(*,*) 'step too big'
                                step = (ilast*step)/dimj
                                if(tracePrint) write(*,*) sqrt(pos(x,ilast)**2+pos(z,ilast)**2)
                        endif
        
                        if(ilast == dimj)then
                                if(bl .and. abs(sqrt(pos(x,ilast)**2+pos(z,ilast)**2)-boundary)<1.0e-8)then
                                        isopen = .false.
                                        if(debug)then
                                          write(*,*) 'found endpoint'
                                          write(*,*) 'r =',sqrt(pos(x,i)**2+pos(z,i)**2)
                                          write(6,*)' xo =',xo,'zo =',zo
                                          write(6,*)' x =',pos(x,i),' z =',pos(z,i)
                                          write(6,*)' i=',i
                                        end if
                                        exit
                                else
                                        step = step - abs(sqrt(pos(x,ilast)**2+pos(z,ilast)**2)-boundary)/dimj
                                endif
                                if(xl .or. yl .or. zl)then
                                        isopen = .true.
                                        exit
                                endif
                        
                        endif
        
                        if(ilast == dimj+1)then
                                if(tracePrint) write(*,*) 'step too small'
                                step = step + abs(sqrt(pos(x,ilast)**2+pos(z,ilast)**2)-boundary)/dimj
                                if(tracePrint) write(*,*) sqrt(pos(x,dimj)**2+pos(z,dimj)**2)
                        endif
                endif
         enddo
       
        

        
end subroutine find_equator_tracerRK4a
