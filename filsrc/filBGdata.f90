!***************** Background Data Routines *************************

subroutine saveGridToFile(t,a,fileName,gridi,gridj)
       implicit none
       integer :: gridi,gridj
       real,dimension(gridi,gridj) :: a
       integer :: i,j,bufsize
       character(len=20) :: fileName
       real :: t
       
       write(*,*) 'file grid size',gridi,gridj
       
       inquire(Iolength=bufsize) a
       
       open(unit=12, file= trim(fileName), action="write", status="replace",Recl=bufsize)
         write(12,*) t
         write(12,*) gridi
         write(12,*) gridj
         do j = 1,gridj
          do i = 1,gridi
              if(i==gridi)then
                      write(12,"(f12.4)") a(i,j)
              else
                     write(12,"(f12.4,a)",advance='no') a(i,j),' '
              end if
          end do
         end do
       close(unit=12)
end subroutine

subroutine readTecFileRCME(fileName)
       use constants
       use potential
       use parametersBackground
       implicit none
       integer :: i,j,k,io
       real :: c0,c1,c2,c3,c4,c5
       real :: xin,yin,zin,bxin,byin,bzin,pin,rhoin,vxin,vyin,vzin,fcpin
       character(len=*) :: fileName
       character(len=100) :: lineRead
       integer :: ogi,ogj,ogk,yslice
       real :: dipx,dipy,dipz
       logical :: binary_exists
       
       write(*,*) 'reading file:'
       write(*,*) fileName
       
       xminl = 0.0
       xmaxl = 0.0
       yminl = 0.0
       ymaxl = 0.0
       zminl = 0.0
       zmaxl = 0.0
       inquire(file=fileName//'.bin',exist=binary_exists)
       if(binary_exists)then
              open(unit=11, file= fileName//'.bin',form="unformatted",status="old")
                     read(11) ogi
                     read(11) ogj
                     read(11) ogk
                     allocate(oPg(ogi,ogj,ogk))
                     allocate(oBxg(ogi,ogj,ogk))
                     allocate(oByg(ogi,ogj,ogk))
                     allocate(oBzg(ogi,ogj,ogk))
                     allocate(oxg(ogi))
                     allocate(oyg(ogj))
                     allocate(ozg(ogk))
                     read(11) oBxg
                     read(11) oByg
                     read(11) oBzg
                     read(11) oPg
                     read(11) oxg
                     read(11) oyg
                     read(11) ozg
              close(11)
       else
              open(unit=11, file= fileName, status="old")
              do i = 1,8
                     read(11,*) lineRead
              enddo
              
              do       
                     read(11,*,IOSTAT=io) i,j,k,xin,yin,zin,bxin,byin,bzin,pin,rhoin,vxin,vyin,vzin,fcpin
                     if(io>0)then
                            write(*,*) 'error in reading tec file'
                            exit
                     elseif(io<0)then
                            exit
                     else
                            gridi = max(i,gridi)
                            gridj = max(j,gridj)
                            gridk = max(k,gridk)
                     endif
              enddo
              close(11)
              
              open(unit=11, file= fileName, status="old")
              do i = 1,8
                     read(11,*) lineRead
              enddo
              allocate(oPg(gridi,gridj,gridk))
              allocate(oBxg(gridi,gridj,gridk))
              allocate(oByg(gridi,gridj,gridk))
              allocate(oBzg(gridi,gridj,gridk))
              allocate(oxg(gridi))
              allocate(oyg(gridj))
              allocate(ozg(gridk))
              do
                     !VARIABLES=          I J K  X  Y  Z(Re) Bx  By   Bz(nT)P(nP) rho Vx  Vy   Vz    FCP
                     read(11,*,IOSTAT=io) i,j,k,xin,yin,zin,bxin,byin,bzin,pin,rhoin,vxin,vyin,vzin,fcpin
 ! NOTE - conversion to nT and nPa                    
                     bxin = bxin*1.0E9
                     byin = byin*1.0E9
                     bzin = bzin*1.0E9
                     pin = pin*1.0E9
                     
                     if(io>0)then
                            write(*,*) 'error in reading tec file'
                            exit
                     elseif(io<0)then
                            exit
                     else
                            oxg(i) = xin
                            oyg(j) = yin
                            ozg(k) = zin
                            oBxg(i,j,k) = bxin
                            oByg(i,j,k) = byin
                            oBzg(i,j,k) = bzin
                            oPg(i,j,k) = pin
                     endif
              enddo
              close(11)
              ogi = gridi
              ogj = gridj
              ogk = gridk
              open(unit=11, file= fileName//'.bin',form="unformatted" ,status="new")
                     write(11) ogi
                     write(11) ogj
                     write(11) ogk
                     write(11) oBxg
                     write(11) oByg
                     write(11) oBzg
                     write(11) oPg
                     write(11) oxg
                     write(11) oyg
                     write(11) ozg
              close(11)
       endif
       

       yslice = (ogj+1)/2

       gridi = ogi
       gridj = ogk-2 + ogk-1
       
       
       allocate(Pg(gridi,gridj))
       allocate(Bxg(gridi,gridj))
       allocate(Byg(gridi,gridj))
       allocate(Bzg(gridi,gridj))
       allocate(xg(gridi))
       allocate(zg(gridj))
       
       xg = oxg
       
       do j = 1,gridj/2
              zg(j) = -ozg(ogk-j+1)
       enddo
       zg(gridj/2+1) = ozg(2)
       
       zg(gridj/2+2:gridj) = ozg(3:ogk)
       
       do j = 1,gridj/2
              Pg(:,j) = oPg(:,yslice,ogk-j+1)
       enddo
       Pg(:,gridj/2+1) = oPg(:,yslice,2)
       
       Pg(:,gridj/2+2:gridj) = oPg(:,yslice,3:ogk)
       
       do j = 1,gridj/2
              Bxg(:,j) = -oBxg(:,yslice,ogk-j+1)
       enddo
       Bxg(:,gridj/2+1) = oBxg(:,yslice,2)
       
       Bxg(:,gridj/2+2:gridj) = oBxg(:,yslice,3:ogk)
       
       do j = 1,gridj/2
              Bzg(:,j) = oBzg(:,yslice,ogk-j+1)
       enddo
       Bzg(:,gridj/2+1) = oBzg(:,yslice,2)
       
       Bzg(:,gridj/2+2:gridj) = oBzg(:,yslice,3:ogk)
       
       xminl = xg(1)
       xmaxl = xg(gridi)
       zminl = zg(1)
       zmaxl = zg(gridj)
       
       write(*,*) xminl,xmaxl
       write(*,*) zminl,zmaxl
       
       !write(*,*) gridj/2,gridj/2+1,gridj/2+2
       !write(*,*) xg(1),xg(gridi)
       !write(*,*) zg(1),zg(gridj/2),zg(gridj/2+1),zg(gridj/2+2),zg(gridj)
       !write(*,*) 'done loading, stopping'
       !write(*,*) Bzg(gridi,1),Bzg(gridi,gridj/2),Bzg(gridi,gridj/2+1),Bzg(gridi,gridj/2+2),Bzg(gridi,gridj)
       !write(*,*) Bxg(gridi,1),Bxg(gridi,gridj/2),Bxg(gridi,gridj/2+1),Bxg(gridi,gridj/2+2),Bxg(gridi,gridj)
       !write(*,*) Pg(gridi,1),Pg(gridi,gridj/2),Pg(gridi,gridj/2+1),Pg(gridi,gridj/2+2),Pg(gridi,gridj)
       !stop
       
       allocate(Ptotg(gridi,gridj))
       
       Ptotg = (Bxg**2+Bzg**2)/(2.0*mu0) + Pg
       
       
       deallocate(oPg)
       deallocate(oBxg)
       deallocate(oByg)
       deallocate(oBzg)
       deallocate(oxg)
       deallocate(oyg)
       deallocate(ozg)
end subroutine readTecFileRCME


subroutine readTecFileRCME2d(fileName)
       use constants
       use potential
       use parametersBackground
       implicit none
       integer :: i,j,k,io
       real :: c0,c1,c2,c3,c4,c5
       real :: xin,yin,zin,bxin,byin,bzin,pin,rhoin,vxin,vyin,vzin,fcpin
       character(len=*) :: fileName
       character(len=100) :: lineRead
       integer :: ogi,ogj,ogk,yslice
       real :: dipx,dipy,dipz
       logical :: binary_exists
       
       write(*,*) 'reading file:'
       write(*,*) fileName
       
       xminl = 0.0
       xmaxl = 0.0
       yminl = 0.0
       ymaxl = 0.0
       zminl = 0.0
       zmaxl = 0.0
       inquire(file=fileName//'.bin',exist=binary_exists)
       if(binary_exists)then
              open(unit=11, file= fileName//'.bin',form="unformatted",status="old")
                     read(11) ogi
                     read(11) ogj
                     read(11) ogk
                     allocate(oPg(ogi,ogj,ogk))
                     allocate(oBxg(ogi,ogj,ogk))
                     allocate(oByg(ogi,ogj,ogk))
                     allocate(oBzg(ogi,ogj,ogk))
                     allocate(oxg(ogi))
                     allocate(oyg(ogj))
                     allocate(ozg(ogk))
                     read(11) oBxg
                     read(11) oByg
                     read(11) oBzg
                     read(11) oPg
                     read(11) oxg
                     read(11) oyg
                     read(11) ozg
              close(11)
       else
              open(unit=11, file= fileName, status="old")
              do i = 1,8
                     read(11,*) lineRead
              enddo
              !VARIABLES="I" "K" "X(Re)" "Z(Re)"  "Bx" "By" "Bz" "P" 
              do       
                     read(11,*,IOSTAT=io) i,k,xin,zin,bxin,byin,bzin,pin
                     if(io>0)then
                            write(*,*) 'error in reading tec file'
                            exit
                     elseif(io<0)then
                            exit
                     else
                            gridi = max(i,gridi)
                            gridk = max(k,gridk)
                     endif
              enddo
              close(11)
              gridj = 1
              j = 1
              open(unit=11, file= fileName, status="old")
              do i = 1,8
                     read(11,*) lineRead
              enddo
              allocate(oPg(gridi,gridj,gridk))
              allocate(oBxg(gridi,gridj,gridk))
              allocate(oByg(gridi,gridj,gridk))
              allocate(oBzg(gridi,gridj,gridk))
              allocate(oxg(gridi))
              allocate(oyg(gridj))
              allocate(ozg(gridk))
              do
                     !VARIABLES="I" "K" "X(Re)" "Z(Re)"  "Bx" "By" "Bz" "P" 
                     read(11,*,IOSTAT=io) i,k,xin,zin,bxin,byin,bzin,pin
                     
                     bxin = bxin*1.0E9
                     byin = byin*1.0E9
                     bzin = bzin*1.0E9
                     pin = pin*1.0E9
                     
                     if(io>0)then
                            write(*,*) 'error in reading tec file'
                            exit
                     elseif(io<0)then
                            exit
                     else
                            oxg(i) = xin
                            ozg(k) = zin
                            oBxg(i,j,k) = bxin
                            oByg(i,j,k) = byin
                            oBzg(i,j,k) = bzin
                            oPg(i,j,k) = pin
                     endif
              enddo
              close(11)
              oyg = 0.0
              ogi = gridi
              ogj = gridj
              ogk = gridk
              open(unit=11, file= fileName//'.bin',form="unformatted" ,status="new")
                     write(11) ogi
                     write(11) ogj
                     write(11) ogk
                     write(11) oBxg
                     write(11) oByg
                     write(11) oBzg
                     write(11) oPg
                     write(11) oxg
                     write(11) oyg
                     write(11) ozg
              close(11)
       endif
       

       yslice = 1

       gridi = ogi
       gridj = ogk-2 + ogk-1
       
       
       allocate(Pg(gridi,gridj))
       allocate(Bxg(gridi,gridj))
       allocate(Byg(gridi,gridj))
       allocate(Bzg(gridi,gridj))
       allocate(xg(gridi))
       allocate(zg(gridj))
       
       xg = oxg
       
       do j = 1,gridj/2
              zg(j) = -ozg(ogk-j+1)
       enddo
       zg(gridj/2+1) = ozg(2)
       
       zg(gridj/2+2:gridj) = ozg(3:ogk)
       
       do j = 1,gridj/2
              Pg(:,j) = oPg(:,yslice,ogk-j+1)
       enddo
       Pg(:,gridj/2+1) = oPg(:,yslice,2)
       
       Pg(:,gridj/2+2:gridj) = oPg(:,yslice,3:ogk)
       
       do j = 1,gridj/2
              Bxg(:,j) = -oBxg(:,yslice,ogk-j+1)
       enddo
       Bxg(:,gridj/2+1) = oBxg(:,yslice,2)
       
       Bxg(:,gridj/2+2:gridj) = oBxg(:,yslice,3:ogk)
       
       do j = 1,gridj/2
              Bzg(:,j) = oBzg(:,yslice,ogk-j+1)
       enddo
       Bzg(:,gridj/2+1) = oBzg(:,yslice,2)
       
       Bzg(:,gridj/2+2:gridj) = oBzg(:,yslice,3:ogk)
       
       xminl = xg(1)
       xmaxl = xg(gridi)
       zminl = zg(1)
       zmaxl = zg(gridj)
       
       write(*,*) xminl,xmaxl
       write(*,*) zminl,zmaxl
       
       !write(*,*) gridj/2,gridj/2+1,gridj/2+2
       !write(*,*) xg(1),xg(gridi)
       !write(*,*) zg(1),zg(gridj/2),zg(gridj/2+1),zg(gridj/2+2),zg(gridj)
       !write(*,*) 'done loading, stopping'
       !write(*,*) Bzg(gridi,1),Bzg(gridi,gridj/2),Bzg(gridi,gridj/2+1),Bzg(gridi,gridj/2+2),Bzg(gridi,gridj)
       !write(*,*) Bxg(gridi,1),Bxg(gridi,gridj/2),Bxg(gridi,gridj/2+1),Bxg(gridi,gridj/2+2),Bxg(gridi,gridj)
       !write(*,*) Pg(gridi,1),Pg(gridi,gridj/2),Pg(gridi,gridj/2+1),Pg(gridi,gridj/2+2),Pg(gridi,gridj)
       !stop
       
       allocate(Ptotg(gridi,gridj))
       
       Ptotg = (Bxg**2+Bzg**2)/(2.0*mu0) + Pg
       
       
       deallocate(oPg)
       deallocate(oBxg)
       deallocate(oByg)
       deallocate(oBzg)
       deallocate(oxg)
       deallocate(oyg)
       deallocate(ozg)
end subroutine readTecFileRCME2d

subroutine modelGrid
       use constants
       use potential
       use parametersBackground
       implicit none
       integer :: i,j,k
       real :: Bx,Bz,pBack,Aback

       !switch off using grid while analytic model is calculated
       useGrid = .false.

       xminl = -40.0
       xmaxl = 0.0
       zminl = -4.0
       zmaxl = 4.0
       gridi = 400
       gridj = 80

       allocate(Pg(gridi,gridj))
       allocate(Ag(gridi,gridj))
       allocate(Bxg(gridi,gridj))
       allocate(Bzg(gridi,gridj))
       allocate(xg(gridi))
       allocate(zg(gridj))

       do i = 1,gridi
              xg(i) = xminl+(xmaxl-xminl)*(i-1)/(gridi-1)
       enddo
       do j = 1,gridj
              zg(j) = zminl+(zmaxl-zminl)*(j-1)/(gridj-1)
       enddo
       do i = 1,gridi
        do j = 1,gridj
              Bxg(i,j) = Bx(xg(i),zg(j),-1)
              Bzg(i,j) = Bz(xg(i),zg(j),-1)
              Pg(i,j) = pBack(xg(i),zg(j),-1)
              Ag(i,j) = Aback(xg(i),zg(j),-1)
        enddo
       enddo
       allocate(Ptotg(gridi,gridj))
       Ptotg = (Bxg**2+Bzg**2)/(2.0*mu0) + Pg

       !back to using the grid
       useGrid = .true.
end subroutine modelGrid

subroutine readAM03Fric
       use constants
       use potential
       use parametersBackground
       use dnams
       implicit none
       integer :: i,j,k,io
       real :: c0,c1,c2,c3,c4,c5
       real :: xin,yin,zin,bxin,byin,bzin,pin,rhoin,vxin,vyin,vzin,fcpin
       character(len=100) :: instr
       character(len=100) :: out
       integer :: ogi,ogj,ogk,yslice
       real :: tslice
       real :: dipx,dipy,dipz
       logical :: useGaussSmoother = .true.
       
       tslice = am03ModelTime
       
       write(*,*) 'reading file: '//trim(am03Model)
       
       xminl = 0.0
       xmaxl = 0.0
       yminl = 0.0
       ymaxl = 0.0
       zminl = 0.0
       zmaxl = 0.0
       
       open(unit=11, file= trim(am03ModelDir)//trim(am03Model), status="old")
       do i = 1,2
              read(11,*) out
              !write(*,*) out
       enddo
       
       do       
              !VARIABLES="I" "J" "X" "Z" "Bx" "Bz" "P"  
              read(11,*,IOSTAT=io) i,j,xin,zin,bxin,bzin,pin
              if(io>0)then
                     write(*,*) 'error in reading tec file'
                     exit
              elseif(io<0)then
                     exit
              else
                     gridi = max(i,gridi)
                     gridj = max(j,gridj)
                     xminl = min(c0,xminl)
                     xmaxl = max(c0,xmaxl)
                     zminl = min(c1,zminl)
                     zmaxl = max(c1,zmaxl)
              endif
       enddo
       close(11)
       
       open(unit=11, file= trim(am03ModelDir)//trim(am03Model), status="old")
       do i = 1,8
              read(11,*) out
       enddo
       allocate(Pg(gridi,gridj))
       allocate(Bxg(gridi,gridj))
       allocate(Bzg(gridi,gridj))
       allocate(xg(gridi))
       allocate(zg(gridj))
       do
              
              read(11,*,IOSTAT=io) i,j,xin,zin,bxin,bzin,pin
              
              if(io>0)then
                     write(*,*) 'error in reading tec file'
                     exit
              elseif(io<0)then
                     exit
              else
                     xg(i) = xin
                     zg(j) = zin
                     Bxg(i,j) = bxin
                     Bzg(i,j) = bzin
                     Pg(i,j) = pin
                     
              endif
       enddo
       close(11)
       
       xminl = xg(1)
       xmaxl = xg(gridi)
       zminl = zg(1)
       zmaxl = zg(gridj)
       
       
       if(useGaussSmoother)then
              !call gaussSmoother(Bxg,xg,zg,gridi,gridj,Bxg)
              !call gaussSmoother(Bzg,xg,zg,gridi,gridj,Bzg)
              !call gaussSmoother(Pg,xg,zg,gridi,gridj,Pg)
              call gaussSmoother2(Bxg,xg,zg,gridi,gridj)
              call gaussSmoother2(Bzg,xg,zg,gridi,gridj)
              call gaussSmoother2(Pg,xg,zg,gridi,gridj)
       endif
       
       if(.true.)then
       c0 = 0.55
       c1 = 0.15
       Pg = c0**2*Pg + c1
       Bxg = c0*Bxg
       Bzg = c0*Bzg
       endif
       
       
       allocate(Ptotg(gridi,gridj))
       Ptotg = (Bxg**2+Bzg**2)/(2.0*mu0) + Pg
       
       write(*,*) 'starting with all_data.dat'
       
       open(unit=11, file= trim(am03ModelDir)//"all_data.dat", status="old")
       i = 1
       do       
              read(11,'(a)',IOSTAT=io) instr
              if(io>0)then
                     write(*,*) 'error in reading file'
                     exit
              elseif(io<0)then
                     exit
              elseif(scan(instr,"h")>0 .or. scan(instr,"r")>0)then
                     cycle
              else
                     read(instr,*) c0,xin,yin,zin,pin,c1,bxin,byin,bzin
              endif
              if(c0 > tslice-0.01 .and. c0 < tslice+0.01)then
                     i = i+1
              endif
              if(c0>tslice+0.01)then
                     exit
              endif
              
       enddo
       close(11)
       
       

       plasGridi = i-1
       
       allocate(plasGrid(plasGridi))
       allocate(plasGridx(plasGridi))
       allocate(plasGridz(plasGridi))
       
       open(unit=11, file= trim(am03ModelDir)//"all_data.dat", status="old")
       i = 1
       do       
              read(11,'(a)',IOSTAT=io) instr
              if(io>0)then
                     write(*,*) 'error in reading file'
                     exit
              elseif(io<0)then
                     exit
              elseif(scan(instr,"h")>0 .or. scan(instr,"r")>0)then
                     cycle
              else
                     read(instr,*) c0,xin,yin,zin,pin,c1,bxin,byin,bzin
              endif
              if(c0 > tslice-0.01 .and. c0 < tslice+0.01)then
                     plasGrid(i) = pin
                     plasGridx(i) = xin
                     plasGridz(i) = zin
                     i = i+1
              endif
              if(c0>tslice+0.01)then
                     exit
              endif
       enddo
       close(11)
       
       write(*,*) 'done with all_data.dat'
       
       !write(*,*) gridj/2,gridj/2+1,gridj/2+2
!       write(*,*) xg(1),xg(gridi)
!       write(*,*) zg(1),zg(gridj/2),zg(gridj/2+1),zg(gridj/2+2),zg(gridj)
!       write(*,*) Bzg(gridi,1),Bzg(gridi,gridj/2),Bzg(gridi,gridj/2+1),Bzg(gridi,gridj/2+2),Bzg(gridi,gridj)
!       write(*,*) Bxg(gridi,1),Bxg(gridi,gridj/2),Bxg(gridi,gridj/2+1),Bxg(gridi,gridj/2+2),Bxg(gridi,gridj)
!       write(*,*) Pg(gridi,1),Pg(gridi,gridj/2),Pg(gridi,gridj/2+1),Pg(gridi,gridj/2+2),Pg(gridi,gridj)
!       stop
end subroutine readAM03Fric

subroutine gaussSmoother(a,x,z,gi,gj,c)
! Kernel smoothing routine for eliminating grid scale oscillations.  AMS 2015
! For a symmetrical kernel of two arguments $K(x_1,x_2) = K(x_2,x_1)$ the smoothing routine calculates
! a smoothed interpolant $\hat y$ of a data set $y$ at a point $x_o$
! $\hat y(x_o) = \frac{\sum\limits_{i=1}^N K(x_o,x_i)y(x_i) }{\sum\limits_{i=1}^N K(x_o,x_i)}$
! This procedure may be commutatively applied to each dimension of an N-dimensional data set
! a=original grid, x=xaxis grid, z= z axis grid, gi=x grid size, gj = z grid size, c=smoothed output
!
       implicit none
       integer,intent(in) :: gi,gj
       real,dimension(gi,gj),intent(in) :: a
       real,dimension(gi),intent(in) :: x
       real,dimension(gj),intent(in) :: z
       real,dimension(gi,gj) :: b,c
       integer :: i,j,k
       real :: c0,c1
       real :: kern !kernel function
       real :: w !smoothing weight
       
       !w = 0.25 !0.25 smoothed out grid scale oscillations and preserved force balance
       w = (x(2)-x(1)) !scales with different grid resolutions, set to max grid spacing
       
       !smooth over x-direction
       do i = 1,gi
         do j = 1,gj
                c0 = 0.0; c1 = 0.0;
                !perform numerator sum
              do k = 1,gi
                     c0 = c0 + a(k,j)*kern(x(i),x(k),w)
              enddo
              !perform denominator sum
              do k = 1,gi
                     c1 = c1 + kern(x(i),x(k),w)
              enddo
              b(i,j) = c0/c1
         enddo
       enddo
       !smooth over z-direction
       do i = 1,gi
         do j = 1,gj
              c0 = 0.0; c1 = 0.0;
              !perform numerator sum
              do k = 1,gj
                     c0 = c0 + b(i,k)*kern(z(j),z(k),w)
              enddo
              !perform denominator sum
              do k = 1,gj
                     c1 = c1 + kern(z(j),z(k),w)
              enddo
              c(i,j) = c0/c1
         enddo
       enddo
       
end subroutine

subroutine gaussSmoother2(a,x,z,gi,gj)
! For a symmetrical kernel of two arguments $K(x_1,x_2) = K(x_2,x_1)$ the smoothing routine calculates
! a smoothed interpolant $\hat y$ of a data set $y$ at a point $x_o$
! $\hat y(x_o) = \frac{\sum\limits_{i=1}^N K(x_o,x_i)y(x_i) }{\sum\limits_{i=1}^N K(x_o,x_i)}$
! This procedure may be commutatively applied to each dimension of an N-dimensional data set
!
! a=original grid, x=xaxis grid, z= z axis grid, gi=x grid size, gj = z grid size, c=smoothed output
!
     implicit none
     integer,intent(in) :: gi,gj
     real,dimension(gi,gj),intent(inout) :: a
     real,dimension(gi),intent(in) :: x
     real,dimension(gj),intent(in) :: z
     real,dimension(gi,gj) :: b,c
     integer :: i,j,k
     real :: c0,c1
     real :: kern !kernel function
     real :: w !smoothing weight
     real :: kf ! kernal function result
     real, parameter :: sigma = 3.0 ! how many sigmas to do smoothing
     
     !w = 0.25 !0.25 smoothed out grid scale oscillations and preserved force balance
     w = (x(2)-x(1)) !scales with different grid resolutions, set to max grid spacing
     
     !smooth over x-direction
     do i = 1,gi
       do j = 1,gj
         c0 = 0.0; c1 = 0.0;
         !perform numerator sum
       do k = 1,gi
         if(abs(x(i)-x(k))<sigma*w)then ! for speed, limit the range to sigma
         kf = kern(x(i),x(k),w)
         c0 = c0 + a(k,j)*kf
       !perform denominator sum
         c1 = c1 + kf
         end if
       enddo
       b(i,j) = c0/c1
       enddo
     enddo
     !smooth over z-direction
     do i = 1,gi
       do j = 1,gj
       c0 = 0.0; c1 = 0.0;
       !perform numerator sum
       do k = 1,gj
         if(abs(z(j)-z(k))<sigma*w)then ! for speed, limit the range to sigma
         kf = kern(z(j),z(k),w)
         c0 = c0 + b(i,k)*kf
       !perform denominator sum
         c1 = c1 + kf
       end if
       enddo
       c(i,j) = c0/c1
       enddo
     enddo
! finally return the new smoothed value
     a = c
     return
end subroutine gaussSmoother2

real function kern(x0,x1,w)
       implicit none
       real,intent(in) :: x0,x1,w
       real,save :: norm = 1/(sqrt(8.0*atan(1.0))) !normalization const.
       !gaussian kernel for smoothing routine
       kern = (norm/w)*exp(-(x0-x1)**2/(2*w**2))
end function

subroutine readAM03
       use constants
       use potential
       use parametersBackground
       use dnams
       implicit none
       integer :: i,j,k,io
       real :: c0,c1,c2,c3,c4,c5
       real :: xin,yin,zin,bxin,byin,bzin,pin,rhoin,vxin,vyin,vzin,fcpin
       character(len=100) :: instr
       integer :: ogi,ogj,ogk,yslice
       real :: dipx,dipy,dipz
       logical :: openPoint
       real :: tslice,upSC
       real :: Bx,Bz
       real,allocatable,dimension(:,:) :: tempArray
       real,allocatable,dimension(:) :: tempArray1
       integer,allocatable,dimension(:,:) :: bcGrid
       integer :: bcGridi
       character(len=20) :: realtostr,filetime

       upSC = upscalingFactor

       write(*,*) 'reading file rd_gse.dat'

       xminl = 0.0
       xmaxl = 0.0
       yminl = 0.0
       ymaxl = 0.0
       zminl = 0.0
       zmaxl = 0.0

       filetime = adjustl(realtostr(am03ModelTime))
       open(unit=11, file= trim(am03ModelDir)//"rd_gse.dat", status="old")
       allocate(oPg(100,100,100))
       allocate(oBxg(100,100,100))
       allocate(oByg(100,100,100))
       allocate(oBzg(100,100,100))
       allocate(oxg(100))
       allocate(oyg(100))
       allocate(ozg(100))
       i = 100
       oxg(i) = 5.0
       j = 1
       do       
              read(11,'(a)',IOSTAT=io) instr
              if(io>0)then
                     write(*,*) 'error in reading tec file'
                     exit
              elseif(io<0)then
                     exit
              elseif(scan(instr,"h")>0 .or. scan(instr,"r")>0)then
                     cycle
              else
                     read(instr,*) c0,xin,yin,zin,bxin,byin,bzin
              endif
              if(xin<=5.0 .and. tslice-0.01<c0 .and. c0<tslice+0.01 )then
                     if(oxg(i)>xin)then
                            oxg(i-1) = xin
                            i = i-1
                     else
                            oxg(i) = xin
                     endif
                     if(j>1)then
                      if(ozg(j-1)>zin)then
                            j = 1
                      endif
                     endif
                     ozg(j) = zin
                     oBxg(i,j,1) = bxin
                     oBzg(i,j,1) = bzin
                     j = j+1
              endif
              if(c0>tslice+0.01)then
                     exit
              endif       
       enddo
       close(11)

       write(*,*) 'done with rd_gse.dat'

       ogi = 100-i
       ogj = j - 1
       gridi = ogi
       gridj = ogj
       allocate(Bxg(gridi,gridj))
       allocate(Byg(gridi,gridj))
       allocate(Bzg(gridi,gridj))
       allocate(xg(gridi))
       allocate(zg(gridj))
       xg = oxg(i:100)
       zg = ozg(1:ogj)
       Bxg = oBxg(i:100,1:ogj,1)
       Bzg = oBzg(i:100,1:ogj,1)
       xminl = xg(1)
       xmaxl = xg(gridi)
       zminl = zg(1)
       zmaxl = zg(gridj)

       !********** bug hunt ************
       allocate(tempArray(gridi,gridj))
       do i = 1,gridi
        do j = 1,gridj
              tempArray(i,j) = Bxg(i,j)
        enddo
       enddo
       deallocate(Bxg)
       allocate(Bxg(gridi,gridj))
       Bxg = tempArray
       deallocate(tempArray)
       allocate(tempArray(gridi,gridj))
       do i = 1,gridi
        do j = 1,gridj
              tempArray(i,j) = Bzg(i,j)
        enddo
       enddo
       deallocate(Bzg)
       allocate(Bzg(gridi,gridj))
       Bzg = tempArray
       deallocate(tempArray)
       allocate(tempArray1(gridi))
       do i = 1,gridi
              tempArray1(i) = xg(i)
       enddo
       deallocate(xg)
       allocate(xg(gridi))
       xg = tempArray1
       deallocate(tempArray1)
       allocate(tempArray1(gridj))
       do j = 1,gridj
              tempArray1(j) = zg(j)
       enddo
       deallocate(zg)
       allocate(zg(gridj))
       zg = tempArray1
       deallocate(tempArray1)
       !***************************

       xminl = xg(1)
       xmaxl = xg(gridi)
       zminl = zg(1)
       zmaxl = zg(gridj)
       write(*,*) xminl,xmaxl,zminl,zmaxl
       allocate(Bxg_d(gridi,gridj,3))
       allocate(Bzg_d(gridi,gridj,3))
       call cdiff(Bxg,Bxg_d,xg,zg,gridi,gridj)
       call cdiff(Bzg,Bzg_d,xg,zg,gridi,gridj)

       if(upSC>1.0)then
              deallocate(oxg)
              deallocate(ozg)
              deallocate(oBxg)
              deallocate(oBzg)
              allocate(oxg(ceiling(upSC*gridi)))
              allocate(ozg(ceiling(upSC*gridj)))
              allocate(oBxg(ceiling(upSC*gridi),ceiling(upSC*gridj),1))
              allocate(oBzg(ceiling(upSC*gridi),ceiling(upSC*gridj),1))
              do i = 1,ceiling(upSC*gridi)
                     oxg(i) = xminl+(xmaxl-xminl)*(i-1)/(ceiling(upSC*gridi)-1)
              enddo
              do j = 1,ceiling(upSC*gridj)
                     ozg(j) = zminl+(zmaxl-zminl)*(j-1)/(ceiling(upSC*gridj)-1)
              enddo
              do i = 1,ceiling(upSC*gridi)
                do j = 1,ceiling(upSC*gridj)
                     oBxg(i,j,1) = Bx(oxg(i),ozg(j),-1)
                     oBzg(i,j,1) = Bz(oxg(i),ozg(j),-1)
                enddo
              enddo
              gridi = ceiling(upSC*gridi)
              gridj = ceiling(upSC*gridj)
              deallocate(xg)
              deallocate(zg)
              allocate(xg(gridi))
              allocate(zg(gridj))
              xg = oxg
              zg = ozg
              deallocate(Bxg)
              deallocate(Bzg)
              allocate(Bxg(gridi,gridj))
              allocate(Bzg(gridi,gridj))
              Bxg = oBxg(:,:,1)
              Bzg = oBzg(:,:,1)
              deallocate(Bxg_d)
              deallocate(Bzg_d)
              allocate(Bxg_d(gridi,gridj,3))
              allocate(Bzg_d(gridi,gridj,3))
              call cdiff(Bxg,Bxg_d,xg,zg,gridi,gridj)
              call cdiff(Bzg,Bzg_d,xg,zg,gridi,gridj)
       endif

       allocate(Pg(gridi,gridj))
       write(*,*) 'starting with all_data.dat'
       open(unit=11, file= trim(am03ModelDir)//"all_data.dat", status="old")
       i = 1
       do       
              read(11,'(a)',IOSTAT=io) instr
              if(io>0)then
                     write(*,*) 'error in reading file'
                     exit
              elseif(io<0)then
                     exit
              elseif(scan(instr,"h")>0 .or. scan(instr,"r")>0)then
                     cycle
              else
                     read(instr,*) c0,xin,yin,zin,pin,c1,bxin,byin,bzin
              endif
              if(c0 > tslice-0.01 .and. c0 < tslice+0.01)then
                     i = i+1
              endif
              if(c0>tslice+0.01)then
                     exit
              endif
              
       enddo
       close(11)

       plasGridi = i-1
       allocate(plasGrid(plasGridi))
       allocate(plasGridx(plasGridi))
       allocate(plasGridz(plasGridi))
       open(unit=11, file= trim(am03ModelDir)//"all_data.dat", status="old")
       i = 1
       do       
              read(11,'(a)',IOSTAT=io) instr
              if(io>0)then
                     write(*,*) 'error in reading file'
                     exit
              elseif(io<0)then
                     exit
              elseif(scan(instr,"h")>0 .or. scan(instr,"r")>0)then
                     cycle
              else
                     read(instr,*) c0,xin,yin,zin,pin,c1,bxin,byin,bzin
              endif
              if(c0 > tslice-0.01 .and. c0 < tslice+0.01)then
                     plasGrid(i) = pin
                     plasGridx(i) = xin
                     plasGridz(i) = zin
                     i = i+1
              endif
              if(c0>tslice+0.01)then
                     exit
              endif
       enddo
       close(11)
       write(*,*) 'done with all_data.dat'

       open(unit=11, file= trim(am03ModelDir)//"radialPres.dat", status="replace")
       do i = 1,plasGridi
              write(11,*) plasGridx(i),plasGridz(i),plasGrid(i)
       enddo
       close(11)
       open(unit=11, file= trim(am03ModelDir)//"radialBz.dat", status="replace")
       do i = 1,gridi
       j=43
         write(11,*) xg(i),zg(j),Bxg(i,j),Bzg(i,j),Bxg_d(i,j,1),Bzg_d(i,j,1),Bxg_d(i,j,2),Bzg_d(i,j,2),Bxg_d(i,j,3),Bzg_d(i,j,3)
       enddo
       close(11)
       open(unit=11, file= trim(am03ModelDir)//"radialBbcint.dat", status="replace")
       do i = 1,300
              c0 = -5.0-29.0*(i-1)/(300-1)
              c1 = 0.0
              write(11,*) c0,Bx(c0,c1,-1),Bz(c0,c1,-1)
       enddo
       close(11)
       open(unit=11, file= trim(am03ModelDir)//"Binterp_debug.dat", status="replace")
       write(11,*) 'TITLE="debug"'
       write(11,*) 'VARIABLES="I" "J" "X" "Z" "Bx" "Bz"'
       !write(11,*) 'ZONE T="IT=00003410"'
       do i = 1,300
        do j = 1,300
              c0 = -5.0-29.0*(i-1)/(300-1)
              c1 = -10.0+20.0*(j-1)/(300-1)
              write(11,*) i,j,c0,c1,log(abs(Bx(c0,c1,-1))),log(abs(Bz(c0,c1,-1)))
        enddo
       enddo       
       close(11)
       open(unit=11, file= trim(am03ModelDir)//"Binterp_debug2.dat", status="replace")
       write(11,*) 'TITLE="debug"'
       write(11,*) 'VARIABLES="I" "J" "X" "Z" "Bx" "Bz" "Bxdx" "Bzdx" "Bxdz" "Bzdz" "Bxdxdz" "Bzdxdz"'
       do i = 1,gridi
        do j = 1,gridj
          write(11,*) i,j,xg(i),zg(j),Bxg(i,j),Bzg(i,j),Bxg_d(i,j,1),&
                  Bzg_d(i,j,1),Bxg_d(i,j,2),Bzg_d(i,j,2),Bxg_d(i,j,3),Bzg_d(i,j,3)
        enddo
       enddo
       close(11)
       allocate(testMesh(gridi,gridj))
       allocate(testMesh_d(gridi,gridj,3))
       do i = 1,gridi
         do j = 1,gridj
              testMesh(i,j) = sqrt(xg(i)**2+zg(j)**2)
         enddo
       enddo
       call cdiff(testMesh,testMesh_d,xg,zg,gridi,gridj)
       open(unit=11, file= trim(am03ModelDir)//"Binterp_debug3.dat", status="replace")
       write(11,*) 'TITLE="debug"'
       write(11,*) 'VARIABLES="I" "J" "X" "Z" "r" "rdx" "rdz" "rdxdz"'
       do i = 1,gridi
        do j = 1,gridj
              write(11,*) i,j,xg(i),zg(j),testMesh(i,j),testMesh_d(i,j,1),testMesh_d(i,j,2),testMesh_d(i,j,3)
        enddo
       enddo
       close(11)
       bcGridi = 300
       allocate(bcGrid(bcGridi,8))
       do i = 1,bcGridi
              c0 = plasGridx(1)+(plasGridx(plasGridi)-plasGridx(1))*(i-1)/(bcGridi-1)
              do j = 1,plasGridi-1
                     if (plasGridx(j)<=c0 .and. c0<=plasGridx(j+1)) then
                            exit
                     endif
              enddo
              c1 = plasGridz(j)+(plasGridz(j+1)-plasGridz(j))*(c0-plasGridx(j))/(plasGridx(j+1)-plasGridx(j))
              call bcNbPoints(c0,c1,bcGrid(i,:))
       enddo
       open(unit=11, file= trim(am03ModelDir)//"marina_bc.dat", status="replace")
       write(11,*) bcGridi
       do i = 1,bcGridi
              write(11,*) i,bcGrid(i,:)
       enddo
       close(11)
       write(*,*) 'calculating pressure from field line tracing'
       
       do i = 1,gridi
         do j = 1,gridj
              if(i==1 .or. j==1 .or. i==gridi .or. j==gridj)then
                     Pg(i,j) = plasGrid(1)
              else
                     call tracerRK4pres(i,j,xg(i),zg(j),plasGrid,plasGridx,plasGridz,Pg(i,j),openPoint)
                     !if (.not. openPoint) then
                            !write(*,*) 'point not open,  P=',Pg(i,j),'i=',i,'j=',j,xg(i),zg(j)
                     !else
                            !write(*,*) 'point was set to P= ',plasGrid(1),'  i=',i,'j=',j,xg(i),zg(j),"*******************"
                     !endif
              endif
         enddo
         write(*,*) (100.0*i)/gridi,'% done'
       enddo
       
       open(unit=11, file= trim(am03ModelDir)//"marina_"//filetime(1:5)//".dat", status="replace")
       write(11,*) 'TITLE="marina data"'
       write(11,*) 'VARIABLES="I" "J" "X" "Z" "Bx" "Bz" "P"'
       !write(11,*) 'ZONE T="IT=00003410"'
       do i = 1,gridi
        do j = 1,gridj
              write(11,*) i,j,xg(i),zg(j),Bxg(i,j),Bzg(i,j),Pg(i,j)
        enddo
       enddo
       close(11)
       allocate(Ptotg(gridi,gridj))
       Ptotg = (Bxg**2+Bzg**2)/(2.0*mu0) + Pg
       deallocate(oPg)
       deallocate(oBxg)
       deallocate(oByg)
       deallocate(oBzg)
       deallocate(oxg)
       deallocate(oyg)
       deallocate(ozg)
       deallocate(Bxg_d)
       deallocate(Bzg_d)
       deallocate(testMesh)
       deallocate(testMesh_d)
end subroutine readAM03

subroutine bcNbPoints(x,z,bc)
       use potential
       implicit none
       integer,dimension(8) :: bc
       real :: x,z
       integer :: i,j
       
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
       bc(1) = i
       bc(2) = j
       bc(3) = i+1
       bc(4) = j
       bc(5) = i+1
       bc(6) = j+1
       bc(7) = i
       bc(8) = j+1
end subroutine

subroutine setupDensGrid
       use constants
       use potential
       use parametersBackground
       use initCond
       implicit none
       integer :: i,j
       real :: Nps_tm2003_new
       allocate(Ng(gridi,gridj))
       do i = 1,gridi
        do j = 1,gridj
              if(use_gallagherDensity)then
                     if(sqrt(xg(i)**2+zg(j)**2)<boundary)then
                            dBack = 1.0
                     else
                            call gallagher(xg(i),0.0,zg(j),dBack)
                     endif
              endif
              if(use_tm2003_density)then
                     dBack = Nps_tm2003_new(abs(xg(i)),0.0,sw_density,sw_velocity,bzimf)
              endif
              
              Ng(i,j) = dBack
        enddo
       enddo
end subroutine


subroutine calculateBG_grads
       use constants
       use potential
       use parametersBackground
       implicit none
       integer :: i,j
       
       allocate(dx_Bxg(gridi,gridj))
       allocate(dx_Bzg(gridi,gridj))
       allocate(dx_Ptotg(gridi,gridj))
       allocate(dx_Pg(gridi,gridj))
       allocate(dz_Bxg(gridi,gridj))
       allocate(dz_Bzg(gridi,gridj))
       allocate(dz_Ptotg(gridi,gridj))
       allocate(dz_Pg(gridi,gridj))
       
       allocate(Ag_d(gridi,gridj,3))
       allocate(Bxg_d(gridi,gridj,3))
       allocate(Bzg_d(gridi,gridj,3))
       allocate(Pg_d(gridi,gridj,3))
       allocate(Ng_d(gridi,gridj,3))
       allocate(Ptotg_d(gridi,gridj,3))
       allocate(testMesh_d(gridi,gridj,3))
       
       call cdiff(Ag,Ag_d,xg,zg,gridi,gridj)
       call cdiff(Bxg,Bxg_d,xg,zg,gridi,gridj)
       call cdiff(Bzg,Bzg_d,xg,zg,gridi,gridj)
       call cdiff(Pg,Pg_d,xg,zg,gridi,gridj)
       call cdiff(Ng,Ng_d,xg,zg,gridi,gridj)
       call cdiff(Ptotg,Ptotg_d,xg,zg,gridi,gridj)
       call cdiff(testMesh,testMesh_d,xg,zg,gridi,gridj)
       
       do i = 2,gridi-1
         do j = 2,gridj-1
              dx_Bxg(i,j) = (Bxg(i+1,j)-Bxg(i-1,j))/(xg(i+1)-xg(i-1))
              dx_Bzg(i,j) = (Bzg(i+1,j)-Bzg(i-1,j))/(xg(i+1)-xg(i-1))
              dx_Ptotg(i,j) = (Ptotg(i+1,j)-Ptotg(i-1,j))/(xg(i+1)-xg(i-1))
              dx_Pg(i,j) = (Pg(i+1,j)-Pg(i-1,j))/(xg(i+1)-xg(i-1))
              dz_Bxg(i,j) = (Bxg(i,j+1)-Bxg(i,j-1))/(zg(j+1)-zg(j-1))
              dz_Bzg(i,j) = (Bzg(i,j+1)-Bzg(i,j-1))/(zg(j+1)-zg(j-1))
              dz_Ptotg(i,j) = (Ptotg(i,j+1)-Ptotg(i,j-1))/(zg(j+1)-zg(j-1))
              dz_Pg(i,j) = (Pg(i,j+1)-Pg(i,j-1))/(zg(j+1)-zg(j-1))
         enddo
       enddo
       !******************************************************************
       i = 1
         do j = 1,gridj
              dx_Bxg(i,j) = (Bxg(i+1,j)-Bxg(i,j))/(xg(i+1)-xg(i))
              dx_Bzg(i,j) = (Bzg(i+1,j)-Bzg(i,j))/(xg(i+1)-xg(i))
              dx_Ptotg(i,j) = (Ptotg(i+1,j)-Ptotg(i,j))/(xg(i+1)-xg(i))
              dx_Pg(i,j) = (Pg(i+1,j)-Pg(i,j))/(xg(i+1)-xg(i))
         enddo
       
       i = gridi
         do j = 1,gridj
              dx_Bxg(i,j) = (Bxg(i,j)-Bxg(i-1,j))/(xg(i)-xg(i-1))
              dx_Bzg(i,j) = (Bzg(i,j)-Bzg(i-1,j))/(xg(i)-xg(i-1))
              dx_Ptotg(i,j) = (Ptotg(i,j)-Ptotg(i-1,j))/(xg(i)-xg(i-1))
              dx_Pg(i,j) = (Pg(i,j)-Pg(i-1,j))/(xg(i)-xg(i-1))
         enddo
         
       j = 1
         do i = 2,gridi-1
              dx_Bxg(i,j) = (Bxg(i+1,j)-Bxg(i-1,j))/(xg(i+1)-xg(i-1))
              dx_Bzg(i,j) = (Bzg(i+1,j)-Bzg(i-1,j))/(xg(i+1)-xg(i-1))
              dx_Ptotg(i,j) = (Ptotg(i+1,j)-Ptotg(i-1,j))/(xg(i+1)-xg(i-1))
              dx_Pg(i,j) = (Pg(i+1,j)-Pg(i-1,j))/(xg(i+1)-xg(i-1))
         enddo
       
       j = gridj
         do i = 2,gridi-1
              dx_Bxg(i,j) = (Bxg(i+1,j)-Bxg(i-1,j))/(xg(i+1)-xg(i-1))
              dx_Bzg(i,j) = (Bzg(i+1,j)-Bzg(i-1,j))/(xg(i+1)-xg(i-1))
              dx_Ptotg(i,j) = (Ptotg(i+1,j)-Ptotg(i-1,j))/(xg(i+1)-xg(i-1))
              dx_Pg(i,j) = (Pg(i+1,j)-Pg(i-1,j))/(xg(i+1)-xg(i-1))
         enddo
       
       !******************************************************************
       j = 1
         do i = 1,gridi
              dz_Bxg(i,j) = (Bxg(i,j+1)-Bxg(i,j))/(zg(j+1)-zg(j))
              dz_Bzg(i,j) = (Bzg(i,j+1)-Bzg(i,j))/(zg(j+1)-zg(j))
              dz_Ptotg(i,j) = (Ptotg(i,j+1)-Ptotg(i,j))/(zg(j+1)-zg(j))
              dz_Pg(i,j) = (Pg(i,j+1)-Pg(i,j))/(zg(j+1)-zg(j))
         enddo
       
       j = gridj
         do i = 1,gridi
              dz_Bxg(i,j) = (Bxg(i,j)-Bxg(i,j-1))/(zg(j)-zg(j-1))
              dz_Bzg(i,j) = (Bzg(i,j)-Bzg(i,j-1))/(zg(j)-zg(j-1))
              dz_Ptotg(i,j) = (Ptotg(i,j)-Ptotg(i,j-1))/(zg(j)-zg(j-1))
              dz_Pg(i,j) = (Pg(i,j)-Pg(i,j-1))/(zg(j)-zg(j-1))
         enddo
       
       i = 1
         do j = 2,gridj-1
              dz_Bxg(i,j) = (Bxg(i,j+1)-Bxg(i,j-1))/(zg(j+1)-zg(j-1))
              dz_Bzg(i,j) = (Bzg(i,j+1)-Bzg(i,j-1))/(zg(j+1)-zg(j-1))
              dz_Ptotg(i,j) = (Ptotg(i,j+1)-Ptotg(i,j-1))/(zg(j+1)-zg(j-1))
              dz_Pg(i,j) = (Pg(i,j+1)-Pg(i,j-1))/(zg(j+1)-zg(j-1))
         enddo

         i = gridi
         do j = 2,gridj-1
              dz_Bxg(i,j) = (Bxg(i,j+1)-Bxg(i,j-1))/(zg(j+1)-zg(j-1))
              dz_Bzg(i,j) = (Bzg(i,j+1)-Bzg(i,j-1))/(zg(j+1)-zg(j-1))
              dz_Ptotg(i,j) = (Ptotg(i,j+1)-Ptotg(i,j-1))/(zg(j+1)-zg(j-1))
              dz_Ptotg(i,j) = (Pg(i,j+1)-Pg(i,j-1))/(zg(j+1)-zg(j-1))
         enddo
end subroutine calculateBG_grads

subroutine cdiff(a,b,x,z,gridi,gridj)
       implicit none
       real,dimension(gridi,gridj) :: a
       real,dimension(gridi,gridj,3) :: b
       real,dimension(gridi) :: x
       real,dimension(gridj) :: z
       integer :: i,j,gridi,gridj
       b = 0.0
       do i = 2,gridi-1       
        do j = 2,gridj-1
              b(i,j,1) = (a(i+1,j)-a(i-1,j))/(x(i+1)-x(i-1))
              b(i,j,2) = (a(i,j+1)-a(i,j-1))/(z(j+1)-z(j-1))
              b(i,j,3) = (a(i+1,j+1)-a(i-1,j+1)-a(i+1,j-1)+a(i-1,j-1))/((x(i+1)-x(i-1))*(z(j+1)-z(j-1)))
        enddo
       enddo
end subroutine cdiff

subroutine readMatlabFile(fileName)
       use constants
       use potential
       use parametersBackground, only: A0
       implicit none
       integer :: i,j,k,io
       real :: c0,c1,c2,c3,c4,c5
       character(len=20) :: fileName
       character(len=20) :: out
       
       write(*,*) 'reading file:'
       write(*,*) fileName
       
       xminl = 0.0
       xmaxl = 0.0
       zminl = 0.0
       zmaxl = 0.0
       
       open(unit=11, file= fileName, status="old")
       do i = 1,8
              read(11,*) out
              !write(*,*) out
       enddo
       
       do       
              read(11,*,IOSTAT=io) i,k,c0,c1,c2,c3
              if(io>0)then
                     write(*,*) 'error in reading tec file'
                     exit
              elseif(io<0)then
                     exit
              else
                     gridi = max(i,gridi)
                     gridj = max(k,gridj)
                     xminl = min(c0,xminl)
                     xmaxl = max(c0,xmaxl)
                     zminl = min(c1,zminl)
                     zmaxl = max(c1,zmaxl)
              endif
       enddo
       close(11)
       
       write(*,*) 'xrange:',xminl,xmaxl
       write(*,*) 'zrange:',zminl,zmaxl
       
       open(unit=11, file= fileName, status="old")
       do i = 1,8
              read(11,*) out
       enddo
       allocate(Pg(gridi,gridj))
       allocate(Ag(gridi,gridj))
       allocate(xg(gridi))
       allocate(zg(gridj))
       do       
              read(11,*,IOSTAT=io) i,k,c0,c1,c2,c3
              
              c2 = c2*A0*A0/mu0
              c3 = -c3*A0
              
              if(io>0)then
                     write(*,*) 'error in reading tec file'
                     exit
              elseif(io<0)then
                     exit
              else
                     xg(i) = c0
                     zg(k) = c1
                     Pg(i,k) = c2
                     Ag(i,k) = c3
              endif
       enddo
       close(11)
end subroutine readMatlabFile

subroutine setupGrid
       use constants
       use potential
       use initCond
       use parametersBackground
       implicit none
       integer :: i,j
       real :: Aback,pBack
       xminl = 0.0
       xmaxl = 0.0
       zminl = 0.0
       zmaxl = 0.0
       
       gridi = 100
       gridj = 100

       allocate(Pg(gridi,gridj))
       allocate(Ag(gridi,gridj))
       allocate(xg(gridi))
       allocate(zg(gridj))
       
       do i = 1,gridi
              xg(i) = boundary-1.0 + (1.0+apex)*(i-1)/(gridi-1)
       end do
       
       do j = 1,gridj
              zg(j) = -ht + 2.0*ht*(j-1)/(gridj-1)
       end do
       
       do i = 1,gridi
         do j = 1,gridj
              Ag(i,j) = Aback(xg(i),zg(j),-1)
         end do
       end do
       
       do i = 1,gridi
         do j = 1,gridj
              Pg(i,j) = pBack(xg(i),zg(j),-1)
         end do
       end do
       
       xminl = xg(1)
       do i = 1,gridi
              xminl = min(xg(i),xminl)
       end do
       zminl = zg(1)
       do j = 1,gridj
              zminl = min(zg(j),zminl)
       end do
       
       xmaxl = xg(1)
       do i = 1,gridi
              xmaxl = max(xg(i),xmaxl)
       end do
       zmaxl = zg(1)
       do j = 1,gridj
              zmaxl = max(zg(j),zmaxl)
       end do
       
end subroutine setupGrid

subroutine addGridToFile(t,a,fileName,gridi,gridj)
       implicit none
       integer :: gridi,gridj
       real,dimension(gridi,gridj) :: a
       integer :: i,j,bufsize
       character(len=20) :: fileName
       real :: t
       inquire(Iolength=bufsize) a
       open(unit=12, file= fileName, position="append", status="old", Recl=bufsize)
         write(12,*) t
         write(12,*) gridi
         write(12,*) gridj
         do j = 1,gridj
          do i = 1,gridi
              if(i==gridi)then
                      write(12,"(f12.4)") a(i,j)
              else
                     write(12,"(f12.4,a)",advance='no') a(i,j),' '
              end if
          end do
         end do
       close(unit=12)
end subroutine

subroutine resetGrid
       use constants
       use potential
       use initCond
       use parametersBackground
       implicit none
       integer :: i,j
       real :: Aback,pBack
       xminl = 0.0
       xmaxl = 0.0
       zminl = 0.0
       zmaxl = 0.0
       
       gridi = 100
       gridj = 100

       do i = 1,gridi
              xg(i) = boundary-1.0 + (1.0+apex)*(i-1)/(gridi-1)
       end do
       do j = 1,gridj
              zg(j) = -ht + 2.0*ht*(j-1)/(gridj-1)
       end do
       do i = 1,gridi
         do j = 1,gridj
              Ag(i,j) = Aback(xg(i),zg(j),-1)
         end do
       end do
       do i = 1,gridi
         do j = 1,gridj
              Pg(i,j) = pBack(xg(i),zg(j),-1)
         end do
       end do
       xminl = xg(1)
       do i = 1,gridi
              xminl = min(xg(i),xminl)
       end do
       zminl = zg(1)
       do j = 1,gridj
              zminl = min(zg(j),zminl)
       end do
       
       xmaxl = xg(1)
       do i = 1,gridi
              xmaxl = max(xg(i),xmaxl)
       end do
       zmaxl = zg(1)
       do j = 1,gridj
              zmaxl = max(zg(j),zmaxl)
       end do
       
end subroutine resetGrid



subroutine readTsyFric
       use constants
       use potential
       use parametersBackground
       use dnams
       implicit none
       integer :: i,j,k,io
       real :: c0,c1,c2,c3,c4,c5
       real :: xin,yin,zin,bxin,byin,bzin,pin,rhoin,vxin,vyin,vzin,fcpin
       character(len=100) :: instr
       character(len=100) :: out
       integer :: ogi,ogj,ogk,yslice
       real :: tslice
       real :: dipx,dipy,dipz,Bx,Bz,pgxmax,pwolf,pgxmin
       logical :: useGaussSmoother = .true.
       
       
       write(*,*) 'reading file: '//trim(tsyModelDir)//"/"//trim(tsyModel)
       
       xminl = 0.0
       xmaxl = 0.0
       yminl = 0.0
       ymaxl = 0.0
       zminl = 0.0
       zmaxl = 0.0
       
       open(unit=11, file= trim(tsyModelDir)//"/"//trim(tsyModel), status="old")
       do i = 1,2
              read(11,*) out
              !write(*,*) out
       enddo
       
       do       
              !VARIABLES="I" "J" "X" "Z" "Bx" "Bz" "P"  
              read(11,*,IOSTAT=io) i,j,xin,zin,bxin,bzin,pin
              if(io>0)then
                     write(*,*) 'error in reading tec file'
                     exit
              elseif(io<0)then
                     exit
              else
                     gridi = max(i,gridi)
                     gridj = max(j,gridj)
                     xminl = min(c0,xminl)
                     xmaxl = max(c0,xmaxl)
                     zminl = min(c1,zminl)
                     zmaxl = max(c1,zmaxl)
              endif
       enddo
       close(11)
       
       open(unit=11, file= trim(tsyModelDir)//"/"//trim(tsyModel), status="old")
       do i = 1,8
              read(11,*) out
       enddo
       allocate(Pg(gridi,gridj))
       allocate(Bxg(gridi,gridj))
       allocate(Bzg(gridi,gridj))
       allocate(xg(gridi))
       allocate(zg(gridj))
       do
              
              read(11,*,IOSTAT=io) i,j,xin,zin,bxin,bzin,pin
              
              if(io>0)then
                     write(*,*) 'error in reading tec file'
                     exit
              elseif(io<0)then
                     exit
              else
                     xg(i) = xin
                     zg(j) = zin
                     Bxg(i,j) = bxin
                     Bzg(i,j) = bzin
                     Pg(i,j) = pin
                     
              endif
       enddo
       close(11)
       
       xminl = xg(1)
       xmaxl = xg(gridi)
       zminl = zg(1)
       zmaxl = zg(gridj)
       
       
       if(useGaussSmoother)then
              !call gaussSmoother(Bxg,xg,zg,gridi,gridj,Bxg)
              !call gaussSmoother(Bzg,xg,zg,gridi,gridj,Bzg)
              !call gaussSmoother(Pg,xg,zg,gridi,gridj,Pg)
              call gaussSmoother2(Bxg,xg,zg,gridi,gridj)
              call gaussSmoother2(Bzg,xg,zg,gridi,gridj)
              call gaussSmoother2(Pg,xg,zg,gridi,gridj)
       endif
       
       
       allocate(Ptotg(gridi,gridj))
       Ptotg = (Bxg**2+Bzg**2)/(2.0*mu0) + Pg

       write(*,*) 'allocating plasGrid'
       allocate(plasGrid(plasGridi))
       allocate(plasGridx(plasGridi))
       allocate(plasGridz(plasGridi))
       pgxmin = txmin+1.0
       pgxmax = -1.0
       do i = 1,plasGridi
              plasGridx(i) = pgxmin+(pgxmax-pgxmin)*(i-1)/(plasGridi-1)
       enddo
       write(*,*) 'finding sheet'

       allocate(Bxg_d(gridi,gridj,3))
       allocate(Bzg_d(gridi,gridj,3))
       call cdiff(Bxg,Bxg_d,xg,zg,gridi,gridj)
       call cdiff(Bzg,Bzg_d,xg,zg,gridi,gridj)
       
       do i = 1,plasGridi
              c1=100.0
              c2=0.0
              do j = 1,10000
                     c0 = zminl + 5.0 + (zmaxl-zminl - 10.0)*(j-1)/(10000-1)
                     if(c1>abs(Bx(plasGridx(i),c0,-1)))then
                            c1 = abs(Bx(plasGridx(i),c0,-1))
                            c5 = c0
                     endif
              enddo
              plasGridz(i) = c5
       enddo

       do i = 1,plasGridi
              plasGrid(i) = &
!              min(pModelWang2013(sqrt(plasGridx(i)**2+plasGridz(i)**2),pi,kpindex,0),pressureCap)
              min(pwolf(sqrt(plasGridx(i)**2+plasGridz(i)**2)),pressureCap)
       enddo
       
       deallocate(Bxg_d)
       deallocate(Bzg_d)
end subroutine readTsyFric


subroutine makeTsyData
       use constants
       use potential
       use parametersBackground
       use dnams
       implicit none
       integer :: i,j,k,l
       real :: c0,c1,c2,c3,c4,c5,dxBz,dzBx,Bx,Bz,pwolf
       real :: dipx,dipy,dipz
       logical :: openPoint
       real :: pgxmax,pgxmin,pgzmax,pgzmin
       real,allocatable,dimension(:) :: plasGridCopy
       real :: smoother1D
       real*4,dimension(10) :: parmod
       real*4 :: bxr4,byr4,bzr4,tiltr4,g1r4,g2r4,byimfr4,bzimfr4,dstr4,swpressr4,xr4,yr4,zr4
       real*4 :: r0 = 0.0

       plasGridi = 400

       xminl = txmin
       xmaxl = txmax
       yminl = tyslice
       ymaxl = tyslice
       zminl = tzmin
       zmaxl = tzmax

       gridi = tgridx
       gridj = tgridz

       allocate(Bxg(gridi,gridj))
       allocate(Byg(gridi,gridj))
       allocate(Bzg(gridi,gridj))
       allocate(xg(gridi))
       allocate(zg(gridj))
       
       do i = 1,gridi
              xg(i) = xminl+(xmaxl-xminl)*(i-1)/(gridi-1)
       enddo
       do j = 1,gridj
              zg(j) = zminl+(zmaxl-zminl)*(j-1)/(gridj-1)
       enddo
       do i = 1,gridi
         do j = 1,gridj
                tiltr4=real(tilt,4)
                g1r4=real(tsyg1,4)
                g2r4=real(tsyg2,4)
                byimfr4=real(byimf,4)
                bzimfr4=real(bzimf,4)
                dstr4=real(dst,4)
                swpressr4=real(swpress,4)
                xr4=real(xg(i),4)
                yr4=real(tyslice,4)
                zr4=real(zg(j),4)
                parmod = (/swpressr4,dstr4,byimfr4,bzimfr4,g1r4,g2r4,r0,r0,r0,r0/)
              select case(tsyVersion)
              case('t96')
                     call T96_01(0,parmod,tiltr4,xr4,yr4,zr4,bxr4,byr4,bzr4)
                     call Dipole_modified_fil(dipole_moment, tilt, xg(i), tyslice, zg(j), dipx, dipy, dipz)
              case('t89')
                     call T89D_SP(kpindex,parmod,tiltr4,xr4,yr4,zr4,bxr4,byr4,bzr4)
                     call Dipole_modified_fil(dipole_moment, tilt, xg(i), tyslice, zg(j), dipx, dipy, dipz)
              case('t01')
                     call T01_01(0,parmod,tiltr4,xr4,yr4,zr4,bxr4,byr4,bzr4)
                     call Dipole_modified_fil(dipole_moment, tilt, xg(i), tyslice, zg(j), dipx, dipy, dipz)
              case default
                     write(*,*) 'Error, no tsy model match'
                     stop
              end select
              Bxg(i,j) = bxr4
              Byg(i,j) = byr4
              Bzg(i,j) = bzr4
              Bxg(i,j) = Bxg(i,j)+dipx
              Byg(i,j) = Byg(i,j)+dipy
              Bzg(i,j) = Bzg(i,j)+dipz
         enddo
       enddo

       allocate(Bxg_d(gridi,gridj,3))
       allocate(Bzg_d(gridi,gridj,3))
       call cdiff(Bxg,Bxg_d,xg,zg,gridi,gridj)
       call cdiff(Bzg,Bzg_d,xg,zg,gridi,gridj)

       allocate(Pg(gridi,gridj))

       write(*,*) 'allocating plasGrid'
       allocate(plasGrid(plasGridi))
       allocate(plasGridx(plasGridi))
       allocate(plasGridz(plasGridi))
       pgxmin = txmin+1.0
       pgxmax = -1.0
       do i = 1,plasGridi
              plasGridx(i) = pgxmin+(pgxmax-pgxmin)*(i-1)/(plasGridi-1)
       enddo
       write(*,*) 'finding sheet'
       
       do i = 1,plasGridi
              c1=100.0
              c2=0.0
              do j = 1,10000
                     c0 = zminl + 5.0 + (zmaxl-zminl - 10.0)*(j-1)/(10000-1)
                     if(c1>abs(Bx(plasGridx(i),c0,-1)))then
                            c1 = abs(Bx(plasGridx(i),c0,-1))
                            c5 = c0
                     endif
              enddo
              plasGridz(i) = c5
       enddo

       j = 0
       do i = 1,plasGridi
              plasGrid(i) = & 
!              min(pressureCap,pModelWang2013(sqrt(plasGridx(i)**2+plasGridz(i)**2),pi,kpindex,0))
              min(pressureCap,pwolf(sqrt(plasGridx(i)**2+plasGridz(i)**2)))
              if (plasGrid(i)==pressureCap .and. j==0) j=i
       enddo

       if (j>0) then
              write(*,*) 'smoothing discontinuity in plasGrid'
              allocate(plasGridCopy(plasGridi))
              plasGridCopy = plasGrid
              k = plasGridi-j
              do i = j-k,plasGridi
                            plasGrid(i) = smoother1D(plasGridCopy,plasGridi,i,j,k)
              enddo
              deallocate(plasGridCopy)
       endif

       ! do i =1,plasGridi
       !        write(*,*) abs(plasGridx(i)), plasGrid(i)
       ! enddo
       ! stop

       write(*,*) 'calculating pressure from field line tracing'
       
       do i = 1,gridi
         do j = 1,gridj
              if(i==1 .or. j==1 .or. i==gridi .or. j==gridj)then
                     Pg(i,j) = plasGrid(1)
              else
                     call tracerRK4pres(i,j,xg(i),zg(j),plasGrid,plasGridx,plasGridz,Pg(i,j),openPoint)
                     !if (.not. openPoint) then
                            !write(*,*) 'point not open,  P=',Pg(i,j),'i=',i,'j=',j,xg(i),zg(j)
                     !else
                            !write(*,*) 'point was set to P= ',plasGrid(1),'  i=',i,'j=',j,xg(i),zg(j),"*******************"
                     !endif
              endif
         enddo
         write(*,*) (100.0*i)/gridi,'% done'
       enddo

       call fric2d_wrapper(Bxg,Bzg,Pg,xg,zg,gridi,gridj,tilt,Bxg,Bzg,Pg,fricDir)

       open(unit=11, file= trim(tsyModelDir)//"/"//trim(tsyModel), status="replace")
       write(11,*) 'TITLE="tsy data"'
       write(11,*) 'VARIABLES="I" "J" "X" "Z" "Bx" "Bz" "P"'
       do i = 1,gridi
        do j = 1,gridj
              write(11,*) i,j,xg(i),zg(j),Bxg(i,j),Bzg(i,j),Pg(i,j)
        enddo
       enddo
       close(11)
       call gaussSmoother2(Bxg,xg,zg,gridi,gridj)
       call gaussSmoother2(Bzg,xg,zg,gridi,gridj)
       call gaussSmoother2(Pg,xg,zg,gridi,gridj)
       allocate(Ptotg(gridi,gridj))
       Ptotg = (Bxg**2+Bzg**2)/(2.0*mu0) + Pg
       deallocate(Bxg_d)
       deallocate(Bzg_d)
end subroutine makeTsyData

real function smoother1D(plasGridCopy,plasGridi,i,j,k)
       implicit none
       integer :: plasGridi,i,j,k,iii,jjj
       real,dimension(plasGridi) :: plasGridCopy
       real :: c0,c1,c2,c3
       do jjj = 1,3
         do iii=j-k+1,plasGridi-1
              plasGridCopy(iii) = 0.5*plasGridCopy(iii+1)+0.5*plasGridCopy(iii-1)
         enddo
       enddo
       smoother1D = plasGridCopy(i)
end function
