
!********************************* initial condition subroutines *********************************************
! These routines are meant as a scratch space to make plots and setup the initial conditions
! using repeated field line traces and approximation methods.
!
! State is set up for its initial configuration at the end after calculations have been
! made. Before that the state, oldstate, and newstate arrays are repeatedly used as dummy arrays.
!*************************************************************************************************************

subroutine initialConditions1

  use stateMod
  use stateIndex
  use initCond
  use parametersBackground
  use constants
  use filConstants
  use potential
  use cont
  use errorMod
  use dnams
  use variables1, only: rebuild_background,write_full_filament
  use boundmod
  implicit none
  real :: pBack,Bx,Bz
  real :: l,ptotal,ptot
  real :: x0,x1,z0,z1
  real :: dxdz
  real :: Vbg1,Vbg2,Kbg1,Kbg2,filK,fildsob3
  integer :: i,j,k
  real :: Aback,nBack
  real :: c0,c1,c2,u,a,b0,c00,c11,c22,c33,c44
  logical :: notDone
  integer :: bufsize
  real,dimension(x:z,dimj) :: pos
  logical :: isopen
  integer :: ilast,ifirst
  real :: dipx,dipy,dipz
  real :: bilinearInterp
  real,allocatable,dimension(:) :: apexi,apex_z,springx,springz,filBetai,filMassi,midPres,midBz,midV,midN
  real :: dk
  real :: filBeta
  real :: dxBx,dxBz,dzBx,dzBz,dxPtot,dzPtot
  real :: massFil,filVol
  real,allocatable,dimension(:) :: accxx,accxz,acczz,acczx,springxpm,wb_linear,springzpm,springMass
  real,dimension(x:z,dimj) :: acc1,acc2
  real :: springH
  integer :: potGrid,mid1,mid2,mid3,mid4
  real,allocatable,dimension(:) :: potposx,potposz
  real,allocatable,dimension(:,:) :: accxGrid,acczGrid
  real :: apexPotx0,apexPotz0,PotAmp,Nps_tm2003_new,baodensity
  real,allocatable,dimension(:) :: xn,zn,s0,sn,pn,dn,tn,bn,pdgn
  real :: interp1d

  if(firstbgSetup .or. rebuild_background)then
        if(rebuild_background .and. .not. firstbgSetup)then
               if (allocated(Ag)) deallocate(Ag)
               if (allocated(Bxg)) deallocate(Bxg)
               if (allocated(Byg)) deallocate(Byg)
               if (allocated(Bzg)) deallocate(Bzg)
               if (allocated(Pg)) deallocate(Pg)
               if (allocated(Ng)) deallocate(Ng)
               if (allocated(Ptotg)) deallocate(Ptotg)
               if (allocated(testMesh)) deallocate(testMesh)
               if (allocated(oAg)) deallocate(oAg)
               if (allocated(oBxg)) deallocate(oBxg)
               if (allocated(oBzg)) deallocate(oBzg)
               if (allocated(oByg)) deallocate(oByg)
               if (allocated(oPg)) deallocate(oPg)
               if (allocated(xg)) deallocate(xg)
               if (allocated(yg)) deallocate(yg)
               if (allocated(zg)) deallocate(zg)
               if (allocated(oxg)) deallocate(oxg)
               if (allocated(oyg)) deallocate(oyg)
               if (allocated(ozg)) deallocate(ozg)
               if (allocated(dx_Bxg)) deallocate(dx_Bxg)
               if (allocated(dz_Bxg)) deallocate(dz_Bxg)
               if (allocated(dx_Bzg)) deallocate(dx_Bzg)
               if (allocated(dz_Bzg)) deallocate(dz_Bzg)
               if (allocated(dx_Ptotg)) deallocate(dx_Ptotg)
               if (allocated(dz_Ptotg)) deallocate(dz_Ptotg)
               if (allocated(Ag_d)) deallocate(Ag_d)
               if (allocated(Bxg_d)) deallocate(Bxg_d)
               if (allocated(Byg_d)) deallocate(Byg_d)
               if (allocated(Bzg_d)) deallocate(Bzg_d)
               if (allocated(Pg_d)) deallocate(Pg_d)
               if (allocated(Ng_d)) deallocate(Ng_d)
               if (allocated(Ptotg_d)) deallocate(Ptotg_d)
               if (allocated(testMesh_d)) deallocate(testMesh_d)
               if (allocated(PVG_bg)) deallocate(PVG_bg)
               if (allocated(dsob3_bg)) deallocate(dsob3_bg)
               if (allocated(plasGrid)) deallocate(plasGrid)
               if (allocated(plasGridx)) deallocate(plasGridx)
               if (allocated(plasGridz)) deallocate(plasGridz)

        endif

        if(.not.firstbgSetup) plotPVG = .false.
        if(firstbgSetup) pressureScale_in = pressureScale

        if(useGrid)then
               if(rcmdata)then
                      write(*,*) 'setting up rcm-e data'
                      if (tecPlot3D) then
                             call readTecFileRCME(trim(tecPlotDir)//'/'//trim(tecPlotFileName))
                      else
                             call readTecFileRCME2d(trim(tecPlotDir)//'/'//trim(tecPlotFileName))
                      endif
                      allocate(Ag(gridi,gridj))
                      Ag = 0.0
                      allocate(testMesh(gridi,gridj))
                      do i = 1,gridi
                        do j = 1,gridj
                             testMesh(i,j) = sqrt(xg(i)**2+zg(j)**2)
                        enddo
                      enddo
                      call setupDensGrid
                      write(*,*) 'calculating grads'
                      call calculateBG_grads
                      write(*,*) 'done with grads'
               elseif(useAM03ModelGrid) then
                      write(*,*) 'using kubysh data'
                      if(buildKubyshData)then
                             call readAM03
                             write(*,*) 'friction code input ready'
                      else
                             call readAM03Fric
                      endif
                      allocate(Ag(gridi,gridj))
                      Ag = 0.0
                      allocate(testMesh(gridi,gridj))
                      do i = 1,gridi
                        do j = 1,gridj
                             testMesh(i,j) = sqrt(xg(i)**2+zg(j)**2)
                        enddo
                      enddo
                      call setupDensGrid
                      write(*,*) 'calculating grads'
                      call calculateBG_grads
                      write(*,*) 'done with grads'
               elseif(useTsyModelGrid) then
                      write(*,*) 'using Tsyganenko model'
                      if(buildTsyData)then
                             call makeTsyData
                             write(*,*) 'friction code input ready'
                      else
                             call readTsyFric
                      endif
                      allocate(Ag(gridi,gridj))
                      Ag = 0.0
                      allocate(testMesh(gridi,gridj))
                      do i = 1,gridi
                        do j = 1,gridj
                             testMesh(i,j) = sqrt(xg(i)**2+zg(j)**2)
                        enddo
                      enddo
                      call setupDensGrid
                      write(*,*) 'calculating grads'
                      call calculateBG_grads
                      write(*,*) 'done with grads'
               elseif(useAnalyticModelGrid)then
                      write(*,*) 'setting up model grid'
                      call modelGrid
                      allocate(testMesh(gridi,gridj))
                      do i = 1,gridi
                        do j = 1,gridj
                             testMesh(i,j) = sqrt(xg(i)**2+zg(j)**2)
                        enddo
                      enddo
                      call setupDensGrid
                      write(*,*) 'calculating grads'
                      call calculateBG_grads
                      write(*,*) 'done with grads'
               else
                      write(*,*) 'Error: no grid method selected but useGrid is TRUE'
                      stop
               endif
        else
               write(*,*) 'setting up model grid, using analytic'
               call modelGrid
               allocate(testMesh(gridi,gridj))
               do i = 1,gridi
                 do j = 1,gridj
                      testMesh(i,j) = sqrt(xg(i)**2+zg(j)**2)
                 enddo
               enddo
               call setupDensGrid
               write(*,*) 'calculating grads'
               call calculateBG_grads
               write(*,*) 'done with grads'
               useGrid = .false.
        endif
        firstbgSetup = .false.
  endif

  write(*,*) 'grid loaded'
  if(write_full_filament)then
        inquire(Iolength=bufsize) xg
        open(unit=12, file=trim(datapath)//".xg", action="write", status="replace",Recl=bufsize)
          do i=1,gridi
               write(12,*) xg(i)
          end do
        close(unit=12)
        ! write z coordinates to file
        inquire(Iolength=bufsize) zg
        open(unit=12, file=trim(datapath)//".zg", action="write", status="replace",Recl=bufsize)
          do j=1,gridj
               write(12,*) zg(j)
          end do
        close(unit=12)
        write(6,*) 'grid written'
  endif

  write(*,*) 'using grid:',useGrid

  !*****************************Plot PVG***********************************************
  if(plotPVG)then
        write(*,*) 'plotting PVG...'
        PVGgrid = plotPVG_grid
        allocate(PVG_bg(PVGgrid))
        allocate(dsob3_bg(PVGgrid))
        allocate(apexi(PVGgrid))
        allocate(springx(PVGgrid))
        allocate(springxpm(PVGgrid))
        allocate(wb_linear(PVGgrid))
        allocate(springz(PVGgrid))
        allocate(springzpm(PVGgrid))
        allocate(springMass(PVGgrid))
        allocate(apex_z(PVGgrid))
        allocate(filBetai(PVGgrid))
        allocate(filMassi(PVGgrid))
        allocate(accxx(PVGgrid))
        allocate(accxz(PVGgrid))
        allocate(acczx(PVGgrid))
        allocate(acczz(PVGgrid))
        allocate(midPres(PVGgrid))
        allocate(midV(PVGgrid))
        allocate(midBz(PVGgrid))
        allocate(midN(PVGgrid))
        
        do k = 1,PVGgrid
               apexi(k) = plotPVG_xmin + (plotPVG_xmax-plotPVG_xmin)*(PVGgrid-k)/(PVGgrid-1)
               
               if(rcmdata .or. useAnalyticModelGrid)then
                      ifirst = dimj/2+1
                      
                      call tracerRK4(apexi(k),0.0,0.0,1,isopen,pos(x:z,:),ifirst,ilast)
                      if(ilast /= dimj)then
                             write(*,*) 'error in tracerRK4'
                             stop
                      endif
                      
                      state(x:z,ifirst:ilast) = pos(x:z,ifirst:ilast)
                      call tracerRK4(apexi(k),0.0,0.0,-1,isopen,pos(x:z,:),ifirst,ilast)
                      if(ilast /= dimj)then
                             write(*,*) 'error in tracerRK4'
                             stop
                      endif
                      
                      do i = 1,(dimj/2)
                             state(x:z,i) = pos(x:z,dimj-i+1)
                      enddo
               else
                      x0 = apexi(k)
                      do i = 1,plasGridi
                        call line_intersects(x0,zg(1),x0,zg(gridj),plasGridx(i),plasGridz(i),plasGridx(i+1),&
                                             plasGridz(i+1),isopen,x1,z0,1,i,.true.)
                             if (isopen) then
                                    write(*,*) x0,z0
                                    exit
                             endif
                      enddo
                      apex_z(k)=z0
                      ifirst = dimj/2+1
                      call tracerRK4(x0,0.0,z0,1,isopen,pos(x:z,:),ifirst,ilast)
                      if(ilast /= dimj)then
                             write(*,*) 'error in tracerRK4'
                             stop
                      endif
                      
                      state(x:z,ifirst:ilast) = pos(x:z,ifirst:ilast)
                      call tracerRK4(x0,0.0,z0,-1,isopen,pos(x:z,:),ifirst,ilast)
                      if(ilast /= dimj)then
                             write(*,*) 'error in tracerRK4'
                             stop
                      endif
                      
                      do i = 1,(dimj/2)
                             state(x:z,i) = pos(x:z,dimj-i+1)
                      enddo
               endif

               state(vx:vz,:) = 0.0
               do i = 1,dimj
                      state(p,i) = pBack(state(x,i),state(z,i),i)
               end do
               if(useBackgroundDensity)then
                      if(use_gallagherDensity)then
                        call gallagher(state(x,dimj/2+1),0.0,state(z,dimj/2+1),dBack)
                      endif
                      if(use_tm2003_density)then
                        dBack = Nps_tm2003_new(-sqrt(state(x,dimj/2+1)**2+state(z,dimj/2+1)**2),&
                                                0.0,sw_density,sw_velocity,bzimf)
                      endif
                      if(use_bao_density)then
                        dBack = baodensity(state(x,dimj/2+1),0.0,state(z,dimj/2+1))
                      endif
                      state(d,:) = dBack
                      state(t,:) = state(p,:)/(kb*state(d,:))
               else
                      state(t,:) = tBack*tkev
                      state(d,:) = state(p,:)/(kb*state(t,:))
               endif
               do i = 1,dimj
                      state(b,i) = sqrt(2.0*mu0*ptotal(state(x,i),state(z,i),i)-2.0*mu0*state(p,i))
               end do
               state(pdg,:) = state(p,:)*(state(d,:))**(-gamma)
               call findMf2(state,state(mf,:))
               

               PVG_bg(k) = filK(state)
               dsob3_bg(k) = fildsob3(state) ! integral ds/b^3
               filBetai(k) = filBeta(state)
               filMassi(k) = massFil(state)
               springz(k) = filVol(state)*dxBz(apexi(k),apex_z(k),-1)*dzBx(apexi(k),apex_z(k),-1)/mu0
               
               midV(k) = filVol(state)
               midBz(k) = Bz(state(x,dimj/2+1),state(z,dimj/2+1),-1)
               midPres(k) = pBack(state(x,dimj/2+1),state(z,dimj/2+1),-1)
               midN(k) = nBack(state(x,dimj/2+1),state(z,dimj/2+1),-1)
               springMass(k) = sum(state(mf,:))
               
               
               !*****************************************************************************
               !pretty damn close for dimj=299
               !springH = 0.5
               !j= dimj/2-85
               !*****************************************************************************
               !*****************************************************************************
               springH = 0.5
               j= dimj/2-floor(fmsection*dimj)
               !*****************************************************************************
               oldstate = state
               
               state = 0.0
               do i = 1,dimj
                      state(x,i) = sin(Pi*(i-1)/(dimj-1))*springH/2.0
               enddo
               newstate = oldstate+state
               !call algebra1(newstate,0.0,0.0,newstate,-1)
               !call accel4(newstate,acc1,1,dimj)
               acc1 = 0.0
               newstate = oldstate-state
               !call algebra1(newstate,0.0,0.0,newstate,-1)
               !call accel4(newstate,acc2,1,dimj)
               acc2 = 0.0

               c0 = 0.0
               c1 = 0.0
               do i = 1+j,dimj-j
                      c0 = c0 + newstate(mf,i)*(60**2)*(acc1(x,i)-acc2(x,i))/springH
               enddo
               do i = 1+j,dimj-j
                      c1 = c1 + newstate(mf,i)
               enddo
               accxx(k) = c0/c1
               
               state = 0.0
               do i = 1,dimj
                      state(z,i) = sin(Pi*(i-1)/(dimj-1))*springH/2.0
               enddo
               newstate = oldstate+state
               ! call algebra1(newstate,0.0,0.0,newstate,-1)
               ! call accel4(newstate,acc1,1,dimj)
               newstate = oldstate-state
               ! call algebra1(newstate,0.0,0.0,newstate,-1)
               ! call accel4(newstate,acc2,1,dimj)
               c0 = 0.0
               c1 = 0.0
               do i = 1+j,dimj-j
                      c0 = c0 + newstate(mf,i)*(60**2)*(acc1(z,i)-acc2(z,i))/springH
               enddo
               do i = 1+j,dimj-j
                      c1 = c1 + newstate(mf,i)
               enddo
               acczz(k) = c0/c1
               
               !*****************************************************************************
               springH = 0.5
               j= dimj/2
               !*****************************************************************************
               
               state = 0.0
               do i = 1,dimj
                      state(z,i) = sin(Pi*(i-1)/(dimj-1))*springH/2.0
               enddo
               newstate = oldstate+state
               ! call algebra1(newstate,0.0,0.0,newstate,-1)
               ! call accel4(newstate,acc1,1,dimj)
               newstate = oldstate-state
               ! call algebra1(newstate,0.0,0.0,newstate,-1)
               ! call accel4(newstate,acc2,1,dimj)
               c0 = 0.0
               c1 = 0.0
               do i = 1+j,dimj-j
                      c0 = c0 + newstate(mf,i)*(60**2)*(acc1(x,i)-acc2(x,i))/springH
               enddo
               do i = 1+j,dimj-j
                      c1 = c1 + newstate(mf,i)
               enddo
               accxz(k) = c0/c1
               
               state = 0.0
               do i = 1,dimj
                      state(x,i) = sin(Pi*(i-1)/(dimj-1))*springH/2.0
               enddo
               newstate = oldstate+state
               ! call algebra1(newstate,0.0,0.0,newstate,-1)
               ! call accel4(newstate,acc1,1,dimj)
               newstate = oldstate-state
               ! call algebra1(newstate,0.0,0.0,newstate,-1)
               ! call accel4(newstate,acc2,1,dimj)
               c0 = 0.0
               c1 = 0.0
               do i = 1+j,dimj-j
                      c0 = c0 + newstate(mf,i)*(60**2)*(acc1(z,i)-acc2(z,i))/springH
               enddo
               do i = 1+j,dimj-j
                      c1 = c1 + newstate(mf,i)
               enddo
               acczx(k) = c0/c1
                                    
               state = oldstate
               !*****************************************************************************
        enddo
        do k = 2,PVGgrid-1
               c00 = pBack(apexi(k),apex_z(k),-1)
               c11 = abs((PVG_bg(k+1)-PVG_bg(k-1))/(apexi(k+1)-apexi(k-1)))
               c22 = PVG_bg(k)*Bz(apexi(k),apex_z(k),-1)
               c33 = 1.0+gamma*filBetai(k)/2.0
               springx(k) = Pi*c00*c11/(c22*c33)
               springx(k) = abs(springx(k))
               springxpm(k) = springx(k)/(springMass(k)*mi*Re*Re*cm)
               wb_linear(k) = sqrt(springxpm(k))
        enddo
        
  endif

  if(forceModel)then
        !real,dimension(100) :: potentialx,potentialz
        !real :: apexPotx0,apexPotz0,PotAmp
        potGrid = 21
        allocate(accxGrid(potGrid,potGrid))
        allocate(acczGrid(potGrid,potGrid))
        allocate(potposx(potGrid))
        allocate(potposz(potGrid))
        
        apexPotx0 = apex
        
        if(rcmdata .or. useAnalyticModelGrid)then
               ifirst = dimj/2+1
               apexPotz0 = 0.0
               call tracerRK4(apexPotx0,0.0,0.0,1,isopen,pos(x:z,:),ifirst,ilast)
               if(ilast /= dimj)then
                      write(*,*) 'error in tracerRK4'
                      stop
               endif
               
               state(x:z,ifirst:ilast) = pos(x:z,ifirst:ilast)
               call tracerRK4(apexPotx0,0.0,0.0,-1,isopen,pos(x:z,:),ifirst,ilast)
               if(ilast /= dimj)then
                      write(*,*) 'error in tracerRK4'
                      stop
               endif
               
               do i = 1,(dimj/2)
                      state(x:z,i) = pos(x:z,dimj-i+1)
               enddo
        else
               x0 = -abs(apexPotx0)
               do i = 1,plasGridi
                 call line_intersects(x0,zg(1),x0,zg(gridj),plasGridx(i),plasGridz(i),&
                                      plasGridx(i+1),plasGridz(i+1),isopen,x1,z0,1,i,.true.)
                      if (isopen) then
                             write(*,*) x0,z0
                             exit
                      endif
               enddo
               apexPotz0=z0
               ifirst = dimj/2+1
               call tracerRK4(x0,0.0,z0,1,isopen,pos(x:z,:),ifirst,ilast)
               if(ilast /= dimj)then
                      write(*,*) 'error in tracerRK4'
                      stop
               endif
               
               state(x:z,ifirst:ilast) = pos(x:z,ifirst:ilast)
               call tracerRK4(x0,0.0,z0,-1,isopen,pos(x:z,:),ifirst,ilast)
               if(ilast /= dimj)then
                      write(*,*) 'error in tracerRK4'
                      stop
               endif
               
               do i = 1,(dimj/2)
                      state(x:z,i) = pos(x:z,dimj-i+1)
               enddo
        endif
        state(vx:vz,:) = 0.0
        do i = 1,dimj
               state(p,i) = pBack(state(x,i),state(z,i),i)
        end do
        if(useBackgroundDensity)then
               if(use_gallagherDensity)then
                      call gallagher(state(x,dimj/2+1),0.0,state(z,dimj/2+1),dBack)
               endif
               if(use_tm2003_density)then
                      dBack = Nps_tm2003_new(-sqrt(state(x,dimj/2+1)**2+state(z,dimj/2+1)**2),&
                              0.0,sw_density,sw_velocity,bzimf)
               endif
               if(use_bao_density)then
                      dBack = baodensity(state(x,dimj/2+1),0.0,state(z,dimj/2+1))
               endif
               state(d,:) = dBack
               state(t,:) = state(p,:)/(kb*state(d,:))
        else
               state(t,:) = tBack*tkev
               state(d,:) = state(p,:)/(kb*state(t,:))
        endif
        do i = 1,dimj
               state(b,i) = sqrt(2.0*mu0*ptotal(state(x,i),state(z,i),i)-2.0*mu0*state(p,i))
        end do
        state(pdg,:) = state(p,:)*(state(d,:))**(-gamma)
        call findMf2(state,state(mf,:))
        
        oldstate = state
        mid1 = dimj/2-floor(fmsection*dimj)
        PotAmp = 2.0
        do j = 1,potGrid
          do k = 1,potGrid
               c00 = 2.0*PotAmp*(j-1)/(potGrid-1) - PotAmp
               c11 = 2.0*PotAmp*(k-1)/(potGrid-1) - PotAmp
               potposx(j) = c00
               potposz(k) = c11
               state = 0.0
               do i = 1,dimj
                      state(x,i) = sin(Pi*(i-1)/(dimj-1))*c00
               enddo
               do i = 1,dimj
                      state(z,i) = sin(Pi*(i-1)/(dimj-1))*c11
               enddo
               newstate = oldstate+state
               call algebra1(newstate,0.0,0.0,newstate,-1)
               call accel4(newstate,acc1,1,dimj)
               c0 = 0.0
               c1 = 0.0
               do i = 1+mid1,dimj-mid1
                      c0 = c0 + newstate(mf,i)*(60**2)*acc1(x,i)
               enddo
               do i = 1+mid1,dimj-mid1
                      c1 = c1 + newstate(mf,i)
               enddo
               accxGrid(j,k) = c0/c1
               c0 = 0.0
               c1 = 0.0
               do i = 1+mid1,dimj-mid1
                      c0 = c0 + newstate(mf,i)*(60**2)*acc1(z,i)
               enddo
               do i = 1+mid1,dimj-mid1
                      c1 = c1 + newstate(mf,i)
               enddo
               acczGrid(j,k) = c0/c1
               !*****************************************************************************
          enddo
        enddo
        
        state = oldstate
        open(unit=11, file=trim(datapath)//".accGrid.dat", status="replace")
        write(11,*) potGrid
        do j = 1,potGrid
          do k = 1,potGrid
               write(11,*) j,k,potposx(j),potposz(k),accxGrid(j,k),acczGrid(j,k)
          enddo
        enddo
        close(11)
        
        call littleSolver2d(potposx,potposz,accxGrid,acczGrid,potGrid,x0,z0)
        
  endif

  if(plotPVG)then
          write(6,*)' writing out pv^gamma information to *.PVG.dat'

          open(unit=11, file=trim(datapath)//".PVG.dat", status="replace")
          do i = 1,PVGgrid
                 write(11,*) apexi(i), apex_z(i),PVG_bg(i),pBack(apexi(i),apex_z(i),-1),Bz(apexi(i),apex_z(i),-1),dsob3_bg(i)
          enddo
          close(11)

        ! open(unit=11, file=trim(datapath)//".PVGp.dat", status="replace")
        ! do i = 2,PVGgrid-1
        !        write(11,*) apexi(i),(PVG_bg(i+1)-PVG_bg(i-1))/(apexi(i+1)-apexi(i-1)) 
        ! enddo
        ! close(11)

        ! open(unit=11, file=trim(datapath)//".midN.dat", status="replace")
        ! do i = 1,PVGgrid
        !        write(11,*) apexi(i),midN(i)
        ! enddo
        ! close(11)

         open(unit=11, file=trim(datapath)//".mass.dat", status="replace")
         do i = 2,PVGgrid-1
                write(11,*) apexi(i),filMassi(i)
         enddo
         close(11)

        ! open(unit=11, file=trim(datapath)//".springz.dat", status="replace")
        ! do i = 2,PVGgrid-1
        !        write(11,*) apexi(i),springz(i) 
        ! enddo
        ! close(11)

        ! open(unit=11, file=trim(datapath)//".springx.dat", status="replace")
        ! do i = 2,PVGgrid-1
        !        write(11,*) apexi(i),springx(i)
        ! enddo
        ! close(11)

        ! open(unit=11, file=trim(datapath)//".periodvsx.dat", status="replace")
        ! do i = 2,PVGgrid-1
        !        write(11,*) apexi(i),2.0*Pi*nano*sqrt(mi*nano*cm/Re)*sqrt(filMassi(i)/springx(i))/60.0, &
        !        & 2.0*Pi*nano*sqrt(mi*nano*cm/Re)*sqrt(filMassi(i)/springz(i))/60.0
        ! enddo
        ! close(11)

        ! open(unit=11, file=trim(datapath)//".cdiffSpring.dat", status="replace")
        ! do i = 1,PVGgrid
        !        write(11,*) apexi(i),accxx(i),acczz(i),accxz(i),acczx(i)
        ! enddo
        ! close(11)

        ! open(unit=11, file=trim(datapath)//".cdiffPeriods.dat", status="replace")
        ! do i = 1,PVGgrid
        !        c0 = abs(accxx(i))/2.0+abs(acczz(i))/2.0 &
        !        & +sqrt((abs(accxx(i))-abs(acczz(i)))**2+4.0*(abs(accxz(i))/2.0+abs(acczx(i))/2.0)**2)
        !        c1 = abs(accxx(i))/2.0+abs(acczz(i))/2.0 &
        !        & -sqrt((abs(accxx(i))-abs(acczz(i)))**2+4.0*(abs(accxz(i))/2.0+abs(acczx(i))/2.0)**2)
        !        !write(11,*) apexi(i),accxx(i),acczz(i),accxz(i),acczx(i)
        !        write(11,*) apexi(i),2.0*Pi/sqrt(c0),2.0*Pi/sqrt(c1),2.0*Pi/sqrt(abs(accxx(i))),2.0*Pi/sqrt(abs(acczz(i)))
        ! enddo
        ! close(11)

        ! open(unit=11, file=trim(datapath)//".cdiffSpring.dat", status="replace")
        ! do i = 1,PVGgrid
        !        write(11,*) apexi(i),-accxx(i),-acczz(i),-accxz(i),-acczx(i)
        ! enddo
        ! close(11)

        ! open(unit=11, file=trim(datapath)//".cdiffPsi.dat", status="replace")
        ! do i = 1,PVGgrid
        !        write(11,*) apexi(i),1/sqrt(1+(accxx(i)-acczz(i))**2/(4.0*(accxz(i)/2.0+acczx(i)/2.0)**2))
        ! enddo
        ! close(11)

        open(unit=11, file=trim(datapath)//".midPointData.dat", status="replace")
        do i = 2,PVGgrid-1
               write(11,*) apexi(i),apex_z(i),midPres(i),midBz(i),midV(i),midN(i),PVG_bg(i),&
                 (PVG_bg(i+1)-PVG_bg(i-1))/(apexi(i+1)-apexi(i-1)), &
               & filBetai(i),springxpm(i),wb_linear(i)
        enddo
        close(11)

        ! open(unit=11, file= trim(tsyModelDir)//"/"//"withN_"//trim(tsyModel), status="replace")
        ! write(11,*) 'TITLE="tsy data"'
        ! write(11,*) 'VARIABLES="I" "J" "X" "Z" "Bx" "Bz" "P" "N"'
        ! do i = 1,gridi
        !  do j = 1,gridj
        !        write(11,*) i,j,xg(i),zg(j),Bxg(i,j),Bzg(i,j),Pg(i,j),Ng(i,j)
        !  enddo
        ! enddo
        ! close(11)
  endif

  !****************************************************************************

  write(*,*) 'tracing field line'

  if (use_dPdK)then
        if(dPdK_fixed)then
               if(rcmdata .or. useAnalyticModelGrid)then
                      ifirst = dimj/2+1
                      if(useKpk)then
                             call tracerRK4(Kpk,0.0,0.0,1,isopen,pos(x:z,:),ifirst,ilast)
                      else
                             call tracerRK4(apex,0.0,0.0,1,isopen,pos(x:z,:),ifirst,ilast)
                      endif
                      
                      if(ilast /= dimj)then
                             write(*,*) 'error in tracerRK4'
                             stop
                      endif
                      state(x:z,ifirst:ilast) = pos(x:z,ifirst:ilast)
                      if(useKpk)then
                             call tracerRK4(Kpk,0.0,0.0,-1,isopen,pos(x:z,:),ifirst,ilast)
                      else
                             call tracerRK4(apex,0.0,0.0,-1,isopen,pos(x:z,:),ifirst,ilast)
                      endif
                      if(ilast /= dimj)then
                             write(*,*) 'error in tracerRK4'
                             stop
                      endif
                      do i = 1,(dimj/2)
                             state(x:z,i) = pos(x:z,dimj-i+1)
                      enddo
                      if(fac > 0.0)call remap_state(state(x,:),state(z,:),dimj,fac)

               else
                      x0 = -abs(apex)
                      do i = 1,plasGridi
                        call line_intersects(x0,zg(1),x0,zg(gridj),plasGridx(i),&
                                             plasGridz(i),plasGridx(i+1),plasGridz(i+1),isopen,x1,z0,1,i,.true.)
                             if (isopen) then
                                    !write(*,*) x0,z0
                                    exit
                             endif
                      enddo
                      ifirst = dimj/2+1
                      call tracerRK4(x0,0.0,z0,1,isopen,pos(x:z,:),ifirst,ilast)
                      if(ilast /= dimj)then
                             write(*,*) 'error in tracerRK4'
                             stop
                      endif
               
                      state(x:z,ifirst:ilast) = pos(x:z,ifirst:ilast)
                      call tracerRK4(x0,0.0,z0,-1,isopen,pos(x:z,:),ifirst,ilast)
                      if(ilast /= dimj)then
                             write(*,*) 'error in tracerRK4'
                             stop
                      endif
                      do i = 1,(dimj/2)
                             state(x:z,i) = pos(x:z,dimj-i+1)
                      enddo

                     if(fac > 0.0)call remap_state(state(x,:),state(z,:),dimj,fac) 

               endif
               state(vx:vz,:) = 0.0
               do i = 1,dimj
                      state(p,i) = pBack(state(x,i),state(z,i),i)
               end do
               
               if(useBackgroundDensity)then
                      if(use_gallagherDensity)then
                             call gallagher(state(x,dimj/2+1),0.0,state(z,dimj/2+1),dBack)
                      endif
                      if(use_tm2003_density)then
                       dBack = Nps_tm2003_new(-sqrt(state(x,dimj/2+1)**2+state(z,dimj/2+1)**2),&
                               0.0,sw_density,sw_velocity,bzimf)
                      endif
                      if(use_bao_density)then
                       dBack = baodensity(state(x,dimj/2+1),0.0,state(z,dimj/2+1))
                      endif
                      state(d,:) = dBack
                      state(t,:) = state(p,:)/(kb*state(d,:))
               else
                      state(t,:) = tBack*tkev
                      state(d,:) = state(p,:)/(kb*state(t,:))
               endif
               
               do i = 1,dimj
                      state(b,i) = sqrt(2.0*mu0*ptotal(state(x,i),state(z,i),i)-2.0*mu0*state(p,i))
               end do
               state(pdg,:) = state(p,:)*(state(d,:))**(-gamma)
               call findMf2(state,state(mf,:))
               
               K_peak = filK(state)
               
               
               dK = filK(state)*(1.0-pressureScale_in)
        else
               PVGgrid = 3
               if(.not. plotPVG)then
                      allocate(PVG_bg(PVGgrid))
                      allocate(apexi(PVGgrid))
               endif
               apexi(1) = apex-ampli
               apexi(2) = apex
               apexi(3) = apex+ampli
               do k = 1,3
                      if(rcmdata .or. useAnalyticModelGrid)then
                             ifirst = dimj/2+1
                             call tracerRK4(apexi(k),0.0,0.0,1,isopen,pos(x:z,:),ifirst,ilast)
                             if(ilast /= dimj)then
                                    write(*,*) 'error in tracerRK4'
                                    stop
                             endif
                             state(x:z,ifirst:ilast) = pos(x:z,ifirst:ilast)
                             call tracerRK4(apexi(k),0.0,0.0,-1,isopen,pos(x:z,:),ifirst,ilast)
                             if(ilast /= dimj)then
                                    write(*,*) 'error in tracerRK4'
                                    stop
                             endif
                             do i = 1,(dimj/2)
                                    state(x:z,i) = pos(x:z,dimj-i+1)
                             enddo
                      else
                             x0 = -abs(apexi(k))
                             do i = 1,plasGridi
                               call line_intersects(x0,zg(1),x0,zg(gridj),plasGridx(i),&
                                                    plasGridz(i),plasGridx(i+1),plasGridz(i+1),isopen,x1,z0,1,i,.true.)
                                    if (isopen) then
                                           !write(*,*) x0,z0
                                           exit
                                    endif
                             enddo
                             ifirst = dimj/2+1
                             call tracerRK4(x0,0.0,z0,1,isopen,pos(x:z,:),ifirst,ilast)
                             if(ilast /= dimj)then
                                    write(*,*) 'error in tracerRK4'
                                    stop
                             endif
                             state(x:z,ifirst:ilast) = pos(x:z,ifirst:ilast)
                             call tracerRK4(x0,0.0,z0,-1,isopen,pos(x:z,:),ifirst,ilast)
                             if(ilast /= dimj)then
                                    write(*,*) 'error in tracerRK4'
                                    stop
                             endif
                             do i = 1,(dimj/2)
                                    state(x:z,i) = pos(x:z,dimj-i+1)
                             enddo
                            if(fac > 0.0)call remap_state(state(x,:),state(z,:),dimj,fac) 
                      endif
                      state(vx:vz,:) = 0.0
                      do i = 1,dimj
                             state(p,i) = pBack(state(x,i),state(z,i),i)
                      end do
                      if(useBackgroundDensity)then
                             if(use_gallagherDensity)then
                               call gallagher(state(x,dimj/2+1),0.0,state(z,dimj/2+1),dBack)
                             endif
                             if(use_tm2003_density)then
                               dBack = Nps_tm2003_new(-sqrt(state(x,dimj/2+1)**2+state(z,dimj/2+1)**2),&
                                       0.0,sw_density,sw_velocity,bzimf)
                             endif
                             if(use_bao_density)then
                               dBack = baodensity(state(x,dimj/2+1),0.0,state(z,dimj/2+1))
                             endif
                             state(d,:) = dBack
                             state(t,:) = state(p,:)/(kb*state(d,:))
                      else
                             state(t,:) = tBack*tkev
                             state(d,:) = state(p,:)/(kb*state(t,:))
                      endif
                      do i = 1,dimj
                             state(b,i) = sqrt(2.0*mu0*ptotal(state(x,i),state(z,i),i)-2.0*mu0*state(p,i))
                      end do
                      state(pdg,:) = state(p,:)*(state(d,:))**(-gamma)
                      call findMf2(state,state(mf,:))
                      PVG_bg(k) = filK(state)
               enddo
               
               if(PVG_bg(1)>PVG_bg(2).and.PVG_bg(2)>PVG_bg(3))then
                      dK = PVG_bg(1)-PVG_bg(2)
                      apex = apexi(1)
                      write(*,*) 'slope is neg'
               elseif(PVG_bg(1)<PVG_bg(2).and.PVG_bg(2)<PVG_bg(3))then
                      dK = PVG_bg(3)-PVG_bg(2)
                      apex = apexi(3)
                      write(*,*) 'slope is pos'
               elseif(PVG_bg(1)>PVG_bg(2).and.PVG_bg(2)<PVG_bg(3))then
                      dK = PVG_bg(3)-PVG_bg(2)
                      apex = apexi(3)
                      write(*,*) 'point is min'
               elseif(PVG_bg(1)<PVG_bg(2).and.PVG_bg(2)>PVG_bg(3))then
                      dK = PVG_bg(2)-PVG_bg(1)
                      apex = apexi(1)
                      write(*,*) 'point is max'
               else
                      dK = PVG_bg(2)*(1-pressureScale_in)
               endif
        endif
        
        !write(*,*) filK(state)
        
        if(rcmdata .or. useAnalyticModelGrid)then
               ifirst = dimj/2+1
               call tracerRK4(apex,0.0,0.0,1,isopen,pos(x:z,:),ifirst,ilast)
               if(ilast /= dimj)then
                      write(*,*) 'error in tracerRK4'
                      stop
               endif
               
               state(x:z,ifirst:ilast) = pos(x:z,ifirst:ilast)
               call tracerRK4(apex,0.0,0.0,-1,isopen,pos(x:z,:),ifirst,ilast)
               if(ilast /= dimj)then
                      write(*,*) 'error in tracerRK4'
                      stop
               endif
               
               do i = 1,(dimj/2)
                      state(x:z,i) = pos(x:z,dimj-i+1)
               enddo
        else
               x0 = -abs(apex)
               do i = 1,plasGridi
                 call line_intersects(x0,zg(1),x0,zg(gridj),plasGridx(i),plasGridz(i),&
                                      plasGridx(i+1),plasGridz(i+1),isopen,x1,z0,1,i,.true.)
                      if (isopen) then
                             !write(*,*) x0,z0
                             exit
                      endif
               enddo
               ifirst = dimj/2+1
               call tracerRK4(x0,0.0,z0,1,isopen,pos(x:z,:),ifirst,ilast)
               if(ilast /= dimj)then
                      write(*,*) 'error in tracerRK4'
                      stop
               endif
               
               state(x:z,ifirst:ilast) = pos(x:z,ifirst:ilast)
               call tracerRK4(x0,0.0,z0,-1,isopen,pos(x:z,:),ifirst,ilast)
               if(ilast /= dimj)then
                      write(*,*) 'error in tracerRK4'
                      stop
               endif
               
               do i = 1,(dimj/2)
                      state(x:z,i) = pos(x:z,dimj-i+1)
               enddo
              if(fac > 0.0)call remap_state(state(x,:),state(z,:),dimj,fac) 
        endif
        state(vx:vz,:) = 0.0
        do i = 1,dimj
               state(p,i) = pBack(state(x,i),state(z,i),i)
        end do
        if(useBackgroundDensity)then
               if(use_gallagherDensity)then
                 call gallagher(state(x,dimj/2+1),0.0,state(z,dimj/2+1),dBack)
               endif
               if(use_tm2003_density)then
                 dBack = Nps_tm2003_new(-sqrt(state(x,dimj/2+1)**2+state(z,dimj/2+1)**2),0.0,sw_density,sw_velocity,bzimf)
               endif
               if(use_bao_density)then
                 dBack = baodensity(state(x,dimj/2+1),0.0,state(z,dimj/2+1))
               endif
               state(d,:) = dBack
               state(t,:) = state(p,:)/(kb*state(d,:))
        else
               state(t,:) = tBack*tkev
               state(d,:) = state(p,:)/(kb*state(t,:))
        endif
        do i = 1,dimj
               state(b,i) = sqrt(2.0*mu0*ptotal(state(x,i),state(z,i),i)-2.0*mu0*state(p,i))
        end do
        state(pdg,:) = state(p,:)*(state(d,:))**(-gamma)
        call findMf2(state,state(mf,:))
        
        
        c0 = filK(state)
        write(*,*) 'Background K=PV^Gamma:',c0
        pressureScale = 1.0
        oldState = state
        !iterate to solve for proper pressure depletion
        do while(.true.)
               c1 = filK(state)
               write(*,*) 'filK',c1,'filK/bgK',c1/c0
               if(abs(dK-abs(c0-c1))<1.0e-10)then
                      exit
               elseif(c0-c1<dK)then !elseif(c0-c1>=0)then
                      pressureScale = pressureScale - abs(dK-abs(c0-c1))/(c1*(1+gamma*filBeta(state)/2.0))
               elseif(c0-c1>dK)then !elseif(c0-c1<0)then
                      pressureScale = pressureScale + abs(dK-abs(c0-c1))/(c1*(1+gamma*filBeta(state)/2.0))
               endif
               state(p,:) = pressureScale*oldState(p,:)
               state(vx:vz,:) = 0.0
               if(useBackgroundDensity)then
                      if(use_gallagherDensity)then
                             call gallagher(state(x,dimj/2+1),0.0,state(z,dimj/2+1),dBack)
                      endif
                      if(use_tm2003_density)then
                             dBack = Nps_tm2003_new(-sqrt(state(x,dimj/2+1)**2+state(z,dimj/2+1)**2),&
                                     0.0,sw_density,sw_velocity,bzimf)
                      endif
                      if(use_bao_density)then
                             dBack = baodensity(state(x,dimj/2+1),0.0,state(z,dimj/2+1))
                      endif
                      state(d,:) = dBack
                      state(t,:) = state(p,:)/(kb*state(d,:))
               else
                      state(t,:) = tBack*tkev
                      state(d,:) = state(p,:)/(kb*state(t,:))
               endif
               do i = 1,dimj
                      state(b,i) = sqrt(2.0*mu0*ptotal(state(x,i),state(z,i),i)-2.0*mu0*state(p,i))
               end do
               state(pdg,:) = state(p,:)*(state(d,:))**(-gamma)
               call findMf2(state,state(mf,:))
        enddo
  endif



  if(rcmdata .or. useAnalyticModelGrid)then
        ifirst = dimj/2+1
        if (trace_scale_step)then
               call tracerRK4a(apex,0.0,0.0,1,isopen,pos(x:z,:),ifirst,ilast)
        else
               call tracerRK4(apex,0.0,0.0,1,isopen,pos(x:z,:),ifirst,ilast)
        endif
        if(ilast /= dimj)then
               write(*,*) 'error in tracerRK4'
               stop
        endif
        write(*,*) 'isopen',isopen
        write(*,*) 'tracing field line'
        state(x:z,ifirst:ilast) = pos(x:z,ifirst:ilast)
        if (trace_scale_step)then
               call tracerRK4a(apex,0.0,0.0,-1,isopen,pos(x:z,:),ifirst,ilast)
        else
               call tracerRK4(apex,0.0,0.0,-1,isopen,pos(x:z,:),ifirst,ilast)
        endif
        if(ilast /= dimj)then
               write(*,*) 'error in tracerRK4'
               stop
        endif
        write(*,*) 'isopen',isopen
        do i = 1,(dimj/2)
               state(x:z,i) = pos(x:z,dimj-i+1)
        enddo

        if(fac > 0.0)call remap_state(state(x,:),state(z,:),dimj,fac) 
  else
        x0 = -abs(apex)
        do i = 1,plasGridi
         call line_intersects(x0,zg(1),x0,zg(gridj),plasGridx(i),plasGridz(i),&
                              plasGridx(i+1),plasGridz(i+1),isopen,x1,z0,1,i,.true.)
               if (isopen) then
                      write(*,*) x0,z0
                      exit
               endif
        enddo
        ifirst = dimj/2+1
        if (trace_scale_step) then
               call tracerRK4a(x0,0.0,z0,1,isopen,pos(x:z,:),ifirst,ilast)
        else
               call tracerRK4(x0,0.0,z0,1,isopen,pos(x:z,:),ifirst,ilast)
        endif
        if(ilast /= dimj)then
               write(*,*) 'error in tracerRK4'
               stop
        endif
        write(*,*) 'isopen',isopen
        write(*,*) 'tracing field line'
        state(x:z,ifirst:ilast) = pos(x:z,ifirst:ilast)
        if (trace_scale_step) then
               call tracerRK4a(x0,0.0,z0,-1,isopen,pos(x:z,:),ifirst,ilast)
        else
               call tracerRK4(x0,0.0,z0,-1,isopen,pos(x:z,:),ifirst,ilast)
        endif
        if(ilast /= dimj)then
               write(*,*) 'error in tracerRK4'
               stop
        endif
        write(*,*) 'isopen',isopen
        do i = 1,(dimj/2)
               state(x:z,i) = pos(x:z,dimj-i+1)
        enddo
       if(fac > 0.0)call remap_state(state(x,:),state(z,:),dimj,fac) 
  endif

  write(6,*) 'setting up...'

  ! Pressure underpopulated by pressureScale < 1, (nPa)


  do i = 1,dimj
        state(p,i) = pBack(state(x,i),state(z,i),-1)
  end do
  state(vx:vz,:) = 0.0
  if(useBackgroundDensity)then
        if(use_gallagherDensity)then
               call gallagher(state(x,dimj/2+1),0.0,state(z,dimj/2+1),dBack)
        endif
        if(use_tm2003_density)then
               dBack = Nps_tm2003_new(-sqrt(state(x,dimj/2+1)**2+state(z,dimj/2+1)**2),0.0,sw_density,sw_velocity,bzimf)
        endif
        if(use_bao_density)then
               dBack = baodensity(state(x,dimj/2+1),0.0,state(z,dimj/2+1))
        endif
        state(d,:) = dBack
        state(t,:) = state(p,:)/(kb*state(d,:))
  else
        state(t,:) = tBack*tkev
        state(d,:) = state(p,:)/(kb*state(t,:))
  endif
  !find the B filament using P_tot balance equation (nT)
  do i = 1,dimj
        state(b,i) = sqrt(2.0*mu0*ptotal(state(x,i),state(z,i),-1)-2.0*mu0*state(p,i))
  end do
  ! p/d^gamma for updating density (nPa*cm^3gamma)
  state(pdg,:) = state(p,:)*(state(d,:))**(-gamma)
  ! find the mass of each element by integrating dm = Rho/B ds to adjacent elements 
  ! [imposed to be constant] ( 1/(nT*cm^2) ) inverse flux
  call findMf2(state,state(mf,:))

  time = 0.0

  do i = 1,dimj
        state(p,i) = pressureScale*pBack(state(x,i),state(z,i),-1)
  end do
  state(vx:vz,:) = 0.0
  if(useBackgroundDensity)then
        if(use_gallagherDensity)then
               call gallagher(state(x,dimj/2+1),0.0,state(z,dimj/2+1),dBack)
        endif
        if(use_tm2003_density)then
               dBack = Nps_tm2003_new(-sqrt(state(x,dimj/2+1)**2+state(z,dimj/2+1)**2),0.0,sw_density,sw_velocity,bzimf)
        endif
        if(use_bao_density)then
               dBack = baodensity(state(x,dimj/2+1),0.0,state(z,dimj/2+1))
        endif
        if(use_gallagherDensity)then
               do i = 1,dimj
                      call gallagher(state(x,i),0.0,state(z,i),state(d,i))
               enddo
        else
               state(d,:) = dBack
        endif
        state(t,:) = state(p,:)/(kb*state(d,:))
  else
        state(t,:) = tBack*tkev
        state(d,:) = state(p,:)/(kb*state(t,:))
  endif
  !find the B filament using P_tot balance equation (nT)
  do i = 1,dimj
        state(b,i) = sqrt(2.0*mu0*ptotal(state(x,i),state(z,i),-1)-2.0*mu0*state(p,i))
  end do
  ! p/d^gamma for updating density (nPa*cm^3gamma)
  state(pdg,:) = state(p,:)*(state(d,:))**(-gamma)
  ! find the mass of each element by integrating dm = Rho/B ds to adjacent elements 
  ! [imposed to be constant] ( 1/(nT*cm^2) ) inverse flux
  call findMf2(state,state(mf,:))

  write(*,*) 'done'

  if(plotPVG .or. use_dPdK .and. .not. dPdK_fixed)then
        deallocate(PVG_bg)
        deallocate(dsob3_bg)
        deallocate(apexi)
  endif

  if(.false.)then
        do i = 1,dimj
               state(x,i) = state(x,i)-sin(Pi*(i-1)/(dimj-1))*0.0
        enddo
        do i = 1,dimj
               state(z,i) = state(z,i)+sin(Pi*(i-1)/(dimj-1))*0.5
        enddo

       if(fac > 0.0)call remap_state(state(x,:),state(z,:),dimj,fac) 
  endif

  oldState = state

  s_0 = state

  !write(*,*) pressureScale,filK(state),filK(oldState),filK(state)/filK(oldState),filK(state)/K_peak,apex
  !write(*,*) state(p,dimj/2+1),state(t,dimj/2+1)/tkev,state(d,dimj/2+1)
  !write(*,*) filBeta(state),filBeta(oldState)
  !write(*,*) 2*mu0*pBack(-12.51,0.0,-1)/(Bx(-12.51,0.0,-1)**2+Bz(-12.51,0.0,-1)**2)

  !write(*,*) alpha,A0,ht

  call filEigen

end subroutine initialConditions1

subroutine remap_state(x,z,dimj,fac)
        use dnams,only : datapath
        implicit none
        integer, intent(in) ::dimj
        real :: x(dimj),z(dimj)
        real :: s0(dimj),sn(dimj),xn(dimj),zn(dimj)
        real, intent(in) :: fac
        integer :: j
        real :: pi,a0,interp1d

        if(fac == 0.0)return

        pi = acos(-1.0)

         s0(1) = 0.0
         do j=2,dimj
          s0(j) = s0(j-1) + sqrt((x(j)-x(j-1))**2+(z(j)-z(j-1))**2)
         end do
! now compute new stretched grid
         a0 = fac*s0(dimj)/(2.0*pi)

         sn(1) = 0.0
         do j=2,dimj
           sn(j) = s0(j) + a0 * sin(2.0*pi*s0(j)/s0(dimj))
         end do
    ! force first and last point to match, as sometimes it does not
         sn(1) = s0(1)
    !     sn(2) = s0(2)
         sn(dimj) = s0(dimj)    
    !     sn(dimj-1) = s0(dimj-1)

    ! now interpolate locations
      do j=1,dimj
       xn(j) = interp1d(sn(j),s0,x,dimj)
       zn(j) = interp1d(sn(j),s0,z,dimj)
      enddo

      open(1,file=trim(datapath)//'egrid.dat',status='unknown',recl=500)
      do j=1,dimj
       write(1,*)j,s0(j),x(j),z(j),sn(j),xn(j),zn(j)
      end do

      close(1)

      x = xn
      z = zn

      return

end subroutine remap_state
