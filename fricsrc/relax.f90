
subroutine relaxdc (inputfile, iread, outputfile, iwrite, isw, ierr)

! subroutine relaxdc
! 3rd order Adams-Bashforth is used now
! modified for arbitrary grid 3/95 frt
! frictional term (frica) is now varied according to how
! the force varies b/w each iteration
! file relax.inp has some initial settings
! pr(i,j,k) is a switch which is set to zero where the radius is < rmin
! - this ensures that the force at r < rmin =0
! - this one adjusts dt if a instability is detected and
! restarts from a saved time 9/97 frt
! USES OPNE BC at back boundary

! in is now defined as the interation within each call to this subr
! in_tot is the total number if iterations

! 7/97 added computations of total energy w_tot and free energy free_tot
! 4/98 reversed order  if io do loops for speed

! solves momentum equation through jxB using mixture between
! bnew and bold
! avoids taking powers of u by setting u equal to p and solving
! implicitely for the pressure divb=0 all boundaries

! dieser Teil dient ab sofort nur noch der Berechnung
! wie fct2.f, aber mit geeignetem splitting der Glaettung
! streamlined version

!---------------relaxd parameters------------
! inputfile - filename to get initial data from
! iread - record number of inputfile to read
! outputfile - filename to write final data to
! iwrite - record number of outputfile to write
! isw > 0  to compute
! ierr - return an error value
!---------------------------------------------

    USE coredata
    IMPLICIT NONE

    character (LEN=100), intent(in) :: inputfile
    character (LEN=100), intent(in) :: outputfile
    character (LEN=100) :: name
    integer(iprec), intent(in) :: iread, iwrite, isw
    integer(iprec), intent(out) :: ierr

    real(rprec), dimension (zahlx,zahlz) :: w,wb,wp,wv,free
    real(rprec) :: rho_init(zahlx,zahlz)
    real(rprec) :: gsn(zahlx,zahlz)     ! grid-scale-noise estimate

    real(rprec) :: gsnave,gsnave_inner,gsnave_tail    ! estimates of the grid-scale-noise
    real(rprec) :: dt=0.0       ! the length of a timestep
    real(rprec) :: afmax                ! maximum ek (force magnitude)
    real(rprec) :: velmax,velave        ! max/average velocities
    real(rprec) :: ekave_old            ! old average force imbalance
    real(rprec) :: umin
!   real(rprec) :: jxb                  ! magnitude of jxB      (n/u)
!   real(rprec) :: gp                   ! magnitude of grad(P)  (n/u)
!   real(rprec) :: bb                   ! magnitude of B-field  (n/u)
    real(rprec) :: a1=0.0_rprec,a2=0.0_rprec,a3=0.0_rprec ! Adams-Bashforth multi-step weights

    real(rprec) :: w_tot,wb_tot,wp_tot,wv_tot,free_tot
    real(rprec) :: w_tot_old
    real(rprec) :: wv_tot_old
    real(rprec) :: wv_sum               ! Cumulative sum of kinetic energy ??
    real(rprec) :: rho_ave

    real(rprec) :: dmin,dmax
    real(rprec) :: dt_min,dt_max        ! dt should not exceed these bounds
    real(rprec) :: tstep                ! = dt giving CFL <= 1 everywhere

    integer(iprec) :: in,ixk,iyk,izk,ixv,iyv,izv
    integer(iprec) :: iumin,jumin,kumin
    integer(iprec) :: iterv
    integer(iprec) :: iopen

!   real(rprec) :: size         = 10.0
!   real(rprec) :: velminf      = 0.01
!   real(rprec) :: velmaxf      = 1.0
!   real(rprec) :: mu0          = 1.257e-6
!   real(rprec) :: am           = 10.0
!   real(rprec) :: om           = 0.050
    real(rprec) :: wv_sum0      = 1.0e-4
    real(rprec) :: gamma5       = 1.66666667
    integer(iprec) :: iterv0           = 10

    integer(iprec), dimension(2) :: point      ! holds a point location


    call read_relax_inp ()    ! Read friction-code-specific parameters
 

    IF (isw <= 0 .OR. isw > 2) then
      write(*,*) ' You must be using an invalid value for isw!'
      ierr = -1
      return
    else if (isw == 1) then
      iopen = 0
    else
      iopen = 1
      write (*,*) 'Warning: different boundary conditions will be used!'
      STOP 'invalid isw, stopping the code'
    end if


    !---- Begin setting initial conditions and values of intermediate arrays:

    call read_data(inputfile, iread, ierr ) ! Read initial conditions:
    if (ierr < 0) STOP 'Read error at the beginning of relaxdc'



    ! Write out initial conditions to data log file:
    name = trim(dirname)//'/'//trim(friclog)
    IF (in_tot <= 0) THEN  ! start new log data file
       call write_data(name, 1, ierr)
    ELSE                   ! Append to the existing file
       call write_data(name, 0, ierr)
    END IF
    if (ierr < 0) return


    call precalculations()   ! Calculate a few things that we'll need later



    ! Set initial values for momentum (dataset stores RHO and V):
    sx = vx*rho;  sz = vz*rho


    ! Zero initial B-Bdipole inside Earth 
    ! (normally done in push_B but needed by push_mom called first):
!    bxj = bxj*pr2;  bzj = bzj*pr2



    gamma = gammainp
    ! if the wrong gamma is used, then it could cause the code to crash
    if(gamma <= 0.0)then
            write(6,*)' WARNING: gamma le 0, aborting'
            stop
    end if
    write(*,*)' gamma used here:',gamma



    ! press_fac is in case we want to increase the pressure uniformly
    ! by some factor.  It is usually just set to 1.0

    p(:,:) = press_fac*p(:,:)
    call floor_it (zahlx,zahlz,p,pressmin)

 
    ! Floor density to avoid dividing by zero in velocity calculation:
    
        write(*,*) minval(rho)
  
    rho = MAX (rho(:,:), rhomin)
        write(*,*) minval(rho)
    
    !---- End setting initial conditions and values of intermediate arrays-----


    Call Initialize_diagnostic_files 


    CALL Calculate_grid_spacing


    ierr = 0
    wv_sum = 0.0
    iterv = 0
    iend = iend + in_tot    ! reset iend to in_tot+iend
    test = .FALSE.
    ekave_old = 0.0


    ! Set up the timestep parameters

!    dt = 0.9*dmin/vfast
     dt_max = 0.9*dmin/vfast
     dt_min = 0.0001*dt_max

    open(97,file=trim(dirname)//'/'//'time.dat',status='unknown',position='append')
       write(97,*)in_tot,zeit,dt
    close(97)


    ! INITIALIZE DENSITY:

    rho_init = rho

    if (sheathrho == 1) then
      where (rho_init == -1.0)
         rho = huge_rho
      endwhere
    endif

    if (push_rho_opt == -1) then
       ! do nothing, it will be computed at first step
    else if(push_rho_opt == 1) then    ! set rho
!      where(r> 2.0)
       rho = ((bx**2+bz**2) + gamma5*p)/vfast**2
       rho = max(rho, rhomin)
!      end where
    else
       write(*,*)' wrong value of push_rho =',push_rho_opt,' aborting'
       stop
    endif
    write(*,*)' minval(rho) = ', minval(rho)



    in = 0 ! in is an internal counter

 MAIN_LOOP: DO   ! Main computational loop

    in = in + 1
    in_tot = in_tot + 1

    ! Boundary Conditions:
    call set_bc0 (nx,nz,ssx,ssz,rx,rz,ddx,ddz, &
                  p,rho,vx,vz,sx,sz,bxj,bzj, &
                  in,zeit,dt,iopen)
    call floor_it (nx,nz,p,pressmin)
    rho = MAX (rho(:,:), rhomin)
    CALL compute_forces()
    ! find the maximum value of ek, and its location
    afmax = maxval(ek)
    point = maxloc(ek)
    ixk = point(1);  izk = point(2)
!   write(6,'(a6,g12.4,2(1x,i3),2(1x,f12.4))') 'afmax=',afmax,point,ssx(ixk),ssz(izk)
    ! compute average force imbalance (2 algorithms used)
    ekave_old = ekave
    if (force_norm == 1) then   ! force_norm = 1 uses a volume average
        ekave       = volume_avg(ek,global_mask)
        ekave_tail  = volume_avg(ek,tail_mask)
        ekave_inner = volume_avg(ek,inner_mask)
    else                        ! this one uses a grid-based average
        call average (zahlx,zahlz,ek,ekave)
    end if
    if (in == 1) ekave0 = ekave


    CALL Set_frictional_term (fric_alg, fric_min, fric_max, &
                              ekave_old, ekave, zeit, frica)


    if (mod(in,it_print_diagnostic) == 0 .OR. in <= 2) then  ! diagnostic values

        CALL Compute_energy (zahlx,zahlz,wb,wp,wv,w,free)

        ! find the maximum velocity, and its location
        velmax = maxval(vel); point = maxloc(vel)
        ixv = point(1); izv = point(2)

        call average (zahlx,zahlz,vel,velave) ! average v (velave)

        rho_ave = volume_avg(rho, global_mask)

        w_tot_old  = w_tot
        wv_tot_old = wv_tot
            
        wb_tot = volume_int(wb, global_mask)
        wp_tot = volume_int(wp, global_mask)
        wv_tot = volume_int(wv, global_mask)
        free_tot = volume_int(free, global_mask)

        wv_sum = wv_sum + wv_tot
        iterv = iterv + 1
        if (iterv > iterv0) then
            wv_sum=0.0
            iterv = 0
        end if

        w_tot = Volume_int (w, global_mask)

        CALL Compute_grid_noise (zahlx, zahlz, gsn, &
                                 gsnave, gsnave_tail, gsnave_inner)
       
        ! estimate div(B). Subroutine will put zeros in ghostcells:
!        call diver(nx,ny,nz,bxj,byj,bzj,rx,ry,rz,divb)
!       divb = divb*pr2
        where (r < 2)
          divb=0.
        endwhere
        divb_sum = volume_int (divb, global_mask)


        CALL Print_diagnostics ()
         
    end if



    ! AB3 parameter and Courant condition setup (rho calculation)

    if (in == 1) then
      a1 = 1.0_rprec; a2 = 0.0_rprec; a3 = 0.0_rprec
    else if (in == 2) then
      a1 = 1.5_rprec; a2 = -0.5_rprec; a3 = 0.0_rprec
    else if (in == 3) then
      a1 = 23.0_rprec/12.0_rprec; a2 = -16.0_rprec/12.0_rprec; a3 = 5.0_rprec/12.0_rprec
    end if
    


    ! Determine the timestep based on the Courant condition

    if (push_rho_opt == -1) then
       dt = 1.0_rprec
    else
       call courant (zahlx,zahlz,ddx,ddz,vx,vz, &
                     bx,bz,p,rho,gamma, tstep, cfl)
       dt = tstep
    end if

!   write(6,*) 'dt=',dt 
   if (mod(in,it_print_diagnostic) == 0) then
        write(*,'(A,F8.5,A,F8.5)') ' Timestep = ', dt, ',  CFL = ', cfl
        open(97,file=trim(dirname)//'/'//'time.dat',status='unknown',position='append')
        write(97,*)in_tot,zeit,dt, timeut
        close(97)
    end if

    if(dt < dt_min)then
        write(6,*)' dt too small, aborting..., dtmin = ', dt_min, 'dt = ', dt
        write(6,*)' writing to ', trim(dirname)//'/'//trim(friccrash)
        call write_data(trim(dirname)//'/'//trim(friccrash), 1, ierr)
        STOP 'code became unstable, stopping'
    end if
    zeit = zeit + dt
    ! save solution from the previous iteration:
    bxjo = bxj;  bzjo = bzj
    bxo = bx;  bzo = bz
    po = p
    rhoo = rho
    sxo = sx;  szo = sz
    vxo = vx;  vzo = vz
    dsx_2 = dsx_1; dsx_1 = dsx; dsx = kx
    dsz_2 = dsz_1; dsz_1 = dsz; dsz = kz

    ! advance momentum:
    call push_mom(a1, a2, a3, frica, dt)
    IF (limiter_choice == 1) then
       call limit3(sx,rx,rz,real(.08,rprec))
       call limit3(sz,rx,rz,real(.08,rprec))
    ELSE if (limiter_choice == 2) then
       ! make sure you set the bc below the equator before 
       ! and after calling the smoother     
       sx(:,1) =  sx(:,3);  sz(:,1) = -sz(:,3)     
       call smoothf(zahlx,zahlz,ssx,ssz,sx,real(0.25,rprec))
       call smoothf(zahlx,zahlz,ssx,ssz,sz,real(0.25,rprec))
       sx(:,1) =  sx(:,3);  sz(:,1) = -sz(:,3)     
    END IF


    ! Advance density:
    if (push_rho_opt == 1 )then
       call push_rho (a1, a2, a3, dt)
    else 
       call Set_rho_lemon (zahlx,zahlz,ddx,ddz, &
                           sxo,szo,bxo,bzo,po,rho,gamma,dt,cfl)
    end if
    rho = MAX (rho, rhomin)
    CALL Compute_velocity (zahlx,zahlz,sx,sz,rho,vx,vz)
    call push_B(a1, a2, a3, dt)
    call push_p(a1, a2, a3, dt)
    call floor_it(nx,nz,p,pressmin)

   
    ! Now the investigation of the residual force

    if(zeit >= tgamma .AND. gamma > 2.0)then
        gamma = 5./3.
        write(6,*)' .........gamma set to 5/3.......'
    end if

    ! if Wv (kinetic energy) is less than it was it_print_diagnostic iterations ago, then
    ! it has peaked.  Set v = 0 to remove energy from the system.    
    if(stopsign .AND. wv_tot_old >= wv_tot .AND. wv_tot > 0.0)then

        write(6,*)' setting v = 0'
        vx = 0.0;  vz = 0.0
        sx = 0.0;  sz = 0.0
        wv_tot = 0.0

    endif



    ! DECIDE WHETHER TO CONTINUE, STOP, OR PRINT AND PLOT:

    if((mod(in,in_write) == 0 ) .AND. in_write > 0) then  ! write log file
        name = trim(dirname)//'/'//trim(friclog)
        call write_data(name, 0, ierr)
        if (ierr < 0) return
    end if


    if ((in_tot >= iend) .OR. &
        ((wv_sum > 0.) .AND. (wv_sum <= wv_sum0) .AND. (iterv == iterv0)) .OR. &
        (ekave <= ekmin .AND. ekave > TINY(1.0_rprec)) .AND. mod(in,2) == 0) then

        IF (in_tot >= iend) THEN
           WRITE (6,*) ' Code stopped, max number of iterations reached'
        ELSE IF ((wv_sum > 0 .AND. (wv_sum <= wv_sum0))) THEN
           write (6,*) ' Code stopped, max kinetic energy =', wv_sum,'<', wv_sum0
        ELSE 
           write (6,'(//A/)') ' Code stopped, convergence crit met '
        END IF

        call write_data (outputfile, iwrite, ierr)
        if (ierr < 0) STOP 'Write error at the end of relaxdc'

        call diag(bzj,vx,rho,p,ssx,zeit,in_tot)

        write(6,*)' in =',in_tot,' zeit =',zeit,' ekave =',ekave
        write(6,*)' iend =',iend,' tend =',zeit,' ekmin =',ekmin

        if (mod(in,in_write) /= 0 )then
            call write_data(trim(dirname)//'/'//trim(friclog), 0, ierr)
            if (ierr < 0) return
        endif

        EXIT Main_Loop

    else     ! check  to see if the code has blown up

        if (velave > 1.0e2 .OR. ekave > 1.0e2 .OR. velmax > 1.0e4) then
        !   if(dt > dt_min) then
                write(6,'(/A/)')' relaxd code went unstable, aborting'
                write (6,*) 'velave=',velave, ' 1.e2'
                write (6,*) 'ekave=', ekave, '1.0E2'
                write (6,*) 'velmax=', velmax, '1.0E3'
             name = trim(dirname)//'/'//trim(friclog)
             call write_data (name, 0, ierr)
             name = trim(dirname)//'/'//trim(friccrash)
             call write_data (name, 1, ierr)
             STOP
        !   end if
        end if

    endif
 

 END DO Main_Loop 


    return

    CONTAINS

    SUBROUTINE Initialize_diagnostic_files ()

    if (iread > 2) THEN  ! obviously a restart from a data file
       restart = .TRUE.
    else                  ! starting from an initial conditions file made with setup
       restart = .FALSE.
    end if

    if (.NOT. restart) then ! Wipe out old diagnostic files that might exist
                            ! and write out headers
       open(7,file=trim(dirname)//'/'//'relaxd.dat',status='unknown',recl=132)
         WRITE (7,*)'TITLE="FRIC: Monitored Parameters"'
         WRITE (7,*)'VARIABLES = "zeit","ekave","ekave_tail","ekave_inner","velave","in","afmax","fric","velmax", "ixv", &
                     & "izv"' !"
       close(7)

       open(8,file=trim(dirname)//'/'//'energy.dat',status='unknown',recl=132)
        write(8,*) 'in_tot zeit w_tot wb_tot wp_tot wv_tot free_tot rho_ave'
       close(8)

       open(20,file=trim(dirname)//'/'//'diag.txt',form='FORMATTED',status='unknown')
       open(99,file=trim(dirname)//'/'//'pressmin.txt',status='unknown',form='formatted')

       open(97,file=trim(dirname)//'/'//'time.dat',status='unknown')
        write(97,*)'      in       zeit       dt'
       close(97)

       OPEN (98, FILE=trim(dirname)//'/'//'divb.dat', STATUS='REPLACE'); CLOSE (98, STATUS='DELETE')

    end if
    RETURN
    END SUBROUTINE Initialize_diagnostic_files




    SUBROUTINE Print_diagnostics

        open(8,file=trim(dirname)//'/'//'energy.dat',status='unknown', position='append')
        write(8,124)in_tot,zeit,w_tot,wb_tot,wp_tot,wv_tot,free_tot,rho_ave
 124    format(1x,i6,6(1x,g13.6),1x,g16.8)

        close(8)

        call find_min(zahlx,zahlz,p,umin,iumin,kumin)
        write(99,*)in,umin,iumin,kumin

        write(6,1234) in_tot, zeit, afmax, ixk, izk,frica
        write(6,'(1X,a20,ES12.4,a10,2(i3,1x),a16,ES11.3)') &
       'Max v = ',velmax,' at ik = ',ixv,izv,' Average v = ',velave
        write(6,'(1X,a9,ES9.2,a5,2(a2,i3,a2,g12.3),a7,ES11.2)') &
           ' Max v = ',velmax,' at: ','x(',ixv,')=',ssx(ixv),&
                                      'z(',izv,')=',ssz(izv),' <v> = ',velave
        write(6,'(1X,A26,F9.5,A40,ES12.5)') 'Global force imbalance = ',&
            ekave,' Minimum pressure = ',umin
        write(6,'(1X,A26,F9.5,A40,F9.5)') 'Tail force imbalance = ',&
            ekave_tail,'Inner force imbalance = ',ekave_inner
        write(6,'(1X,A26,F9.5,A40,F9.5)') 'Tail_gsn = ',gsnave_tail,'Inner_gsn = ',gsnave_inner
        write(6,'(1X,A26,F9.5)') 'Global_gsn = ',gsnave

        open(7,file=trim(dirname)//'/'//'relaxd.dat',status='unknown',recl=132,position='append')
        write(7,1235)zeit,ekave, ekave_tail, ekave_inner,velave,in_tot,&
                         afmax,frica,velmax,ixv,izv
        close(7)

        write (6,'(2(A,ES12.4))') ' max(div(B)) = ', maxval(divb),&
                 '   integral(div(B)) = ', divb_sum

        write(6,'(A/)') &
        &'--------------------------Written to relaxd.dat---------------------------'
 1234   format(1x,i6,f10.3,1x,e18.8,1x,2(1x,i3),e18.8,1x,g12.4)
 1235   format(1x,f8.1,4(1x,g12.4),1x,i6,1x,g12.4,1x,g12.4,e15.8,2(1x,i3))

        RETURN
    END SUBROUTINE Print_diagnostics



    SUBROUTINE Calculate_grid_spacing

    ! find the minimum and maximum grid spacing
    dmin = max(maxval(rx),  maxval(rz))
    dmax = min(minval(rx),  minval(rz))

    write(6,*)' 1.0/maxval(rx) = ', 1.0/maxval(rx)
    write(6,*)' 1.0/maxval(rz) = ', 1.0/maxval(rz)

    if(dmin > 100*TINY(1.0_rprec))then
        dmin = 1.0/dmin
        write(6,*)' dmin =',dmin
    else
        write(6,*)' dmin undefined, aborting... '
        ierr = -2
        return
    end if

    write(6,*) ' min(minval(rx), minval(rz)) = ', &
               min(minval(rx),  minval(rz))
    if(dmax > 100*TINY(1.0_rprec))then
        dmax = 1.0/dmax
        write(6,*)' dmax =',dmax
    else
        write(6,*)' dmax undefined, aborting... '
        ierr = -2
        return
    end if

    RETURN
    END SUBROUTINE Calculate_grid_spacing



    end subroutine relaxdc
