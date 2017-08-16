MODULE coredata

IMPLICIT NONE
SAVE
    character (LEN=100) :: dirname
    INTEGER, PARAMETER :: iprec = 4 !SELECTED_INT_KIND (9)
    INTEGER, PARAMETER :: rprec = 4 !SELECTED_REAL_KIND (6,37)

! global parameters (data file names, etc.)
    character (LEN=100), PARAMETER :: friclog =  'fric2d.log'
    character (LEN=100), PARAMETER :: friccrash = 'fric2d.crash'
    character (LEN=100), PARAMETER :: relaxinp = 'relax.inp'

    REAL(rprec), PARAMETER :: pi=3.14159265
    REAL(rprec), PARAMETER :: huge_rho = 1.0e15

    real(rprec) :: rmin = 1.1
    real(rprec) :: rmin2 = 1.1
    
! input parameters
    LOGICAL :: stopsign                 ! set v=0 when KE peaks?
    LOGICAL :: restart                  ! Is this a friction code restart?
    integer(iprec) :: it_print_diagnostic  ! write data to stdout and diagnostic files every N iterations
    INTEGER(iprec) :: iend
    INTEGER(iprec) :: fric_alg, force_norm, in_write
    REAL(rprec) :: frica,ekmin,press_fac,fric_min,fric_max
    REAL(rprec) :: cfl                  ! the factor to reduce cour_min
    real(rprec) :: tgamma               ! time to switch to gamma =5/3
    REAL(rprec) :: gamma,gammainp,rhomin,pressmin
    REAL(rprec),PARAMETER :: pfactor = 5.22e-9
    REAL(rprec),PARAMETER :: bfactor = 81.0e-9
    REAL(rprec) :: pressmin0

    INTEGER :: zahlx, zahly, zahlz
    !INCLUDE 'param_fric.h'              ! parameter zahlx = ??, zahly = ??, zahlz = ??

! grid variables
    INTEGER(iprec) :: nx
    INTEGER(iprec) :: ny
    INTEGER(iprec) :: nz
    INTEGER(iprec) :: nx1,nx2,nx3,nx4
    INTEGER(iprec) :: ny1,ny2,ny3,ny4
    INTEGER(iprec) :: nz1,nz2,nz3,nz4
    REAL(rprec),allocatable ::  ssx(:)                  ! grid x-values
    REAL(rprec),allocatable ::  ssz(:)                  ! grid z-values
    REAL(rprec),allocatable ::  rx(:)                   ! central diff d/dx values
    REAL(rprec),allocatable ::  rz(:)                   ! central diff d/dz values
    REAL(rprec),allocatable :: r(:,:)
    REAL(rprec),allocatable :: pr(:,:)
    REAL(rprec),allocatable :: pr2(:,:)
    real(rprec),allocatable :: vol(:,:)       ! volume around each grid point
    real(rprec),allocatable :: xxx(:,:)       ! x values at grid points
    real(rprec),allocatable :: zzz(:,:)       ! z values at grid points
    real(rprec),allocatable :: ddx(:,:)       ! rx values at grid points
    real(rprec),allocatable :: ddz(:,:)       ! rz values at grid points

    logical,allocatable :: global_mask(:,:)
    logical,allocatable :: tail_mask(:,:)
    logical,allocatable :: inner_mask(:,:)

! density variable
    REAL(rprec),allocatable :: rho(:,:)   ! density
    REAL(rprec),allocatable :: rhoo(:,:)  ! density
    real(rprec),allocatable :: divrv(:,:)
    real(rprec),allocatable :: divrv_1(:,:)
    real(rprec),allocatable :: divrv_2(:,:)

! pressure variable
    REAL(rprec),allocatable :: p(:,:)     ! pressure
    REAL(rprec),allocatable :: po(:,:)     ! pressure
    REAL(rprec),allocatable :: dpv(:,:)
    REAL(rprec),allocatable :: divv(:,:)
    real(rprec),allocatable :: dp(:,:)   ! change in u (timestep)
    real(rprec),allocatable :: dp_1(:,:) ! change in u (timestep)
    real(rprec),allocatable :: dp_2(:,:) ! change in u (timestep)

! velocity and momentum variables
    REAL(rprec),allocatable :: sx(:,:)    ! x-momentum
    REAL(rprec),allocatable :: sz(:,:)    ! z-momentum
    REAL(rprec),allocatable :: sxo(:,:)   ! x-momentum
    REAL(rprec),allocatable :: szo(:,:)   ! z-momentum
    REAL(rprec),allocatable :: vx(:,:)    ! x-velocity
    REAL(rprec),allocatable :: vz(:,:)    ! z-velocity
    REAL(rprec),allocatable :: vxo(:,:)   ! x-velocity
    REAL(rprec),allocatable :: vzo(:,:)   ! z-velocity
    real(rprec),allocatable :: dsx(:,:)   ! x-force
    real(rprec),allocatable :: dsz(:,:)   ! z-force
    real(rprec),allocatable :: dsx_1(:,:) ! x-force at previous timestep
    real(rprec),allocatable :: dsz_1(:,:) ! z-force at previous timestep
    real(rprec),allocatable :: dsx_2(:,:) ! x-force two timesteps ago
    real(rprec),allocatable :: dsz_2(:,:) ! z-force two timesteps ago
    REAL(rprec),allocatable :: vel(:,:)   ! velocity magnitude

! magnetic field variables
    REAL(rprec),allocatable :: bx(:,:)    ! x-magnetic field
    REAL(rprec),allocatable :: bz(:,:)    ! z-magnetic field
    REAL(rprec),allocatable :: bxj(:,:)   ! Bx minus dipole
    REAL(rprec),allocatable :: bzj(:,:)   ! Bz minus dipole
    REAL(rprec),allocatable :: bxd(:,:)   ! dipolar Bx
    REAL(rprec),allocatable :: bzd(:,:)   ! dipolar Bz
    real(rprec),allocatable :: dbx(:,:)   ! change in Bx (timestep)
    real(rprec),allocatable :: dbz(:,:)   ! change in Bz (timestep)
    real(rprec),allocatable :: dbx_1(:,:) ! change in Bx (timestep)
    real(rprec),allocatable :: dbz_1(:,:) ! change in Bz (timestep)
    real(rprec),allocatable :: dbx_2(:,:) ! change in Bx (timestep)
    real(rprec),allocatable :: dbz_2(:,:) ! change in Bz (timestep)
    real(rprec),allocatable :: bxo(:,:)   ! old Bx (previous iteration)
    real(rprec),allocatable :: bzo(:,:)   ! old Bz
    real(rprec),allocatable :: bxjo(:,:)  ! old Bx (previous iteration)
    real(rprec),allocatable :: bzjo(:,:)  ! old Bz
    real(rprec),allocatable :: divb(:,:)  ! divergence of B
    real(rprec) :: divb_sum

    REAL(rprec),allocatable :: ajx(:,:)   ! Jx current
    REAL(rprec),allocatable :: ajy(:,:)   ! Jy current
    REAL(rprec),allocatable :: ajz(:,:)   ! Jz current

! Force variables
    REAL(rprec),allocatable :: kx(:,:)    ! x-comp. of force JxB-grad(P)
    REAL(rprec),allocatable :: ky(:,:)    ! y-comp. of force JxB-grad(P)
    REAL(rprec),allocatable :: kz(:,:)    ! z-comp. of force JxB-grad(P)
    real(rprec),allocatable :: ek(:,:)    ! force magnitude (L2 norm)    
    real(rprec) :: ekave                        ! average force imbalance
    real(rprec) :: ekave0                       ! force normalization factor
    real(rprec) :: ekave_tail                   ! average force imbalance in the tail
    real(rprec) :: ekave_inner                  ! average force imbalance in the inner mag

    INTEGER(iprec) :: in_tot                    ! total number of iterations
    INTEGER(iprec) :: sheathrho                 ! give sheath a high density? Input parameter
    INTEGER(iprec) :: push_rho_opt                  ! switch to either push (1) or fix (-1) rho
    INTEGER(iprec) :: pvg_corr_inter            ! Interval to correct pvg (-ve is off) 

    REAL(rprec) :: zeit                         ! time value (with german precision)
    REAL(rprec) :: vis = 0.0                    ! "viscosity" for smoother
    REAL(rprec) :: vis_mom = 0.0                ! "viscosity" for push_mom (experimental, not used routinely)
    REAL(rprec) :: dm, tilt_angle               ! dipole moment and tilt angle, read in from dataset
    REAL(rprec) :: timeut   = -1.0_rprec        ! UT magnetospheric time in seconds 

    INTEGER (iprec) :: limiter_choice = 0       ! Type of limiter to use for momentum: 0--none, 1--limit3, 2--smoothf

    LOGICAL :: use_inner_bc = .FALSE.           ! Flag indicating whether to use inner b.c.
    INTEGER (iprec) :: inner_bnd_index(3,2)     ! array holding locations (indices) of the "inner boundary"
                                                ! 1st index is coord (1-3 is x-z), 2nd index is min/max in
    LOGICAL,allocatable :: ibc_inside_mask(:,:) ! true if grid points are completely inside inner boundary
    LOGICAL,allocatable :: ibc_here_mask(:,:) ! true if grid points are on the inner boundary

    logical :: test
! ------------------pv^gamma tests- added 11/05 frt
     REAL(rprec), PARAMETER :: rs = 1.5     ! location where values are stored
     REAL(rprec), PARAMETER :: vfast = 0.9
     logical :: symmetric

     ! Interfaces:

     INTERFACE
       SUBROUTINE Dipole_modified (dm, tilt, x,z,bx,bz)
          IMPLICIT NONE
          REAL, INTENT (IN) :: dm, tilt, x,z
          REAL, INTENT (OUT) :: bx,bz
       END SUBROUTINE Dipole_modified
     END INTERFACE



CONTAINS

    subroutine read_relax_inp()
      ! subroutine to read input parameters from relax.inp
        use fricoptions
       implicit none

        frica = fric_alpha
        write(6,*)' frica (alpha friction parameter) = ',frica

        in_write = fric_interval
        write(6,*)' Data output interval (in_write)  = ',in_write

        iend = fric_max_steps
        write(6,*)' i ending                         = ',iend

        ekmin = fric_min_avg_force
        write(6,*)' minimum avg. force (ekmin)       = ',ekmin

        it_print_diagnostic = fric_print_diagnostic
        write(6,*)' diagnostic print/write interval  = ',it_print_diagnostic

        press_fac = fric_pressure_factor
        write(6,*)' press factor                     = ',press_fac

        fric_alg = fric_algorithm
        write(6,*)' friction algorithm               = ',fric_alg

        fric_min = fric_min_factor
        fric_max = fric_max_factor
        write(6,*)' min fric factor                  = ',fric_min
        write(6,*)' max fric factor                  = ',fric_max

        force_norm = fric_force_norm
        if(force_norm == 1)then
            write(6,*)' volume based force averaging scheme used'
        else
            write(6,*)' grid based force averaging scheme used'
        end if

        cfl = fric_cfl
        write(6,*)' timestep reduction factor        = ',cfl

        pressmin0 = fric_pressmin
        pressmin = pressmin0/pfactor
        write(6,*)' press min                        = ',pressmin

        rhomin = fric_rhomin
        write(6,*)' rhomin                           = ',rhomin

        gammainp = fric_gammainp
        write(6,*)' Initial gamma                    = ',gammainp

        tgamma = fric_tgamma
        write(6,*)' Time to switch gamma = 5/3       = ',tgamma

        stopsign = fric_stopsign
        write(6,*)' Set v = 0 when K.E. peaks?         ',stopsign

        sheathrho = fric_sheathrho
        write(6,*)' Give sheath a high density?      = ',sheathrho

        push_rho_opt = fric_push_rho_opt
        write(6,*)' push_rho_opt                     = ',push_rho_opt
        if(push_rho_opt == -1)then
        write(*,*)'                                      (rho to be set)'
        elseif(push_rho_opt == 1)then
        write(*,*)'                                      ( rho to be pushed)'
        else
        write(*,*)' incorrect value of push_rho, aborting'
        stop
        endif

        pvg_corr_inter = fric_pvg_corr_inter
        write(6,*)' Interval to correct  pv^gamma    = ',pvg_corr_inter

        limiter_choice = fric_limiter
        IF (limiter_choice == 0) then
           write(6,*) ' Limiter choice:              =  NONE'
        ELSE if (limiter_choice == 1) then
           write(6,*) ' Limiter choice:              =  LIMIT3'
        ELSE if (limiter_choice == 2) then
           write(6,*) ' Limiter choice:              =  SMOOTHF'
        ELSE
           STOP 'invalid value for limiter_choice in relax.inp, stopping'
        END IF

         vis_mom = fric_vis_mom
        write(6,*)' Viscosity for momentum eq. (test,unused) = ',vis_mom

         use_inner_bc = fric_use_inner_bc
        write(6,*) ' Use inner boundary cond?      = ', use_inner_bc
        write(*,*)
        write(*,*) '-----------Finished reading input parameters--------------'
        write(*,*)

    end subroutine read_relax_inp

    

    SUBROUTINE write_data(filename, recnum, ierr)

       ! If recnum == 0, append to the file
       ! If recnum == 1, replace file
       ! If recnum > 0, then write to that record
 

       character (LEN=100) :: filename
       integer, intent(in) :: recnum
       integer, intent(out) :: ierr
       integer :: i, in_temp
       integer, parameter :: LUN_write=91
       LOGICAL :: LogicalFlag

       ierr = 0

       INQUIRE (UNIT=LUN_write, OPENED=LogicalFlag)
       IF (LogicalFlag) STOP 'ERROR IN WRITE_DATA, UNIT IS ALREADY OPEN, STOPPING'


       IF (recnum == 1) then

          write(*,991) '--> Write to '//trim(filename)//' [REPLACE] , REC=', recnum, ' IT=',in_tot
          open (LUN_write,file=trim(filename),form='UNFORMATTED',status='replace')

       ELSE if (recnum == 0) then

          write(*,991) '--> Write to '//trim(filename)//' [APPEND ] , REC=', recnum, ' IT=',in_tot
          open (LUN_write,file=trim(filename),form='unformatted',status='old',&
                position='append')

       ELSE if (recnum > 1) then

          write(*,991) '--> Write to '//trim(filename)//' [POSITION], REC=', recnum, ' IT=',in_tot

          open (LUN_write,file=trim(filename),form='UNFORMATTED',status='old')
          do i=1,recnum-1
              read (LUN_write,end=2001,err=2001) in_temp
              write (6,'(A,I6.6,A,I6.6)') ' WRITE_DATA, POSITIONING, i = ', i, ' REC =',in_temp
          end do

          goto 2003
 2001      ierr = -1
           write(6,*)' Can''t read record number ', recnum, ' from ',filename
           STOP
 2003     continue

       ELSE
          STOP ' IN SUBROUTINE WRITE_DATA, RECORD NUMBER IS NEGATIVE, STOPPING'
       END IF
!        write(*,*) 'writedata',bz(30,5),bz(30,5)*bfactor

       ! Here we actually write the data
       write (LUN_write) in_tot,zeit, timeut, nx,nz,dm,tilt_angle,ssx,ssz,rx,rz, &
                         rho,vx,vz,bx*bfactor,bz*bfactor,p*pfactor, &
                         vis,gamma,kx,kz,ekave0,ekave

       CLOSE (LUN_write)

991 FORMAT (//TR2,A,I5.5,A,I8.8)
    RETURN
    END SUBROUTINE write_data



    subroutine read_data(filename, recnum, ierr)

    ! recnum is the record number to read


        character (LEN=100) :: filename
        integer(iprec), intent(in) :: recnum
        integer(iprec), intent(out) :: ierr
        integer(iprec) :: i, in_temp, ix,iz, istat
        integer(iprec), parameter :: LUN_read = 42
        real(rprec) :: bd1,bd2,bd3

        write (*,'(///A,I6.6)') ' READ_DATA: Opening file='//trim(filename)//' for READING REC=',recnum
        open (lun_read,file=trim(filename),form='unformatted',status='old')
        ierr = 0


        do i=1,recnum-1
            read (lun_read,end=151,err=151,iostat=istat)in_temp
            write(6,*)' i = ', i, ' in =',in_temp
        end do

        read (lun_read,end=151,err=151,iostat=istat) in_tot,zeit, timeut, nx,nz,dm,tilt_angle,ssx,ssz,rx,rz, &
                                       rho,vx,vz,bx,bz,p,vis,gamma,kx,kz,ekave0,ekave
        write(6,'(4(A,I5),A,ES12.5,A,ES12.5)')' i= ',i,' in =',in_tot,' nx= ',nx,' nz = ',nz, &
           ' time =',zeit,' gamma =',gamma

        goto 152
            151 close(lun_read)
            write(6,*)' Can''t read record number ', recnum, ' from ',filename, 'iostat=',istat
            ierr = -1
            return
        152 close(lun_read)

        write(*,*) 'Relaxd  nx,nz = ', nx,nz
        if(nx == 0 .OR. nz == 0)then
            ierr = -1
            return
        end if
        if(nx /= zahlx .OR.  nz /= zahlz)then
            ierr=-1
            write(*,*)'Dimension error on relax,  ', filename
            return
        end if

        ! convert to code normalized units 
        bx = bx/bfactor ;  bz = bz/bfactor
        p = p/pfactor

        do iz=1,nz
        do ix=1,nx
   
          CALL Dipole_modified (dm,tilt_angle,ssx(ix),ssz(iz),bd1,bd3)
  
          ! store the dipole field
          bxd(ix,iz)=  bd1/(bfactor*1.0e9)
          bzd(ix,iz)=  bd3/(bfactor*1.0e9)

        end do
        end do


        ! subtract the dipole
        bxj = bx! - bxd   no dipole case
        bzj = bz! - bzd

    return
    end subroutine read_data





    SUBROUTINE Push_rho (a1,a2,a3,dt)
    IMPLICIT NONE
    integer(iprec) :: ix,iz
    real(rprec), intent(in) :: a1,a2,a3,dt

    divrv_2 = divrv_1
    divrv_1 = divrv

    ! calculate div( rho*v )
    forall (ix=2:nx1,iz=2:nz1)
        divrv(ix,iz) = rx(ix)*(sxo(ix+1,iz) - sxo(ix-1,iz)) &
                        + rz(iz)*(szo(ix,iz+1) - szo(ix,iz-1))
    end forall

    ! update rho
    rho(2:nx1,2:nz1) = rho(2:nx1,2:nz1) &
                      - dt * (a1*divrv  (2:nx1,2:nz1)  + &
                              a2*divrv_1(2:nx1,2:nz1)  + &
                              a3*divrv_2(2:nx1,2:nz1))
    RETURN
    END SUBROUTINE Push_rho







subroutine push_B(a1,a2,a3,dt)
IMPLICIT NONE
integer(iprec) :: ix,iz
real(rprec), intent(in) :: a1,a2,a3,dt
real(rprec) :: fyb(zahlx,zahlz) ! scratch grid variable

! Calculate -E using Ohm's law for a perfect conductor
! (fx,fy,fz) = v x B
    fyb(:,:) = vzo(:,:)*bxo(:,:) - vxo(:,:)*bzo(:,:)
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Bx
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!     write(*,*) 'pushB',maxval(bz)
    dbx_2 = dbx_1
    dbx_1 = dbx

    forall (ix=2:nx1,iz=2:nz1)
        dbx(ix,iz) = -rz(iz)*(fyb(ix,iz+1) - fyb(ix,iz-1)) 
    end forall
               
    bxj = bxjo + dt*(a1*dbx + a2*dbx_1 + a3*dbx_2)

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Bz
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    dbz_2 = dbz_1
    dbz_1 = dbz

    forall (ix=2:nx1,iz=2:nz1)
        dbz(ix,iz) = (fyb(ix+1,iz) - fyb(ix-1,iz))*rx(ix) 
    end forall

    bzj = bzjo + dt*(a1*dbz + a2*dbz_1 + a3*dbz_2)

!    bxj = bxj*pr2;  bzj = bzj*pr2
     
    bx = bxj !+ bxd ! no dipole
    bz = bzj !+ bzd
!   write(*,*) 'pushB',maxval(bz)
end subroutine push_B





subroutine push_mom(a1,a2,a3,fric_alpha_in,dt)
IMPLICIT NONE
integer(iprec) :: ix,iz
real(rprec), intent(in) :: a1,a2,a3,fric_alpha_in,dt
real(rprec) :: fxs(zahlx,zahlz)   ! scratch grid variable
real(rprec) :: fys(zahlx,zahlz)   ! scratch grid variable
real(rprec) :: fzs(zahlx,zahlz)   ! scratch grid variable


!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! push sx (x-momentum)
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

! (fx,fy,fz) = rho*vx*v
    fxs(:,:) = sxo(:,:)*vxo(:,:)
    fzs(:,:) = sxo(:,:)*vzo(:,:)
! jxB - grad(p) - div(fx,fy,fz) - alpha*rho*vx
    forall (ix=2:nx1,iz=2:nz1)
        dsx(ix,iz) = dsx(ix,iz)                              &
                      - (fxs(ix+1,iz) - fxs(ix-1,iz))*rx(ix) &
                      - (fzs(ix,iz+1) - fzs(ix,iz-1))*rz(iz) &
                      - fric_alpha_in*sxo(ix,iz)
        

    end forall
    dsx(:,:) = dsx(:,:)*pr(:,:)
    sx = sxo + dt*(a1*dsx + a2*dsx_1 + a3*dsx_2)

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! push sz (z-momentum)
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    fxs(:,:) = szo(:,:)*vxo(:,:)
    fzs(:,:) = szo(:,:)*vzo(:,:)

! jxB - grad(p) - div(fx,fy,fz) - alpha*rho*vz
    forall (ix=2:nx1,iz=2:nz1)
        dsz(ix,iz) = dsz(ix,iz)                              &
                      - (fxs(ix+1,iz) - fxs(ix-1,iz))*rx(ix) &
                      - (fzs(ix,iz+1) - fzs(ix,iz-1))*rz(iz) &
                      - fric_alpha_in*szo(ix,iz)
    end forall

    dsz(:,:) = dsz(:,:)*pr(:,:)

    sz = szo + dt*(a1*dsz + a2*dsz_1 + a3*dsz_2)

    sx(:,:) = sx(:,:)*pr(:,:)
    sz(:,:) = sz(:,:)*pr(:,:)

end subroutine push_mom

subroutine push_p(a1, a2, a3, dt)
IMPLICIT NONE
integer(iprec) :: ix,iz
real(rprec), intent(in) :: a1,a2,a3,dt
real(rprec) :: fxp(zahlx,zahlz) ! scratch grid variable
real(rprec) :: fzp(zahlx,zahlz) ! scratch grid variable


!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Push u (pressure)
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   dp_2(:,:) = dp_1(:,:)
   dp_1(:,:) = dp(:,:)

   fxp(:,:) = po(:,:)*vxo(:,:)
   fzp(:,:) = po(:,:)*vzo(:,:)

    forall (ix=2:nx1,iz=2:nz1)
        dpv(ix,iz) = (fxp(ix+1,iz) - fxp(ix-1,iz))*rx(ix) &
                      + (fzp(ix,iz+1) - fzp(ix,iz-1))*rz(iz)

        divv(ix,iz) = (vxo(ix+1,iz) - vxo(ix-1,iz))*rx(ix) &
                       + (vzo(ix,iz+1) - vzo(ix,iz-1))*rz(iz)
    end forall

    dp(:,:) = dpv(:,:) + (gamma - 1.0)*po(:,:)*divv(:,:)

    p(:,:) = po(:,:) - dt*(a1*dp + a2*dp_1 + a3*dp_2)

end subroutine push_p




subroutine precalculations()

  IMPLICIT NONE
  integer :: ix,iz

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! set up some variables
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    nx1 = nx-1
    nz1 = nz-1

    nx2 = nx-2
    nz2 = nz-2

    nx3 = nx-3
    nz3 = nz-3

    nx4 = nx-4
    nz4 = nz-4
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Calculate rx,ry, and rz (try numerically)
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    forall(ix=2:nx1)
        rx(ix) = 1.0/(ssx(ix+1) - ssx(ix-1))
    end forall
    forall(iz=2:nz1)
        rz(iz) = 1.0/(ssz(iz+1) - ssz(iz-1))
    end forall

    rx(1) = 0.5/(ssx(2) - ssx(1))
    rx(nx) = 0.5/(ssx(nx) - ssx(nx1))
    rz(1) = 0.5/(ssz(2) - ssz(1))
    rz(nz) = 0.5/(ssz(nz) - ssz(nz1))

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! set up xxx,yyy, and zzz (3-d versions of ssx,ssy,ssz)
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    forall (iz=1:nz)
            xxx(:,iz) = ssx(:)
    end forall
    forall (ix=1:nx)
            zzz(ix,:) = ssz(:)
    end forall    

    forall (iz=1:nz)
            ddx(:,iz) = rx(:)
    end forall
    forall (ix=1:nx)
            ddz(ix,:) = rz(:)
    end forall

    r(:,:) = sqrt(xxx*xxx +  zzz*zzz)

    vol(:,:) = (0.5/ddx) * (0.5/ddz)


!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! set up masks for calculating the force convergence parameters
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

     global_mask = .TRUE.
     tail_mask = .TRUE.
     inner_mask = .TRUE.

     where (r < 1.0 .OR. zzz <= ssz(2)   .OR. zzz > ssz(nz-1) &
                    .OR. xxx < ssx(2) .OR. xxx > ssx(nx-1)) 
        global_mask = .FALSE.
        tail_mask = .FALSE.
        inner_mask = .FALSE.
     endwhere

     where (xxx > -8.0 .OR. zzz > 3.0)
        tail_mask = .FALSE.
     endwhere

     where (r > 8.0 .OR. r < 2.0)
        inner_mask = .FALSE.
     endwhere

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! now setup a switch to make force = 0 for r < rmin
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    where (r > rmin)
        pr = 1.0
    elsewhere
        pr = 0.0
    endwhere
! In one version of the code, pr2 was set to 0.8, which is a
! problem as it effects the currents.  Also, the boundary
! was at 2 re, not rmin.  
! reset to 2*rmin 10/06 frt
    where (r > rmin)
        pr2 = 1.0
    elsewhere
        pr2 = 0.0
    endwhere
! for test no mask set 
!      pr=1.
!      pr2=1.
! exclude region earthward of 3 re
    !where(xxx > 1.0)
     !       pr = 0.0
    !endwhere
end subroutine precalculations




subroutine compute_forces()

    IMPLICIT NONE

    real(rprec) :: lx(zahlx,zahlz)        ! -(jxB)x force
    real(rprec) :: lz(zahlx,zahlz)        ! -(jxB)z force

    real(rprec) :: dpx(zahlx,zahlz)       ! x-pressure gradient
    real(rprec) :: dpz(zahlx,zahlz)       ! z-pressure gradient

    integer(iprec) :: ix,iz
  

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Calculate current before we change stuff (B, in particular)
! Compute the curl of B to find the current (Ampere's Law)
! J = (ajx, ajy, ajz) = curl(Bj)  (where Bj is the field minus dipole)
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= &

 forall (ix=2:nx1,iz=2:nz1)
        ajy(ix,iz) = -rx(ix)*(bzj(ix+1,iz) - bzj(ix-1,iz))   &
                      +  rz(iz)*(bxj(ix,iz+1) - bxj(ix,iz-1))
  end forall
! (lx,ly,lz) = jxB = (ajx,ajy,ajz) x (bx,by,bz)
    forall (ix=2:nx1,iz=2:nz1)
     lx(ix,iz) = ajy(ix,iz)*bz(ix,iz)  
     lz(ix,iz) = - ajy(ix,iz)*bx(ix,iz)
    end forall

! (dpx,dpy,dpz) = grad(P) = pressure gradient
    forall (ix=2:nx1,iz=2:nz1)
       dpx(ix,iz) = -(p(ix+1,iz) - p(ix-1,iz))*rx(ix)
       dpz(ix,iz) = -(p(ix,iz+1) - p(ix,iz-1))*rz(iz)
    end forall
! (kx,ky,kz) = jxB - grad(P) = force   
    do ix=2, nx1;do iz=2,nz1
        kx(ix,iz) = dpx(ix,iz) + lx(ix,iz)
        kz(ix,iz) = dpz(ix,iz) + lz(ix,iz)

!    write(*,*) 'ix,iz',ix,iz,'lx,lz',lx(ix,iz),lz(ix,iz)
!    write(*,*) 'dpx,dpz',dpx(ix,iz),dpz(ix,iz),p(ix,iz+1),p(ix,iz-1)
!    write(*,*)  'kx,kz',kx(ix,iz),kz(ix,iz)
    end do ;end do
      

    forall (ix=3:nx2, iz=2:nz2)
    ! calculate the force magnitude
        ek(ix,iz) = sqrt( kx(ix,iz)*kx(ix,iz)  &
          + kz(ix,iz)*kz(ix,iz) )*pr(ix,iz)
    end forall

end subroutine compute_forces





subroutine trapezoid(imax,x,y,z)
! simple trapezoidal integration routine
! 11/06 frt
IMPLICIT NONE
     integer(iprec),intent(in) :: imax
     real(rprec),intent(in) :: x(imax),y(imax) ! x and y values of integrand
     real(rprec),intent(out) :: z  ! integral
     real(rprec) :: dx(imax)
     integer :: i

     z = 0.0

     do i=1,imax-1
      dx(i) = x(i+1) - x(i)
      if(dx(i) < 0.0)then
              write(6,*)' error in trapezoid, dx <=0 at x =',x(i)
              write(6,*)' imax =',imax,' dx(',i,')=',dx(i)
              return
      endif
     end do

      do i=1,imax-1
        z = z + dx(i)*(y(i)+y(i+1))/2.0
      end do

      return
end subroutine trapezoid

real function volume_avg(f,mask)

   real(rprec) :: f(zahlx,zahlz)
   logical :: mask(zahlx,zahlz)
   real(rprec) :: totalvol(zahlx,zahlz)
   real(rprec) :: total(zahlx,zahlz)

   where (mask)
     total = f
     totalvol = vol
   elsewhere
     total = 0
     totalvol = 0
   endwhere
   
   if (sum(totalvol) /= 0) then
     volume_avg = sum(total*totalvol)/sum(totalvol)
   else
     volume_avg = 9999.0
   end if

end function volume_avg

real function volume_int(f,mask)

   real(rprec) :: f(zahlx,zahlz)
   logical :: mask(zahlx,zahlz)
   real(rprec) :: totalvol(zahlx,zahlz)
   real(rprec) :: total(zahlx,zahlz)

   where (mask)
     total = f
   elsewhere
     total = 0
   endwhere

   volume_int = sum(total)

end function volume_int


END MODULE coredata

!-----------------------------------------------------

      subroutine tracf2(xa,ya,za,fa,xn,yn,zn,fn,er,ep,h,ierr)

!================================================================
! purpose:
! field line tracing program, based on the Voigt first order tracer
! keeps track of changes in flux tube volume to adjust stepsizes
! f90 version 10/03 frt
!
! This routine ensures that your errors will never be greater than
! the "er" error tolerance.  It takes only one step along the
! field line (or against, depending on "ep"), so this routine will
! presumably be called many times

! inputs:
!     xa,ya,za: starting location
!           fa: starting flux tube volume value
!           er: error tolerance for tracer
!           ep: direction (along =1, against -1) along field line
!            h: stepsize estimate (new value will be outputted)
! outputs:
!     xf,yf,zf: ending location
!           vf: end flux tube volume value    
!            h: new stepsize
!         ierr: error flag
!

      
      USE coredata, only: iprec,rprec
      implicit none
      real(rprec) :: xa,ya,za,fa,xn,yn,zn,fn
      real(rprec) :: dx1,dy1,dz1,bf1,h,ep,er
      real(rprec) :: dxt,dyt,dzt,bft
      real(rprec) :: xn1,yn1,zn1,fn1
      real(rprec) :: xn2,yn2,zn2,fn2
      real(rprec) :: xn3,yn3,zn3,fn3
      real(rprec) :: dx2,dy2,dz2,bf2,df2
      real(rprec) :: rtes1,rtes2,rtes3,rtesf,rtes
      real(rprec) :: delx,dely,delz,delf
      real(rprec) :: h1,h2
      REAL(rprec) :: pos(3), bvec(3)
      real(rprec),parameter :: ftv_fact = 1. ! used to be 100
      integer (iprec), parameter :: imax = 50000
      integer(iprec) :: ierr,i

      ! get the normalized field components dx1,dy1,dz1, and the magnitude bf1
      pos(1)=xa
      pos(2)=ya
      pos(3)=za
      CALL Bfield_fric (pos,bvec,ierr)
      bf1 = SQRT (bvec(1)**2+bvec(2)**2+bvec(3)**2)
      if (ierr.lt.0.or.bf1.eq.0.0)return
      dx1 = bvec(1)/bf1
      dy1 = bvec(2)/bf1
      dz1 = bvec(3)/bf1

      DO i=1,imax

         h1 = h * ep ! go forward or backward along the field line?
         h2 = 0.5 * h1 ! h1 = whole step, h2 = half step

         ! First take a whole step
         xn1 = xa + h1*dx1
         yn1 = ya + h1*dy1
         zn1 = za + h1*dz1

         pos(1) = xn1
         pos(2) = yn1
         pos(3) = zn1
         CALL Bfield_fric (pos,bvec,ierr)
         bft = SQRT (bvec(1)**2+bvec(2)**2+bvec(3)**2)
         if (ierr.lt.0.or.bft.eq.0.0)return
         dxt = bvec(1)/bft
         dyt = bvec(2)/bft
         dzt = bvec(3)/bft
         fn1 = fa + 2.0*abs(h1/(bf1+bft))


         ! Now take a half step
          xn2 = xa + h2*dx1
          yn2 = ya + h2*dy1
          zn2 = za + h2*dz1

          pos(1) = xn2
          pos(2) = yn2
          pos(3) = zn2
          CALL Bfield_fric (pos,bvec,ierr)
          bft = SQRT(bvec(1)**2+bvec(2)**2+bvec(3)**2)
          if (ierr.lt.0.or.bft.eq.0.0)return
          dxt = bvec(1)/bft
          dyt = bvec(2)/bft
          dzt = bvec(3)/bft
          fn2 = fa + 2.0*abs(h2/(bf1+bft))


          ! Now take another half step (two halves make a whole?)
          pos(1) = xn2
          pos(2) = yn2
          pos(3) = zn2
          CALL Bfield_fric (pos,bvec,ierr)
          bf2 = SQRT (bvec(1)**2+bvec(2)**2+bvec(3)**2)
          IF (ierr.lt.0.or.bf2.eq.0.0) return
          dx2 = bvec(1)/bf2
          dy2 = bvec(2)/bf2
          dz2 = bvec(2)/bf2

          xn3 = xn2 + h2*dx2
          yn3 = yn2 + h2*dy2
          zn3 = zn2 + h2*dz2


          pos(1) = xn3
          pos(2) = yn3
          pos(3) = zn3
          CALL Bfield_fric (pos,bvec,ierr)
          bft = SQRT (bvec(1)**2+bvec(2)**2+bvec(3)**2)
          if (ierr.lt.0.or.bft.eq.0.0)return
          dxt = bvec(1)/bft
          dyt = bvec(2)/bft
          dzt = bvec(2)/bft
          fn3 = fn2 + 2.0*abs(h2/(bf2+bft))

          ! Compare the outcome of the two half steps with that of the whole step
          delx = xn1 - xn3
          dely = yn1 - yn3
          delz = zn1 - zn3                
          delf = fn1 - fn3                

          ! Use that to calculate an error parameter (rtes)
          rtes1 = delx*delx
          rtes2 = dely*dely
          rtes3 = delz*delz
          rtesf = delf*delf*ftv_fact ! decrease tolerance of ftv by ftv_fact

          rtes = sqrt(rtes1 + rtes2 + rtes3 + rtesf)
!         rtes = rtes/sqrt(xa**2+ya**2+za**2) 


         ! If the error is less than the error tolerance, then we are okay
         ! We can stop refining the step size
          if (rtes <= er) EXIT

!         bb = bf2

          h = 0.5*h   ! Otherwise, we will cut the step size in half and keep trying.

          if (i == imax) then ! problem
            write(*,*)' error in tracf2, i > imax - aborting program'
            write(*,*)' xa ya za =',xa,ya,za
            write(*,*)' er, h =', er,h
            STOP
          end if
 
      END DO


! So we have successfully made one trace step.  Since the calculations are
! already done, we might as well use the info we have.  We are using the
! euler method to trace field lines (which has a first order error term).
! We have a good idea what the error is, so we can subtract it from the
! solution to improve our estimate (to second order, I think).

      xn = xn3 - delx
      yn = yn3 - dely
      zn = zn3 - delz
      fn = fn3 - delf

! check to see if you can increase the stepsize
! if the error is really small, we may be taking unnecessarily small steps.
! h will be returned and used as our initial step size estimate next time.
      if(rtes<(0.25*er))h = 2.0*h

      RETURN
      END SUBROUTINE Tracf2
!----------------------------------------------------------





   subroutine Bfield_fric(r,b,ierr)

!-----------------bfield_fric----------------------
! returns bfieild at a specified point (x,z) in the form
! bx bz from the friction
! code grid
! works in code normalized units
!--------------------------------------------------

    USE coredata, only: bfactor,rprec,iprec, dm, tilt_angle
    implicit none
    real(rprec), intent(IN) :: r(2)
    real(rprec), intent(OUT) ::b(2)

    real(rprec) :: x,y,z,bbx,bbz,hx,hz,bf
    integer(iprec) :: ierr


    INTERFACE
       SUBROUTINE Dipole_modified (dm, tilt, x,z,bbx,bbz)
          IMPLICIT NONE
          REAL, INTENT (IN) :: dm, tilt, x,z
          REAL, INTENT (OUT) :: bbx,bbz
       END SUBROUTINE Dipole_modified
    END INTERFACE

    x = r(1); z = r(2);

! since the grid is good for z > 0 only assume symmetry

    CALL Bgrid (1,x,z,bbx,bbz,ierr)
    if(ierr < 0)return

    CALL Dipole_modified (dm,tilt_angle,x,z,hx,hz)

    bbx =  bbx + (hx)/(bfactor*1.0e9) ! convert to code units
    bbz =  bbz + (hz)/(bfactor*1.0e9)

    b(1) = bbx; b(2) = bbz

    bf = sqrt(bbx**2+bbz**2)

    IF (bf == 0.0)then
        write(6,*)' error in bfield, b=0 at'
        write(6,*)x,y,z
        ierr = -2
        return
    END IF

    RETURN
    END SUBROUTINE Bfield_fric

    SUBROUTINE Bgrid (isw, xp_in,zp_in, bxp, bzp, ierr)

!========================================================
! this subroutine reads in data from tapeXX.dat and
! returns bfield from friction code (delta B or pressure)
! depending on the input flag
! isw = 0 reads in the data
! isw = 1 returns bfield (delta B) in code normalized units
! isw = 2 returns pressure in code normalized units
!========================================================
    use coredata, ONLY: ssx,ssz,bxj,bzj,p,iprec,rprec,zahlx,zahlz

    implicit none

    integer(iprec), intent(IN) :: isw

    integer(iprec) :: idim,kdim
    
    integer(iprec), SAVE :: imax,kmax
    integer(iprec) :: i,k
    real(rprec) :: x(zahlx),z(zahlz)
    real(rprec) :: xtrans,ztrans
    real(rprec) :: xp,zp,bxp,bzp
    real(rprec) :: xp_in,yp_in,zp_in
    real(rprec) :: v1,v2,v3,v4,v5,v6,v7,v8,vt
    real(rprec) :: b1,b2,b3,b4,b5,b6,b7,b8
    real(rprec) :: b,r,rmin,hx,hy,hz
    real(rprec) :: bxmax,bymax,bzmax,pmax,pmin
    real(rprec) :: xbxmax,ybxmax,zbxmax
    real(rprec) :: xbymax,ybymax,zbymax
    real(rprec) :: xbzmax,ybzmax,zbzmax
    real(rprec) :: xpmax,ypmax,zpmax
    real(rprec) :: xpmin,ypmin,zpmin
    integer(iprec) :: ip,jp,kp
    integer(iprec) :: isym
    integer(iprec) :: ierr

    

    data rmin /1.0/
    data isym /1/

 ! convert to code unit
! isym = 1 for the friction code that assumes symmetry

    character(LEN=50) file_in
    idim=zahlx
    kdim=zahlz
    ierr = 0
    bxp = 0.0
    bzp = 0.0

        x = ssx
        z = ssz
        imax = zahlx
        kmax = zahlz


    ! calculate bfield
    ! check to see if the point is outside the grid, if so return null

        xp = xp_in
        zp = zp_in
        if(isym == 1)zp=abs(zp_in)

        if(xp < x(1) .OR. xp > x(imax) .OR.  zp < z(1) .OR. zp > z(kmax)) then
          write(*,*) 'xp = ', xp
          write(*,*) 'zp = ', zp
          write(*,*) 'x(imax) = ', x(imax)
          write(*,*) 'x(1) = ', x(1)
          write(*,*) 'Coordinates are outside the modeling region in Bgrid'
          return
        endif

    ! first locate where you are in the grid

        do i=1,imax-1
            if(x(i) <= xp .AND. x(i+1) >= xp)then
                ip = i
                go to 2
            end if
        end do
        write(6,*)' ip not found for x = ',xp
        2 continue


        do k=1,kmax-1
            if(z(k) <= zp .AND. z(k+1) >= zp)then
                kp = k
                go to 4
            end if
        end do
        write(6,*)' kp not found for z = ',zp
        4 continue

    ! now do the interpolation

    ! (i,j+1,k+1) +-------------+(i+1,j+1,k+1)
    ! / 7  /  8     /|
    ! /----/--------/8|
    ! /  3 /   4    /| |
    ! (i,j+1,k) +------------+4|/|
    ! |     |       | / |
    ! |  3  |   4   |/|5+ (i+1,j,k+1)
    ! |-----+-------|2|/
    ! |     |       | /
    ! |  1  |   2   |/
    ! +-------------+
    ! (i,j,k)     (i+1,j,k)

    ! volume label    vertex
    ! 1         (i  ,j  ,k  )
    ! 2         (i+1,j  ,k  )
    ! 3         (i  ,j+1,k  )
    ! 4         (i+1,j+1,k  )
    ! 5         (i  ,j  ,k+1)
    ! 6         (i+1,j  ,k+1)
    ! 7         (i  ,j+1,k+1)
    ! 8         (i+1,j+1,k+1)
    
        v1 = (xp-x(ip))*  (zp-z(kp))
        v2 = (x(ip+1)-xp)*    (zp-z(kp))
        v3 = (xp-x(ip))*    (zp-z(kp))
        v4 = (x(ip+1)-xp)*  (zp-z(kp))
        v5 = (xp-x(ip))*      (z(kp+1)-zp)
        v6 = (x(ip+1)-xp)*    (z(kp+1)-zp)
        v7 = (xp-x(ip))*    (z(kp+1)-zp)
        v8 = (x(ip+1)-xp)* (z(kp+1)-zp)

        vt = (x(ip+1)-x(ip))*(z(kp+1)-z(kp))


    ! write(6,*)' vol diff ',vt-v1-v2-v3-v4-v5-v6-v7-v8
    ! now for the values
    ! if isw = 2 return the pressure only
        if(isw == 2)then

            b1 = p(ip  ,kp)
            b2 = p(ip+1,kp)
            b3 = p(ip  ,kp)
            b4 = p(ip+1,kp)
            b5 = p(ip  ,kp+1)
            b6 = p(ip+1,kp+1)
            b7 = p(ip  ,kp+1)
            b8 = p(ip+1,kp+1)
        ! remove weight f any points are zero

            if(b1 == 0.0)v8=0.0
            if(b2 == 0.0)v7=0.0
            if(b3 == 0.0)v6=0.0
            if(b4 == 0.0)v5=0.0
            if(b5 == 0.0)v4=0.0
            if(b6 == 0.0)v3=0.0
            if(b7 == 0.0)v2=0.0
            if(b8 == 0.0)v1=0.0

            vt = v1+v2+v3+v4+v5+v6+v7+v8

            if(vt <= 0.0) then
              write(*,*) 'vt <= 0.0 in Bgrid'
              return
            endif

            bxp = (v1*b8+v2*b7+v3*b6+v4*b5+v5*b4+v6*b3+v7*b2+v8*b1)/vt
            bzp = bxp

        ! if(b1.eq.0.0.or.b2.eq.0.0.or.b3.eq.0.0.or.b4.eq.0.0.or.
        ! $   b5.eq.0.0.or.b6.eq.0.0.or.b7.eq.0.0.or.b8.eq.0.0)then
        ! xp = max1(b1,b2,b3,b4,b5,b6,b7,b8)
        ! rite(6,*)' debug point in bgrid, pressure, xp,yp=', xp,yp
        ! end if


        else

            b1 = bxj(ip  ,kp)
            b2 = bxj(ip+1,kp)
            b3 = bxj(ip  ,kp)
            b4 = bxj(ip+1,kp)
            b5 = bxj(ip  ,kp+1)
            b6 = bxj(ip+1,kp+1)
            b7 = bxj(ip  ,kp+1)
            b8 = bxj(ip+1,kp+1)

            bxp = (v1*b8+v2*b7+v3*b6+v4*b5+v5*b4+v6*b3+v7*b2+v8*b1)/vt


            
            b1 = bzj(ip  ,kp)
            b2 = bzj(ip+1,kp)
            b3 = bzj(ip  ,kp)
            b4 = bzj(ip+1,kp)
            b5 = bzj(ip  ,kp+1)
            b6 = bzj(ip+1,kp+1)
            b7 = bzj(ip  ,kp+1)
            b8 = bzj(ip+1,kp+1)

            bzp = (v1*b8+v2*b7+v3*b6+v4*b5+v5*b4+v6*b3+v7*b2+v8*b1)/vt

            if(isym == 1 .AND. zp_in < 0.0)then
                bxp = - bxp
            end if

        end if


    return

    end SUBROUTINE Bgrid


subroutine setfri(fric,zeit,diss,ideal)

! to set frictional coefficient, normalized to one


    USE coredata, only: iprec,rprec
    implicit none
    real(rprec) :: fric,zeit,diss,ideal,h1,h2,h3,h4
    real(rprec),parameter ::pi=3.14159265

! need to calculate the appropriate time frame
! but first perform some relaxation initially

    if (zeit <= 20.) then
        fric = cos(zeit/40.*pi)+0.01
    else
        fric = 0.01
    endif
! h1 = diss+ideal
! h2 = h1*(zeit/h1 - int(zeit/h1))
! if (h2.ge.ideal) then
! h3 = h2-ideal
! fric = sin(pi*h3/diss)
! else
! fric = 0.
! endif
! endif
! fric = 0.
    return
    end subroutine setfri

!===================================================================


    subroutine diag(bz,vx,rho,p,ssx,zeit,in)

    USE coredata, only: iprec,rprec,zahlx,  zahlz
    implicit none

    real(rprec) :: bz(zahlx,zahlz),vx(zahlx,zahlz), &
    rho(zahlx,zahlz), ssx(zahlx), zeit, x &
    ,p(zahlx,zahlz), temp
    integer(iprec) :: in, ix,  iz

! open diag file

! open(20,file='diag.txt',form='FORMATTED',status='unknown')


    write(20,10)in, zeit
    10 format(1x,'Diagnosis: x-axis  interation:', &
    i6,' time:',f20.8)
    write(20,20)
    20 format(1x,'  x      ',1x,'       bzj    ',1x,'         vx    ', &
    1x,'       rho    ',1x,'         p     ')

    iz = 2

    do ix =1,zahlx
        x = ssx(ix)
        write(20,30)x,bz(ix,iz),vx(ix,iz),rho(ix,iz), &
        p(ix,iz)
        30 format(1x,f10.5,4(1x,e15.5))
    end do

! close diag file

! close(20)
    end subroutine diag




    subroutine limit3(u,rx,rz,vis)

! modified by frt 3/95 to include variable y-spacing

! the new smoothing routine that does not affect the boundary, saves
! time and memory and preserves Integral u dV

    USE coredata, only: iprec,rprec,zahlx,zahly,zahlz
    implicit none

    real(rprec) :: u(zahlx,zahlz), rx(zahlx), ry(zahly),rz(zahlz)
    real(rprec) :: vis, d2f1, d2f2, h1, h2, q
! parameter (vis = 0.25)
    integer(iprec) :: ix,  iz


!================= x-Richtung ===================================================

! Try a forall version, for significant speed imporvement
! zero things to make sure the boundaries aren't getting changed

! this won't be the same, because in the loop version, u is updated, &
! and that updated u is used in future calculations 



    do iz = 2, zahlz-1
            do ix = 4, zahlx-2
                d2f1=u(ix-2,iz)-2.0*u(ix-1,iz)+u(ix,iz)
                d2f2=u(ix-1,iz)-2.0*u(ix,iz)+u(ix+1,iz)
                q = (0.5-sign(0.5,d2f1*d2f2))*sign(1.0,d2f2)*vis*min(abs(d2f1),abs(d2f2))
                u(ix,iz)=u(ix,iz) + q
                u(ix-1,iz)=u(ix-1,iz) - q*(rx(ix-1)/rx(ix))
        enddo
    enddo



!================= z-Richtung ===================================================

! here we need to pay attention to the variable grid

    do iz = 4, zahlz-2
            do ix = 2, zahlx-1
                d2f1=u(ix,iz-2)-2.0*u(ix,iz-1)+u(ix,iz)
                d2f2=u(ix,iz-1)-2.0*u(ix,iz)+u(ix,iz+1)
                q = (0.5-sign(0.5,d2f1*d2f2))*sign(1.0,d2f2)*vis*min(abs(d2f1),abs(d2f2))
                u(ix,iz)=u(ix,iz) + q
                u(ix,iz-1)=u(ix,iz-1) - q*(rz(iz-1)/rz(iz))
            enddo
    enddo


    return
    end subroutine limit3

!------------------------------------------------------



!===============================================================================

    subroutine bdp (result,bx,bz,p,rx,rz,ssx,ssz)

    USE coredata, only : rprec,iprec,zahlx,zahly,zahlz
    implicit none

    real(rprec) :: bx(zahlx,zahlz)&
    ,bz(zahlx,zahlz),p(zahlx,zahlz) &
    ,ssx(zahlx),ssz(zahlz)
    real(rprec) :: rx(zahlx), rz(zahlz), result, new, &
    help(zahlx,zahlz), hh
    integer(iprec) :: ix, iz

    do ix=3,zahlx-2
            do iz=3,zahlz-2
                hh = (bx(ix,iz)*(p(ix+1,iz)-p(ix-1,iz))*rx(ix)) &
                +(bz(ix,iz)*(p(ix,iz+1)-p(ix,iz-1))*rz(iz))
                hh = abs(hh) &
                /sqrt(bx(ix,iz)**2+bz(ix,iz)**2)
                help(ix,iz) = hh
            enddo
    enddo


    new = 0.
    call tdint3(help,new,ssx,ssz)
    result = new
    return
    end subroutine bdp

!===============================================================================


! calculates the average value of a function fun
! and returns funave of the abs
    subroutine averageabs(idim,kdim,fun,funave)

    USE coredata, only: iprec,rprec
    implicit none
    integer(iprec) :: idim,kdim,i,j,k
    real(rprec) :: fun(idim,kdim)
    real(rprec) :: funave,sum

    sum = 0.0

    do i=1,idim
            do k=3,kdim

                sum = sum + abs(fun(i,k))

            end do
    end do

    funave = sum/float(idim*(kdim-2))

    return

    end subroutine averageabs
!------------------------------------------------c
! and returns funave
    subroutine average(idim,kdim,fun,funave)
    

    implicit none
    integer :: idim,kdim,i,k
    real :: fun(idim,kdim)
    real :: funave,sum

    sum = 0.0

    do i=1,idim
            do k=3,kdim

                sum = sum + fun(i,k)

            end do
    end do

    funave = sum/float(idim*(kdim-2))

    return

    end subroutine average
!------------------------------------------------c
! calculates the minimun value of a function fun
! and returns   minave i,jk,
    subroutine find_min(idim,kdim,fun,funmin,im,km)

    USE coredata, only: iprec,rprec
    implicit none
    integer(iprec) :: idim,kdim,i,j,k,im,km
    real(rprec) :: fun(idim,kdim)
    real(rprec) :: funmin

    funmin = fun(1,1)

    do k=1,kdim
            do i=1,idim

                if(fun(i,k) <= funmin)then
                    funmin = fun(i,k)
                    im = i
                    km = k
                end if

            end do
    end do


    return

    end subroutine find_min
!------------------------------------------------c
! sets a floor on a variable

    subroutine floor_it(idim,kdim,fun,funmin)

    USE coredata, only: iprec,rprec
    implicit none
    integer(iprec) :: idim,kdim,i,k
    real(rprec) :: fun(idim,kdim)
    real(rprec) :: funmin


    fun(:,:) = max(fun(:,:), funmin)

    return

    end subroutine floor_it
!------------------------------------------------c
! sets a floor on a variable rho

    subroutine floor_rho(idim,kdim,fun)

    USE coredata, only: iprec,rprec
    implicit none
    integer(iprec) :: idim,kdim,i,j,k
    real(rprec) :: fun(idim,kdim)
    real(rprec) :: rhomax,rhomin,funmin0

    data rhomax /0.3/
    data rhomin /0.01/


    do k=1,kdim
            do i=1,idim

                funmin0=rhomin + (rhomax-rhomin)* &
                (float(i-1)/float(idim-1))**4

                fun(i,k) = max(fun(i,k),funmin0)

            end do
    end do


    return

    end subroutine floor_rho
!-------------------------------------
    subroutine zero_it(nx,nz,array)

    USE coredata, only: iprec,rprec
    implicit none

    integer(iprec) :: nx,nz
    integer(iprec) :: i,k
    real(rprec) :: array(nx,nz)

    do k=1,nz
            do i=1,nx

                array(i,k) = 0.0

            end do
    end do

    return
    end subroutine zero_it



    SUBROUTINE Set_bc0 (nx,nz,ssx,ssz,rx,rz,ddx,ddz, &
                     u,rho,vx,vz,sx,sz,bx,bz,in,zeit,dt,iopen)

   ! iopen = 0 - closed bc at back tail
   !       = 1 - open bc at back tail

    USE coredata, only: iprec,rprec,symmetric
    implicit none

    integer(iprec), intent (in) :: nx,nz, iopen, in

    integer(iprec) :: nx1,nz1
    integer(iprec) :: nx2,nz2
    integer(iprec) :: ix,iz
    integer(iprec) :: iv(2,3,8)

    real(rprec) :: ssx(nx),ssz(nz)
    real(rprec) :: ddx(nx,nz),ddz(nx,nz)
    real(rprec) :: rx(nx),rz(nz)
    real(rprec) :: u(nx,nz),rho(nx,nz) &
    ,vx(nx,nz),vz(nx,nz) &
    ,sx(nx,nz),sz(nx,nz) &
    ,bx(nx,nz),bz(nx,nz)

    real(rprec) :: zeit,dt
    real(rprec) :: xmax,ymax,zmax
    real(rprec) :: fac1, fac2

! BOUNDARY CONDITIONS FOR V
! __________________________________________________________________

! 1. index  1=min, 2=max
! 2. index  1# x=const, 2# y=const, 3# z=const
! 3. index  komponente: 1-3: sx, sy, sz, 4: rho
! iv =  1 => symmetric
! iv = -1 => antisymmetric, i.e, hom. Dirichlet B. C.

! x=0
    iv(1,1,1)=-1
    iv(1,1,2)=-1
    iv(1,1,3)=-1
    iv(1,1,4)=1
! x=xmax
    iv(2,1,1)=-1
    iv(2,1,2)=-1
    iv(2,1,3)=-1
    iv(2,1,4)=1
! y=0
    iv(1,2,1)=-1
    iv(1,2,2)=-1
    iv(1,2,3)=-1
    iv(1,2,4)=1
! y=ymaxxc
    iv(2,2,1)=-1
    iv(2,2,2)=-1
    iv(2,2,3)=-1
    iv(2,2,4)=1
! z=0
    if(symmetric)then
    iv(1,3,1)=1
    iv(1,3,2)=1
    iv(1,3,3)=-1
    iv(1,3,4)=1
    else
    iv(1,3,1)=-1
    iv(1,3,2)=-1
    iv(1,3,3)=-1
    iv(1,3,4)=1
    end if
! z=zmax
    iv(2,3,1)=-1
    iv(2,3,2)=-1
    iv(2,3,3)=-1
    iv(2,3,4)=1

! end boundary cond.'s for v
!===================================================================

    nx1 = nx - 1;  nz1 = nz - 1
    nx2 = nx - 2;  nz2 = nz - 2




    ! scalar fields (density, pressure)

    DO iz=2,nz1 
       rho (1,iz) = rho (3,iz)
       u   (1,iz) = u   (3,iz)
       rho (nx,iz)= rho (nx2,iz)
       u   (nx,iz)= u   (nx2,iz)
    END DO


    DO ix=1,nx 
       rho (ix,1)  = rho (ix,3)
       u   (ix,1)  = u   (ix,3)
       rho (ix,nz) = rho (ix,nz2)
       u   (ix,nz) = u   (ix,nz2)
    END DO

    ! vector fields (velocity, momentum):

    fac1 = -1.0; if (iopen == 1) fac1=1.0  ! back tail bc switch for ghostcells
    fac2 =  0.0; if (iopen == 1) fac2=1.0  ! back tail bc switch for 1st rows of real cells

    DO iz=2,nz1;

       sx  ( 1,iz) = fac1*sx  (  3,iz)
       sz  ( 1,iz) = fac1*sz  (  3,iz)
       vx  ( 1,iz) = fac1*vx  (  3,iz)
       vz  ( 1,iz) = fac1*vz  (  3,iz)
       sx  (nx,iz) =  -   sx  (nx2,iz)
       sz  (nx,iz) =  -   sz  (nx2,iz)
       vx  (nx,iz) =  -   vx  (nx2,iz)
       vz  (nx,iz) =  -   vz  (nx2,iz)

       sx(nx1,iz)=0.
       sz(nx1,iz)=0.
       vx(nx1,iz)=0.
       vz(nx1,iz)=0.
       sx(2,iz)=fac2*sx(3,iz)
       sz(2,iz)=fac2*sz(3,iz)
       vx(2,iz)=fac2*vx(3,iz)
       vz(2,iz)=fac2*vz(3,iz)
     END DO




     DO ix=1,nx 
       sx  (ix,1)  =  iv(1,3,1)*sx(ix,3)
       sz  (ix,1)  =  iv(1,3,3)*sz(ix,3)
       vx  (ix,1)  =  iv(1,3,1)*vx(ix,3)
       vz  (ix,1)  =  iv(1,3,3)*vz(ix,3)

       sx  (ix,nz) = -sx(ix,nz2)
       sz  (ix,nz) = -sz(ix,nz2)
       vx  (ix,nz) = -vx(ix,nz2)
       vz  (ix,nz) = -vz(ix,nz2)
     END DO


    DO ix=1,nx
       if(iv(1,3,1) < 0) sx(ix,2)=0.
       if(iv(1,3,3) < 0) sz(ix,2)=0.

       sx(ix,nz1)=0.
       sz(ix,nz1)=0.

       if(iv(1,3,1) < 0) vx(ix,2)=0.
       if(iv(1,3,3) < 0) vz(ix,2)=0.

       vx(ix,nz1)=0.
       vz(ix,nz1)=0.
     END DO



    ! vector field B (div(B)=0 is special case):

    DO iz=2,nz1
       bx  ( 1,iz) = bx  (3,iz)
       bz  ( 1,iz) = bz  (3,iz)
       bx  (nx,iz) = bx  (nx2,iz)
       bz  (nx,iz) = bz  (nx2,iz)
    END DO


     DO ix=1,nx 
     if(symmetric)then
       bx  (ix,1)  = -bx(ix,3)
       bz  (ix,1)  =  bz(ix,3)
     else
       bx  (ix,1)  =  bx(ix,3)
       bz  (ix,1)  =  bz(ix,3)
     end if

       bx  (ix,nz) =  bx(ix,nz2)
       bz  (ix,nz) =  bz(ix,nz2)
    END DO


     Do iz = 2,nz1

          bx (1,iz) = bx(3,iz) + 1.0/rx(2)*(bz(2,iz+1)-bz(2,iz-1))*rz(iz)
         
          bx (nx,iz) = bx(nx2,iz) -1.0/rx(nx-1)*(bz(nx1,iz+1)-bz(nx1,iz-1))*rz(iz)
     enddo


    RETURN
    END SUBROUTINE Set_bc0




!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! subroutine change_dt:  Changes dt by a factor red
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    subroutine change_dt(dt,red,dt_max,dt_min,in_tot,zeit)

    USE coredata, only: iprec,rprec,dirname
    implicit none

    real(rprec) :: dt,red,dt_max,dt_min,zeit
    integer(iprec) :: in_tot

    dt = dt*red

    dt = max(dt,dt_min)
! dt = min(dt,dt_max)

    open(97,file=trim(dirname)//'/'//'time.dat',status='unknown',position='append')
    write(97,*)in_tot,zeit,dt
    close(97)

    return
    end subroutine change_dt



    SUBROUTINE Courant (nx,nz,ddx,ddz, vx,vz,bx,bz,p,rho,gamma,tstep,cfl)

    USE coredata, ONLY: iprec,rprec
    implicit none

    integer(iprec) :: nx,nz
    real(rprec),dimension(nx,nz),intent(in) :: &
                 ddx,  ddz, vx,  vz, bx,  bz, p, rho
    real(rprec),intent(in) :: gamma, cfl
    real(rprec),intent(out) :: tstep


    integer(iprec) :: ix,iz
    real(rprec) :: cour(nx,nz), vfast (nx,nz)


    vfast = SQRT((bx*bx+bz*bz + gamma*p)/rho)

    DO ix=2,nx-1; DO iz=2,nz-1
        vfast(ix,iz) = SQRT((bx(ix,iz)**2 +  bz(ix,iz)**2 + &
                                gamma*p(ix,iz))/rho(ix,iz))
        cour(ix,iz) = 1.0/MAX((ABS(vx(ix,iz))+vfast(ix,iz))*ddx(ix,iz), &
                                 (ABS(vz(ix,iz))+vfast(ix,iz))*ddz(ix,iz))
     END DO; END DO

    tstep = MINVAL(cour(2:nx-1,2:nz-1))*cfl

    RETURN
    END SUBROUTINE Courant



!------------------------------------------------------
! div b cleaner based on the algorixthm of Marder
! from the paper by
! Landgon, Computer Physixcs Comm, 70, 447-450,1992
! 11/97 frt
! uses the simplest centered scheme


    subroutine smoothf(idim,kdim,x,z,f,vis)
    USE coredata, only: iprec,rprec
    implicit none
    integer(iprec) :: idim,kdim
    real(rprec) :: x(idim),z(kdim)
    real(rprec) :: f(idim,kdim)
    real(rprec) :: ftemp(idim,kdim)
    real :: vis

    integer(iprec) :: i,k

    ftemp = f

    do i=2,idim-1
            do k=2,kdim-1

                ftemp(i,k) =((f(i-1,k) + f(i+1,k) + &
                                f(i,k-1) + f(i,k+1))/6.0 + &
                                f(i,k) )/2.0

          enddo
     enddo

        do k=1,kdim
            ftemp(1,k) = ftemp(2,k)
            ftemp(idim,k)=ftemp(idim-1,k)
         enddo


! use 3 for k =1 since k=2 is z=0
    do i=1,idim
            ftemp(i,1) = ftemp(i,3)
            ftemp(i,kdim)=ftemp(i,kdim-1)
     enddo
! now reset f
     f = (1.-vis)*f+vis*ftemp


    return
    end subroutine smoothf


!------------------------------------------------------

!===============================================================================

    subroutine tdint3(f,res,ssx,ssz)

! intgrates over the volume - frt uses trapezoidal rule
! avoid the boundary

    USE coredata, only: iprec,rprec,zahlx,zahlz
    implicit none

    real(rprec) :: f(zahlx,zahlz), twodf(zahlx), &
    onedf(zahlx), rz(zahlz), res, rx(zahlx)
    integer(iprec) :: ix,  iz, nx1,ny1,nz1
    real(rprec) :: ssx(zahlx), ssz(zahlz)

    nx1=zahlx-1
    nz1=zahlz-1

! now integrate in the z-direction from z=0

    do ix=1,zahlx
            twodf(ix) = 0.0
            do iz=2,zahlz-1
                twodf(ix)= twodf(ix)+ &
                (f(ix,iz)+f(ix,iz+1))*(ssz(iz+1)-ssz(iz))*0.5
            enddo
    enddo

! now in the y-direction

! now in the x-direction

    res = 0.0
    do ix=1,zahlx-1
        res = res+(onedf(ix+1)+onedf(ix))*(ssx(ix+1)-ssx(ix))*0.5
    end do

    return

    end subroutine tdint3


      real function ftv_dip(x,y)
      ! function that returns the analytic flux tube volume

      USE coredata, only: iprec,rprec, dm
      IMPLICIT NONE

      real(rprec), intent(in) :: x,y
      real(rprec) ::  l
!     real(rprec), parameter :: be0 = 30574. ! dipole moment
! uses an analytic approximation

      l = sqrt(x**2+y**2)

      ftv_dip = 0.0

      l = sqrt(x**2+y**2)

      if (l > 1.0)then

!       ftv_dip = 32./35.*l**4/be0*sqrt(1.0-1./l)*&
        ftv_dip = 32./35.*l**4/dm*sqrt(1.0-1./l)*&
                  (1.+1./2./l+3./8./l**2+5./16./l**3)

      end if

end function ftv_dip





SUBROUTINE Compute_grid_noise (nx,nz, gsn, gsnave, gsnave_tail, gsnave_inner)

   USE Coredata, ONLY : iprec, rprec, Volume_avg, &
                        p, global_mask, inner_mask, tail_mask, &
                        nx1, nx2, nx3, nx4,  ny3, ny4, nz1, nz2, nz3, nz4
   IMPLICIT NONE
   INTEGER (iprec), INTENT (IN) :: nx, nz
   REAL (rprec), INTENT (OUT) :: gsn (nx,nz), gsnave, gsnave_tail, gsnave_inner

   gsn          = 0.0_rprec
   gsnave       = 0.0_rprec
   gsnave_tail  = 0.0_rprec
   gsnave_inner = 0.0_rprec

   gsn(3:nx2,3:nz2) = ((p(1:nx4,3:nz2) -4*p(2:nx3,3:nz2)  &
           + 6*p(3:nx2,3:nz2) - 4*p(4:nx1,3:nz2) + p(5:nx,3:nz2)) &
           +  (p(3:nx2,3:nz2) -4*p(3:nx2,3:nz2) + 6*p(3:nx2,3:nz2) &
           - 4*p(3:nx2,3:nz2) + p(3:nx2,3:nz2)) &
           +  (p(3:nx2,1:nz4) -4*p(3:nx2,2:nz3) + 6*p(3:nx2,3:nz2) &
           - 4*p(3:nx2,4:nz1) + p(3:nx2,5:nz)))/p(3:nx2,3:nz2)

   gsnave       = Volume_avg (ABS(gsn),global_mask)
   gsnave_inner = Volume_avg (ABS(gsn),inner_mask)
   gsnave_tail  = Volume_avg (ABS(gsn),tail_mask)

   RETURN
END SUBROUTINE Compute_grid_noise






SUBROUTINE Compute_energy (nx, nz, wb, wp, wv, w, free)
   USE Coredata, ONLY : iprec, rprec, r, rmin, &
                        bx,bz,p,rho,vx,vz,vel,gamma, &
                        nx2, ny2, nz2, nx1, ny1, nz1
   IMPLICIT NONE
   INTEGER (iprec), INTENT(IN) :: nx,nz
   REAL (rprec), INTENT (OUT) :: wb(nx,nz), wp(nx,nz), wv(nx,nz),&
                                 w(nx,nz), free(nx,nz)

   INTEGER (iprec) :: ix,  iz


   ! find magnitude of velocity:
   forall (ix=3:nx2, iz=2:nz2)
     vel(ix,iz) = SQRT(vx(ix,iz)**2 +  vz(ix,iz)**2)
   end forall


   where (r > 1.0*rmin)
       wb = 0.5*(bx**2 +  bz**2)
       wp = p/(gamma - 1.0)
       wv = 0.5*rho*vel*vel
   elsewhere
       wb = 0.0
       wp = 0.0
       wv = 0.0
   endwhere

   ! set equal to zero at the boundaries
   ! this will lower the integrated values, but avoid boundary weirdness            
   wb(1:2,  : ) = 0.0
   wb( : ,  1 ) = 0.0
   wb(nx1:nx, :  ) = 0.0
   wb(   :  , nz1:nz) = 0.0

   wp(1:2,  : ) = 0.0
   wp( : ,  1 ) = 0.0
   wp(nx1:nx,  :  ) = 0.0
   wp(   : ,nz1:nz) = 0.0

   wv(1:2,  : ) = 0.0
   wv( : ,  1 ) = 0.0
   wv(nx1:nx,      :  ) = 0.0
   wv(   :  ,   nz1:nz) = 0.0

   free = wb - p
   w = wb + wp + wv

   RETURN
END SUBROUTINE Compute_energy



SUBROUTINE Compute_velocity (nx,nz, sx,sz, rho, vx,vz)

   USE Coredata, ONLY : iprec, rprec

   IMPLICIT NONE
   INTEGER(iprec), INTENT (IN) :: nx,nz
   REAL (rprec), INTENT (IN) :: sx(nx,nz),  sz(nx,nz), rho(nx,nz)
   REAL (rprec), INTENT (OUT) :: vx(nx,nz), vz(nx,nz)


   vx = sx/rho
   vz = sz/rho

   RETURN
END SUBROUTINE Compute_velocity

    SUBROUTINE Set_frictional_term (fric_alg, fric_min, fric_max, ekave_old, ekave, zeit, frica)

    !---------------- set the frictional term----------------

       USE Coredata, ONLY : rprec,iprec
       IMPLICIT NONE
       integer (iprec), intent (in) :: fric_alg
       real (rprec), intent (in) :: fric_min, fric_max, ekave_old, ekave, zeit
       real (rprec), intent (in out) :: frica

       real (rprec) :: fric_mint

       if (fric_alg == 1)then

          if(ekave <= ekave_old)then
              frica = frica*0.5
          else
              frica = frica*2.0
          end if

          ! clamp the values of frica
          frica = max(frica,fric_min)
          frica = min(frica,fric_max)

       else if(fric_alg == 2)then

          ! set values of fric_min
          if(zeit < 40.)then
              fric_mint = fric_max*cos(3.1415*zeit/80.) + fric_min
          else
              fric_mint = fric_min
          end if

          if(ekave < ekave_old)then
               frica = frica*0.5
          else
               frica = frica*2.0
          end if

          ! clamp the values of frica
          frica = max(frica,fric_mint)
          frica = min(frica,fric_max)

       else

          STOP 'invalid fric_alg, aborting'

       end if

       RETURN
    END SUBROUTINE Set_frictional_term



    SUBROUTINE Set_rho_lemon (zahlx,zahlz,ddx,ddz, &
                              sx,sz,bx,bz,p,rho,gamma,tstep,cfl)

       USE Coredata, ONLY : iprec, rprec, huge_rho
    
       IMPLICIT NONE

       integer(iprec),intent(in) :: zahlx,zahlz
       real(rprec),intent(in) :: ddx(zahlx,zahlz)
       real(rprec),intent(in) :: ddz(zahlx,zahlz)
       real(rprec),intent(in) :: bx(zahlx,zahlz)
       real(rprec),intent(in) :: bz(zahlx,zahlz)
       real(rprec),intent(in) :: sx(zahlx,zahlz)
       real(rprec),intent(in) :: sz(zahlx,zahlz)
       real(rprec),intent(in) :: p(zahlx,zahlz)
       real(rprec),intent(in out) :: rho(zahlx,zahlz)
       real(rprec),intent(in) :: gamma, cfl, tstep

       real(rprec) :: vnum(zahlx,zahlz)
       real(rprec) :: qb(zahlx,zahlz)
       real(rprec) :: b2(zahlx,zahlz)
       real(rprec) :: s2(zahlx,zahlz)

       vnum = 0.5*cfl/(tstep*max(ddx,ddz))
       b2 = bx**2+bz**2 + gamma*p
       s2 = sx**2+sz**2
       qb = -(2*sqrt(s2)*vnum + b2)

       where (rho < huge_rho)
          rho = (-qb + sqrt(qb**2 - 4*vnum**2*s2)) / (2*vnum**2)
       endwhere

       RETURN
    END SUBROUTINE Set_rho_lemon

!----------------------------------
       subroutine derivs_fric(x,ndims, f1)
! returns the rhs for the tracer including change in
! flux tube volume, so that f1(4) = ds/B
          USE coredata, only : iprec,rprec
          implicit none
          integer(iprec),intent(in) :: ndims
          real(rprec),intent(in) :: x(ndims)
          real(rprec),intent(out) :: f1(ndims)
          real(rprec) :: dxt,dzt,bft
          real(rprec) :: r(2),b(2)
          integer(iprec) :: ierr
          external bfield_fric

            f1(:) = 0.0
            r(1:2) = x(1:2)
!           write(6,*)' r =',r
          call bfield_fric(r,b,ierr)
          bft = sqrt(dot_product(b,b))
!         bft = sqrt(b(1)**2+b(2)**2+b(3)**2)
          if(ierr .ge. 0 .and. bft > 0.)then
            f1(1) = b(1)/bft
            f1(2) = b(2)/bft
            f1(4) = 1./bft
          end if
!         write(*,*)' rhs derivs x,y,z..=',f1
          return
       end subroutine derivs_fric




subroutine allocate_core_data(gridsizex,gridsizez)
use coredata
implicit none
integer :: gridsizex,gridsizez
zahlx=gridsizex
zahlz=gridsizez
call system('mkdir -p ' // adjustl(trim(dirname)))
    allocate(ssx(zahlx))                  ! grid x-values
    allocate(ssz(zahlz))                  ! grid z-values
    allocate(rx(zahlx))                   ! central diff d/dx values
    allocate(rz(zahlz))                   ! central diff d/dz values
    allocate(r(zahlx,zahlz))
    allocate(pr(zahlx,zahlz))
    allocate(pr2(zahlx,zahlz))
    allocate(vol(zahlx,zahlz))       ! volume around each grid point
    allocate(xxx(zahlx,zahlz))       ! x values at grid points
    allocate(zzz(zahlx,zahlz))       ! z values at grid points
    allocate(ddx(zahlx,zahlz))       ! rx values at grid points
    allocate(ddz(zahlx,zahlz))       ! rz values at grid points

    allocate(global_mask(zahlx,zahlz))
    allocate(tail_mask(zahlx,zahlz))
    allocate(inner_mask(zahlx,zahlz))

! density variable
    allocate(rho(zahlx,zahlz))   ! density
    allocate(rhoo(zahlx,zahlz))  ! density
    allocate(divrv(zahlx,zahlz))
    allocate(divrv_1(zahlx,zahlz))
    allocate(divrv_2(zahlx,zahlz))

! pressure variable
    allocate(p(zahlx,zahlz))    ! pressure
    allocate(po(zahlx,zahlz))    ! pressure
    allocate(dpv(zahlx,zahlz))
    allocate(divv(zahlx,zahlz))
    allocate(dp(zahlx,zahlz))   ! change in u (timestep)
    allocate(dp_1(zahlx,zahlz)) ! change in u (timestep)
    allocate(dp_2(zahlx,zahlz)) ! change in u (timestep)

! velocity and momentum variables
    allocate(sx(zahlx,zahlz))    ! x-momentum
    allocate(sz(zahlx,zahlz))   ! z-momentum
    allocate(sxo(zahlx,zahlz))  ! x-momentum
    allocate(szo(zahlx,zahlz))   ! z-momentum
    allocate(vx(zahlx,zahlz))    ! x-velocity
    allocate(vz(zahlx,zahlz))    ! z-velocity
    allocate(vxo(zahlx,zahlz))   ! x-velocity
    allocate(vzo(zahlx,zahlz))   ! z-velocity
    allocate(dsx(zahlx,zahlz))   ! x-force
    allocate(dsz(zahlx,zahlz))   ! z-force
    allocate(dsx_1(zahlx,zahlz)) ! x-force at previous timestep
    allocate(dsz_1(zahlx,zahlz)) ! z-force at previous timestep
    allocate(dsx_2(zahlx,zahlz)) ! x-force two timesteps ago
    allocate(dsz_2(zahlx,zahlz)) ! z-force two timesteps ago
    allocate(vel(zahlx,zahlz))   ! velocity magnitude

! magnetic field variables
    allocate(bx(zahlx,zahlz))    ! x-magnetic field
    allocate(bz(zahlx,zahlz))    ! z-magnetic field
    allocate(bxj(zahlx,zahlz))   ! Bx minus dipole
    allocate(bzj(zahlx,zahlz))   ! Bz minus dipole
    allocate(bxd(zahlx,zahlz))  ! dipolar Bx
    allocate(bzd(zahlx,zahlz))   ! dipolar Bz
    allocate(dbx(zahlx,zahlz))   ! change in Bx (timestep)
    allocate(dbz(zahlx,zahlz))   ! change in Bz (timestep)
    allocate(dbx_1(zahlx,zahlz)) ! change in Bx (timestep)
    allocate(dbz_1(zahlx,zahlz)) ! change in Bz (timestep)
    allocate(dbx_2(zahlx,zahlz)) ! change in Bx (timestep)
    allocate(dbz_2(zahlx,zahlz)) ! change in Bz (timestep)
    allocate(bxo(zahlx,zahlz))   ! old Bx (previous iteration)
    allocate(bzo(zahlx,zahlz))   ! old Bz
    allocate(bxjo(zahlx,zahlz))  ! old Bx (previous iteration)
    allocate(bzjo(zahlx,zahlz))  ! old Bz
    allocate(divb(zahlx,zahlz))  ! divergence of B

    allocate(ajx(zahlx,zahlz))   ! Jx current
    allocate(ajy(zahlx,zahlz))   ! Jy current
    allocate(ajz(zahlx,zahlz))   ! Jz current

! Force variables
    allocate(kx(zahlx,zahlz))    ! x-comp. of force JxB-grad(P)
    allocate(ky(zahlx,zahlz))    ! y-comp. of force JxB-grad(P)
    allocate(kz(zahlx,zahlz))   ! z-comp. of force JxB-grad(P)
    allocate(ek(zahlx,zahlz))   ! force magnitude (L2 norm)    

    allocate(ibc_inside_mask(zahlx,zahlz)) ! true if grid points are completely inside inner boundary
    allocate(ibc_here_mask(zahlx,zahlz)) ! true if grid points are on the inner boundary

! density variable
    rho=0.0   ! density
    rhoo=0.0  ! density
    divrv=0.0
    divrv_1=0.0
    divrv_2=0.0

! pressure variable
    p=0.0     ! pressure
    po=0.0     ! pressure
    dpv=0.0
    divv=0.0
    dp=0.0   ! change in u (timestep)
    dp_1=0.0 ! change in u (timestep)
    dp_2=0.0 ! change in u (timestep)

! velocity and momentum variables
    sx=0.0    ! x-momentum
    sz=0.0    ! z-momentum
    sxo=0.0   ! x-momentum
    szo=0.0   ! z-momentum
    vx=0.0    ! x-velocity
    vz=0.0    ! z-velocity
    vxo=0.0   ! x-velocity
    vzo=0.0   ! z-velocity
    dsx=0.0   ! x-force
    dsz=0.0   ! z-force
    dsx_1=0.0 ! x-force at previous timestep
    dsz_1=0.0 ! z-force at previous timestep
    dsx_2=0.0 ! x-force two timesteps ago
    dsz_2=0.0 ! z-force two timesteps ago
    vel=0.0   ! velocity magnitude

! magnetic field variables
    bx=0.0    ! x-magnetic field
    bz=0.0    ! z-magnetic field
    bxj=0.0   ! Bx minus dipole
    bzj=0.0   ! Bz minus dipole
    bxd=0.0   ! dipolar Bx
    bzd=0.0   ! dipolar Bz
    dbx=0.0   ! change in Bx (timestep)
    dbz=0.0   ! change in Bz (timestep)
    dbx_1=0.0 ! change in Bx (timestep)
    dbz_1=0.0 ! change in Bz (timestep)
    dbx_2=0.0 ! change in Bx (timestep)
    dbz_2=0.0 ! change in Bz (timestep)
    bxo=0.0   ! old Bx (previous iteration)
    bzo=0.0   ! old Bz
    bxjo=0.0  ! old Bx (previous iteration)
    bzjo=0.0  ! old Bz
    divb=0.0  ! divergence of B

    ajx=0.0   ! Jx current
    ajy=0.0   ! Jy current
    ajz=0.0   ! Jz current

! Force variables
    kx=0.0    ! x-comp. of force JxB-grad(P)
    ky=0.0    ! y-comp. of force JxB-grad(P)
    kz=0.0    ! z-comp. of force JxB-grad(P)
    ek=0.0    ! force magnitude (L2 norm)    
   



end subroutine



subroutine deallocate_core_data
use coredata
implicit none

    deallocate(ssx)                  ! grid x-values
    deallocate(ssz)                  ! grid z-values
    deallocate(rx)                   ! central diff d/dx values
    deallocate(rz)                   ! central diff d/dz values
    deallocate(r)
    deallocate(pr)
    deallocate(pr2)
    deallocate(vol)       ! volume around each grid point
    deallocate(xxx)       ! x values at grid points
    deallocate(zzz)       ! z values at grid points
    deallocate(ddx)       ! rx values at grid points
    deallocate(ddz)       ! rz values at grid points

    deallocate(global_mask)
    deallocate(tail_mask)
    deallocate(inner_mask)

! density variable
    deallocate(rho)   ! density
    deallocate(rhoo)  ! density
    deallocate(divrv)
    deallocate(divrv_1)
    deallocate(divrv_2)

! pressure variable
    deallocate(p)    ! pressure
    deallocate(po)    ! pressure
    deallocate(dpv)
    deallocate(divv)
    deallocate(dp)   ! change in u (timestep)
    deallocate(dp_1) ! change in u (timestep)
    deallocate(dp_2) ! change in u (timestep)

! velocity and momentum variables
    deallocate(sx)    ! x-momentum
    deallocate(sz)   ! z-momentum
    deallocate(sxo)  ! x-momentum
    deallocate(szo)   ! z-momentum
    deallocate(vx)    ! x-velocity
    deallocate(vz)    ! z-velocity
    deallocate(vxo)   ! x-velocity
    deallocate(vzo)   ! z-velocity
    deallocate(dsx)   ! x-force
    deallocate(dsz)   ! z-force
    deallocate(dsx_1) ! x-force at previous timestep
    deallocate(dsz_1) ! z-force at previous timestep
    deallocate(dsx_2) ! x-force two timesteps ago
    deallocate(dsz_2) ! z-force two timesteps ago
    deallocate(vel)   ! velocity magnitude

! magnetic field variables
    deallocate(bx)    ! x-magnetic field
    deallocate(bz)    ! z-magnetic field
    deallocate(bxj)   ! Bx minus dipole
    deallocate(bzj)   ! Bz minus dipole
    deallocate(bxd)  ! dipolar Bx
    deallocate(bzd)   ! dipolar Bz
    deallocate(dbx)   ! change in Bx (timestep)
    deallocate(dbz)   ! change in Bz (timestep)
    deallocate(dbx_1) ! change in Bx (timestep)
    deallocate(dbz_1) ! change in Bz (timestep)
    deallocate(dbx_2) ! change in Bx (timestep)
    deallocate(dbz_2) ! change in Bz (timestep)
    deallocate(bxo)   ! old Bx (previous iteration)
    deallocate(bzo)   ! old Bz
    deallocate(bxjo)  ! old Bx (previous iteration)
    deallocate(bzjo)  ! old Bz
    deallocate(divb)  ! divergence of B

    deallocate(ajx)   ! Jx current
    deallocate(ajy)   ! Jy current
    deallocate(ajz)   ! Jz current

! Force variables
    deallocate(kx)    ! x-comp. of force JxB-grad(P)
    deallocate(ky)    ! y-comp. of force JxB-grad(P)
    deallocate(kz)   ! z-comp. of force JxB-grad(P)
    deallocate(ek)   ! force magnitude (L2 norm)    

    deallocate(ibc_inside_mask) ! true if grid points are completely inside inner boundary
    deallocate(ibc_here_mask) ! true if grid points are on the inner boundary


end subroutine
