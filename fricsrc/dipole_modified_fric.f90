  SUBROUTINE Dipole_modified (dipole_moment, tilt_angle, x, z, bx,  bz)

     ! calculates gsm components of geodipole field with the dipole moment
     ! for a given dipole moment dipole_moment and tilt angle tilt_angle
     ! Inside of sphere of R=1 Re, we replace the dipole field
     ! with a uniform vertical B-field of strenth 2*Dipole_moment

     ! dipole_moment is in Amps*(Re**2)
     ! tilt_angle is in radian
     ! x,y,z - gsm coordinates in re
     ! bx,bz - field components in gsm system in nanotesla
     ! 2d diple /bei
     use coredata, only: iprec,rprec
    implicit none
    real(rprec), intent(in) :: dipole_moment, tilt_angle, x,  z
    real(rprec), intent(OUT) :: bx,  bz

    real(rprec) :: sps,cps, p,u,v,t,q


!    if (dipole_moment < 20000.0) then
!        write(6,*)' ERROR IN DIP, DM seems too small???'
!        stop
!    end if
    sps=sin(tilt_angle)
    cps=cos(tilt_angle)

  1 p=x**2
    u=z**2
    v=2.*z*x

    if (sqrt(p+u) < 1.0) then ! check for inside earth
        bx = 0.0
        bz = - 2.0*dipole_moment
    else
        q=dipole_moment/(p+u)**2
        bx=q*(v)
        bz=q*(u-p)
    end if

  RETURN
  END SUBROUTINE Dipole_modified
