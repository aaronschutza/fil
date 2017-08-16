SUBROUTINE Dipole_modified_fil (dipole_moment, tilt_angle, x, y, z, bx, by, bz)

     ! calculates gsm components of geodipole field with the dipole moment
     ! for a given dipole moment dipole_moment and tilt angle tilt_angle
     ! Inside of sphere of R=1 Re, we replace the dipole field
     ! with a uniform vertical B-field of strenth 2*Dipole_moment

     ! dipole_moment is in Amps*(Re**2)
     ! tilt_angle is in radian
     ! x,y,z - gsm coordinates in re
     ! bx,by,bz - field components in gsm system in nanotesla

    implicit none
    real, intent(in) :: dipole_moment, tilt_angle, x, y, z
    real, intent(OUT) :: bx, by, bz

    real :: sps,cps, p,u,v,t,q


    if (dipole_moment < 20000.0) then
        write(6,*)' ERROR IN DIP, DM seems too small???'
        stop
    end if
    sps=sin(tilt_angle)
    cps=cos(tilt_angle)

  1 p=x**2
    u=z**2
    v=3.*z*x
    t=y**2

    if (sqrt(p+u+t) < 1.0) then ! check for inside earth
        bx = 0.0
        by = 0.0
        bz = - 2.0*dipole_moment
    else
        q=dipole_moment/sqrt(p+t+u)**5
        bx=q*((t+u-2.0*p)*sps-v*cps)
        by=-3.0*y*q*(x*sps+z*cps)
        bz=q*((p+t-2.0*u)*cps-v*sps)
    end if

  RETURN
END SUBROUTINE Dipole_modified_fil