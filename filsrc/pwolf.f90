function pwolf(r)
      implicit none
      real, intent(in) :: r
      real :: pwolf
      ! fit coefficients
      real, parameter :: rmax = 3.2015
      real, parameter :: a1 = 1.1315
      real, parameter :: b1 = 0.943
      real, parameter :: c1 = 3.0
      real, parameter :: a2 = 2.31
      real, parameter :: b2 = 0.38
      real, parameter :: d1 = 0.141
      real, parameter :: d2 = 0.0412
      real, parameter :: d3 = 3.0
      real, parameter :: d4 = 0.3
            
      ! inner ring current pressure
      if(r <rmax)then
           pwolf = 10**(a1-b1*(r-c1)**2)
      else
           pwolf = 10**(a2-b2*r)
      end if
! outer pressure
      pwolf = pwolf + 10**(d1-d2*r)/(1.0+exp((d3-r)/d4))

      return 
end function pwolf
