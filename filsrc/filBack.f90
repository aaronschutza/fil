!********************************* background functions *****************************************************************!

real function pBack(x,z,elnum) ! pressure of background in (nPa)
       use parametersBackground
       use constants
       implicit none
       real,intent(in) :: x
       real,intent(in) :: z
       real :: term1,coef1,interpolate
    integer :: i,j,elnum

       if(useGrid)then    
              pBack = interpolate('Pg',x,z,elnum)
       else
              if (z > ht .or. z < -ht) then
                     pBack=0.0
              else
                     coef1 = A0*cos(Pi*z/(2.0*ht))*exp(-alpha*x)
                     term1 = Pi/(2.0*ht)
                     pBack = coef1*coef1*(term1*term1-(alpha)**2)/(2.0*mu0)
              end if
       endif 
       
    return
    
end function

real function nBack(x,z,elnum) ! number density of background (cm^-3)
       use parametersBackground
       use constants
       implicit none
       real,intent(in) :: x
       real,intent(in) :: z
       real :: term1,coef1,interpolate
    integer :: i,j,elnum

       if(useGrid)then    
              nBack = interpolate('Ng',x,z,elnum)
       else
              nBack = dBack
       endif 
       
    return
    
end function

real function Bx(x,z,elnum) ! x component of background magnetic field (in nT)
    use parametersBackground
       use constants
       implicit none
       real, intent(in) :: x
       real, intent(in) :: z
       real :: alphaN,lambdaN,coef1,Bn,diff,Aback,interpolate
       integer :: n,sign1,sign2,i,j,elnum

    if(useGrid)then  
              Bx = interpolate('Bxg',x,z,elnum)
       else
        if (z > ht .or. z < -ht) then
                     Bx = -(Pi*A0/(2.0*ht))*cos(alpha*(z - ht))*exp(-alpha*x)
        else
                     Bx = -(Pi*A0/(2.0*ht))*sin(Pi*z/(2.0*ht))*exp(-alpha*x)
        end if
       endif
              
end function

real function Bz(x,z,elnum) ! z component of background magnetic field (in nT)
    use parametersBackground
       use constants
       implicit none
       real, intent(in) :: x
       real, intent(in) :: z
       real :: alphaN,lambdaN,coef1,Bn,Aback,interpolate
       integer :: n,sign1,sign2,i,j,elnum
    
       if(useGrid)then
              Bz = interpolate('Bzg',x,z,elnum)
       else
        if (z > ht .or. z < -ht) then
                     Bz = -(Pi*A0/(2.0*ht))*sin(alpha*(z-ht))*exp(-alpha*x)
        else
                     Bz = alpha*A0*cos(Pi*z/(2.0*ht))*exp(-alpha*x)
        end if
       endif
       
end function

real function ptotal(x,z,elnum)
       use parametersBackground
       use constants
       implicit none
       real,intent(in) :: x
       real,intent(in) :: z
       real :: term1,coef1
    integer :: i,j,elnum
       real :: Bx,Bz,pBack,interpolate

       if(useGrid)then
              ptotal = interpolate('Ptotg',x,z,elnum)
       else
              ptotal = (Bx(x,z,elnum)**2+Bz(x,z,elnum)**2)/(2.0*mu0) + pBack(x,z,elnum)
       endif 
    
end function

real function Aback(x,z,elnum) ! vector potential of background medium (nT*Re)
       use parametersBackground
       use constants
       implicit none
       real :: alphaN,lambdaN,coef1,An
       integer :: n,sign1,sign2,i,j,elnum
       real :: x,z
       real :: interpolate
       
       if(useGrid)then
              Aback = interpolate('Ag',x,z,elnum)
       else
              if ( abs(z) <= ht ) then
                     Aback = -A0*cos(Pi*z/(2.0*ht))*exp(-alpha*x)
              else
                     if ( alpha == 0.0 ) then
                            Aback = Pi*A0*(abs(z)-ht)/(2.0*ht)
                     else
                            Aback = Pi*A0/(2.0*alpha*ht) * sin(alpha*(abs(z)-ht))*exp(-alpha*x)
                     end if
              end if
       endif
       
end function Aback

real function dxBx(x,z,elnum)
       use parametersBackground
       use constants
       implicit none
       real,intent(in) :: x
       real,intent(in) :: z
    integer :: i,j,elnum
       real :: interpolate
       
       if(useGrid)then
              dxBx = interpolate('dx_Bxg',x,z,elnum)
       else
              if (z > ht .or. z < -ht) then
                     dxBx = alpha*(Pi*A0/(2.0*ht))*cos(alpha*(z - ht))*exp(-alpha*x)
        else
                     dxBx = alpha*(Pi*A0/(2.0*ht))*sin(Pi*z/(2.0*ht))*exp(-alpha*x)
        end if
       endif
    
end function

real function dxBz(x,z,elnum)
       use parametersBackground
       use constants
       implicit none
       real,intent(in) :: x
       real,intent(in) :: z
    integer :: elnum
       real :: interpolate

       if(useGrid)then
              dxBz = interpolate('dx_Bzg',x,z,elnum)
       else
              if (z > ht .or. z < -ht) then
                     dxBz = alpha*(Pi*A0/(2.0*ht))*sin(alpha*(z-ht))*exp(-alpha*x)
        else
                     dxBz = -alpha*alpha*A0*cos(Pi*z/(2.0*ht))*exp(-alpha*x)
        end if
       endif 
    
end function

real function dzBx(x,z,elnum)
       use parametersBackground
       use constants
       implicit none
       real,intent(in) :: x
       real,intent(in) :: z
    integer :: elnum
       real :: interpolate

       if(useGrid)then
              dzBx = interpolate('dz_Bxg',x,z,elnum)
       else
              if (z > ht .or. z < -ht) then
                     dzBx = alpha*(Pi*A0/(2.0*ht))*sin(alpha*(z - ht))*exp(-alpha*x)
        else
                     dzBx = -(Pi/(2.0*ht))*(Pi*A0/(2.0*ht))*cos(Pi*z/(2.0*ht))*exp(-alpha*x)
        end if
       endif 
    
end function

real function dzBz(x,z,elnum)
       use parametersBackground
       use constants
       implicit none
       real,intent(in) :: x
       real,intent(in) :: z
    integer :: elnum
       real :: interpolate

       if(useGrid)then
              dzBz = interpolate('dz_Bzg',x,z,elnum)
       else
              if (z > ht .or. z < -ht) then
                     dzBz = -alpha*(Pi*A0/(2.0*ht))*cos(alpha*(z-ht))*exp(-alpha*x)
        else
                     dzBz = -(Pi/(2.0*ht))*alpha*A0*sin(Pi*z/(2.0*ht))*exp(-alpha*x)
        end if
       endif 
    
end function

real function dxPtot(x,z,elnum)
       use parametersBackground
       use constants
       implicit none
       real,intent(in) :: x
       real,intent(in) :: z
    integer :: elnum
       real :: interpolate
       real :: c0,c1,Bx,Bz,dxBx,dxBz

       if(useGrid)then
              dxPtot = interpolate('dx_Ptotg',x,z,elnum)
       else
              dxPtot = 0.0
              if (z > ht .or. z < -ht) then
                     c0=0.0
                     c1=0.0
              else
                     c0 = A0*cos(Pi*z/(2.0*ht))*exp(-alpha*x)
                     c1 = -2*alpha*c0**2*((Pi/(2.0*ht))**2 - (alpha)**2)/(2.0*mu0)
              end if
              dxPtot = c1 + Bx(x,z,elnum)*dxBx(x,z,elnum)/mu0 + Bz(x,z,elnum)*dxBz(x,z,elnum)/mu0
       endif
       
end function

real function dzPtot(x,z,elnum)
       use parametersBackground
       use constants
       implicit none
       real,intent(in) :: x
       real,intent(in) :: z
    integer :: elnum
       real :: interpolate
       real :: c0,c1,Bx,Bz,dzBx,dzBz

       if(useGrid)then
              dzPtot = interpolate('dz_Ptotg',x,z,elnum)
       else
              if (z > ht .or. z < -ht) then
                     c0=0.0
                     c1=0.0
              else
                     c0 = A0*exp(-alpha*x)
                     c1 = -(Pi/ht)*sin(Pi*z/(2.0*ht))*cos(Pi*z/(2.0*ht))*c0**2 &
                     & *((Pi/(2.0*ht))**2 - (alpha)**2)/(2.0*mu0)
              end if
              dzPtot = c1 + Bx(x,z,elnum)*dzBx(x,z,elnum)/mu0 + Bz(x,z,elnum)*dzBz(x,z,elnum)/mu0
       endif
       
end function

real function dxP(x,z,elnum)
       use parametersBackground
       use constants
       implicit none
       real,intent(in) :: x
       real,intent(in) :: z
    integer :: elnum
       real :: interpolate
       real :: c0,c1,Bx,Bz,dxBx,dxBz

       if(useGrid)then
              dxP = interpolate('dx_Pg',x,z,elnum)
       else
              dxP = 0.0
              if (z > ht .or. z < -ht) then
                     c0=0.0
                     c1=0.0
              else
                     c0 = A0*cos(Pi*z/(2.0*ht))*exp(-alpha*x)
                     c1 = -2*alpha*c0**2*((Pi/(2.0*ht))**2 - (alpha)**2)/(2.0*mu0)
              end if
               dxP = c1 
       endif
       
end function

real function dzP(x,z,elnum)
       use parametersBackground
       use constants
       implicit none
       real,intent(in) :: x
       real,intent(in) :: z
    integer :: elnum
       real :: interpolate
       real :: c0,c1,Bx,Bz,dzBx,dzBz

       if(useGrid)then
              dzP = interpolate('dz_Pg',x,z,elnum)
       else
              if (z > ht .or. z < -ht) then
                     c0=0.0
                     c1=0.0
              else
                     c0 = A0*exp(-alpha*x)
                     c1 = -(Pi/ht)*sin(Pi*z/(2.0*ht))*cos(Pi*z/(2.0*ht))*c0**2 &
                     & *((Pi/(2.0*ht))**2 - (alpha)**2)/(2.0*mu0)
              end if
               dzP = c1 
       endif
       
end function


real function dzdx(x,k) ! field line derivatives for analytic bg
       use parametersBackground
       use constants
       implicit none
       real :: x,k
       dzdx = 2*alpha*ht*k*exp(alpha*x)/(Pi*A0*sqrt(1-(k/A0*exp(alpha*x))**2))
end function dzdx

real function dxdz(z) ! field line derivatives for analytic bg
       use parametersBackground
       use constants
       implicit none
       real :: z,arg
       if ( alpha == 0.0 ) then
              arg = -1.0
              dxdz = sqrt(arg)
       else
              dxdz = -Pi/(alpha*2.0*ht)*tan(Pi*z/(2.0*ht))
       end if
end function dxdz

real function baodensity(x,y,z)
       use parametersBackground
       use constants
       real :: x,y,z
       real :: Aback,rho0,A1,DA,MMM

       DA = A0/20.0
       A1 = 208.496
       rho0 = 1.374850299401198
       MMM = 1.6364
       baodensity = rho0*(2.0+MMM*tanh((-Aback(x,z,-1)+A1)/DA))


end function
