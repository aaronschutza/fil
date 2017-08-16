!Bicubic interpolation adapted from Numerical Recipes for Fortran 90
!Edited by AMS 2015

SUBROUTINE bcucof(y,y1,y2,y12,d1,d2,c)
!USE nrtype
IMPLICIT NONE
REAL, INTENT(IN) :: d1,d2
REAL, DIMENSION(4), INTENT(IN) :: y,y1,y2,y12
REAL, DIMENSION(4,4), INTENT(OUT) :: c
!Given arrays y,y1,y2, and y12, each of length 4, containing the function, gradients, and
!cross derivative at the four grid points of a rectangular grid cell 
!(numbered counterclockwise from the lower left), and given d1 and d2, the length of the grid 
!cell in the 1- and 2-directions, this routine returns the 4×4 table c that is used by routine
!bcuint for bicubic interpolation.

REAL, DIMENSION(16) :: x
REAL, DIMENSION(16,16) :: wt 
DATA wt /1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,&
8*0,3,0,-9,6,-2,0,6,-4,10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,&
2*0,6,-4,4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,&
2,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2,0,1,-2,1,5*0,&
-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,10*0,-3,3,2*0,2,-2,2*0,&
-1,1,6*0,3,-3,2*0,-2,2,5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,&
-1,0,1,-2,1,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/
x(1:4)=y !Pack a temporary vector x.
x(5:8)=y1*d1
x(9:12)=y2*d2
x(13:16)=y12*d1*d2
x=matmul(wt,x)
!Matrix multiply by the stored table.
c=reshape(x,(/4,4/),order=(/2,1/))
!Unpack the result into the output table.

END SUBROUTINE bcucof

SUBROUTINE bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy,ansy1,ansy2)
!USE nrtype; USE nrutil, ONLY : nrerror
!USE nr, ONLY : bcucof
IMPLICIT NONE
REAL, DIMENSION(4), INTENT(IN) :: y,y1,y2,y12
REAL, INTENT(IN) :: x1l,x1u,x2l,x2u,x1,x2
REAL, INTENT(OUT) :: ansy,ansy1,ansy2
!Bicubic interpolation within a grid square. Input quantities are y,y1,y2,y12(as described in
!bcucof); x1l and x1u, the lower and upper coordinates of the grid square in the 1-direction;
!x2l and x2u likewise for the 2-direction; and x1,x2, the coordinates of the
!desired point for the interpolation. The interpolated function value is returned as ansy,
!and the interpolated gradient values as ansy1 and ansy2. This routine calls bcucof.
INTEGER :: i
REAL :: t,u
REAL, DIMENSION(4,4) :: c
call bcucof(y,y1,y2,y12,x1u-x1l,x2u-x2l,c)
!Get the c’s.
if (x1u == x1l .or. x2u == x2l) write(*,*) 'bcuint: problem with input values - boundary pair equal?'
t=(x1-x1l)/(x1u-x1l)
!Equation (3.6.4).
u=(x2-x2l)/(x2u-x2l)
ansy=0.0
ansy2=0.0
ansy1=0.0
do i=4,1,-1
!Equation (3.6.6).
ansy=t*ansy+((c(i,4)*u+c(i,3))*u+c(i,2))*u+c(i,1)
ansy2=t*ansy2+(3.0*c(i,4)*u+2.0*c(i,3))*u+c(i,2)
ansy1=u*ansy1+(3.0*c(4,i)*t+2.0*c(3,i))*t+c(2,i)
end do
ansy1=ansy1/(x1u-x1l)
ansy2=ansy2/(x2u-x2l)
END SUBROUTINE bcuint


