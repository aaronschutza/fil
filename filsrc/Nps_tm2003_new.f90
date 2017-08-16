       FUNCTION Nps_tm2003_new (xgsm,ygsm,n_sw,v_sw,imf_bz)
       IMPLICIT NONE
       REAL, INTENT (IN) :: xgsm, ygsm, n_sw, v_sw, imf_bz
!
       REAL :: rho, phi, rho_s, v_sw_s, imf_bs_s, imf_bn_s, n_sw_s,t_ps,n_ps,p_ps,p_ps_s
       REAL :: Nps_tm2003_new,imf_clock,imf_b_perp,p_sw,p_sw_s,F_s
       REAL , DIMENSION (3), PARAMETER :: &
          A1 = (/ 1.678, -0.159 ,  0.057 /), &
          A2 = (/-0.1606, 0.608 ,  0.524 /), &
          A3 = (/1.669,   0.5055,  0.0908/), &
          A4 = (/4.820,   0.0796,  0.527 /), &
          A5 = (/2.855,   0.2746,  0.078 /), &
          A6 = (/-0.602,  0.0361,  -4.422/), &
          A7 = (/-0.836,  -0.0342, -1.533/),&
          A8 = (/-2.491,  -0.7935, -1.217/),&
          A9 = (/0.2568,  1.162  , 2.54  /),&
          A10= (/0.2249,  0.4756 , 0.32  /),&
          A11= (/0.1887,  0.7117 , 0.754 /),&
          A12= (/-0.4458, 0.0    , 1.048 /),&
          A13= (/-0.0331, 0.0    , -0.074/),&
          A14= (/-0.0241, 0.0    , 1.015 /),&
          A15= (/-2.689,  0.0    , 0.0   /),&
          A16= (/1.222,   0.0    , 0.0   /)
!
       t_ps = 0.0
       n_ps = 0.0
       Nps_tm2003_new=0.0
!
       rho = SQRT (xgsm**2 + ygsm**2)
       rho_s = rho / 10.0                  ! in Re
       phi   = -ATAN (ygsm/xgsm)
       p_sw = 1.0E-6*(1.67 + 1.67*4*0.04)*n_sw*v_sw**2  !0.04 helium
       p_sw_s = p_sw / 3.0
       n_sw_s = n_sw / 10.0
       v_sw_s = v_sw / 500.0               ! km/s

       IF (imf_bz < 0.0) THEN         ! assume IMF_By=0
           imf_bs_s = -imf_bz
           imf_bn_s = 0.0
           imf_clock= 3.1415926
       ELSE
           imf_bs_s = 0.0
           imf_bn_s = imf_bz
           imf_clock= 0.0000000
       END IF
       imf_bn_s = imf_bn_s/5.0                ! nT
       imf_bs_s = imf_bs_s/5.0                ! nT
       imf_b_perp=abs(imf_bz)       ! since we assume IMF_By=0
       F_s=imf_b_perp*sqrt(sin(0.5*imf_clock))/5.0
!
       Nps_tm2003_new = (A1(2) + A2(2)*n_sw_s**A10(2) + A3(2)*imf_bn_s + A4(2)*v_sw_s*imf_bs_s)*rho_s**A8(2) + &
                        (A5(2)*n_sw_s**A11(2) + A6(2)*imf_bn_s + A7(2)*v_sw_s*imf_bs_s)*rho_s**A9(2)*(sin(phi)**2)
      
!
       RETURN
       END FUNCTION Nps_tm2003_new

