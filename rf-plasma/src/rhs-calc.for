C==============================================================
C     
C     NOTES:
C     - Use /cpp option in Windows Compaq Fortran
C     
C==============================================================

#ifdef WIN32
#include "win32.conf"
#else
#include "linux.conf"
#endif

C=====Make step using theoretical solution of ODE system dx/dt = Ax


c$$$ commented due to use of euler solver
#if 0
      
#ifdef USE_IMSL
      SUBROUTINE EVALEVEC(A,N,EVAL,EVEC)
      USE MSIMSL
      REAL*8 :: A(N,N)
      COMPLEX(8) :: EVAL(N), EVEC(N,N)
      INTENT(IN) :: A,N
      INTENT(OUT) :: EVAL,EVEC
      CALL DEVCRG(N,A,N,EVAL,EVEC,N)
      END SUBROUTINE
#endif
      
      SUBROUTINE EVALEVEC2(AA,N,EVAL,EVEC)
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 :: AA(N,N)
      COMPLEX(8) :: EVAL(N), EVEC(N,N)
      INTENT(IN) :: AA,N
      INTENT(OUT) :: EVAL,EVEC
      
      e11=AA(1,1)
      e41=AA(4,1)
      e42=AA(4,2)
      e43=AA(4,3)
      e44=AA(4,4)
      a=AA(2,6)
      b1=AA(5,1)
      b2=AA(6,1)
      c=AA(2,2)
      d=AA(3,2)
      write(*,*),' e11=',e11,' e41=',e41,' e42=',e42,' e43=',e43,' e44=',
     1     e44,' a=',a,' b1=',b1,' b2=',b2,' c=',c,' d=',d

      t1 = e11**2
      t4 = c**2
      t5 = d**2
      t6 = t1-2*c*e11+t4+t5
      t8 = e11-e44
      t10 = cmplx(0.D0,-1.D0)
      t11 = t10*d
      t12 = cmplx(0.D0,1.D0)
      t13 = t12*d
      t14 = c+t13-e44
      t16 = c+t11-e44
      t18 = b2*c
      t28 = a*d*e44
      t31 = a*c*e44
      t36 = d*b2
      t43 = d*e43
      t46 = e42*c
      t56 = e11*e41
      t59 = b1*a
      t63 = e43*c
      t66 = e42*d
      t86 = e44*(t5+t4)
      EVEC(1,1) = 0
      EVEC(1,2) = 0
      EVEC(1,3) = 0
      EVEC(1,4) = e11*t6*t8
      EVEC(1,5) = 0
      EVEC(1,6) = 0
      EVEC(2,1) = 0
      EVEC(2,2) = t11*t14
      EVEC(2,3) = t13*t16
      EVEC(2,4) = -(2*t18+d*b1-2*e11*b2)*a*t8/2
      EVEC(2,5) = -t28/2
      EVEC(2,6) = -t31
      EVEC(3,1) = 0
      EVEC(3,2) = -t14*d
      EVEC(3,3) = -t16*d
      EVEC(3,4) = -a*(-b1*e11+b1*c-2*t36)*t8/2
      EVEC(3,5) = -t31/2
      EVEC(3,6) = t28
      EVEC(4,1) = 1
      EVEC(4,2) = -t43-e42*(c+t13)+t46
      EVEC(4,3) = -t43-e42*(c+t11)+t46
      EVEC(4,4) = t1*e11*e41-2*t1*e41*c+t56*t4+t56*t5+t59*e43*e11/2-
     1     t59*t63/2-t59*t66/2-t18*e42*a+b2*e42*a*e11+t36*a*e43
      EVEC(4,5) = a*(t63+t66)/2
      EVEC(4,6) = a*(t46-t43)
      EVEC(5,1) = 0
      EVEC(5,2) = 0
      EVEC(5,3) = 0
      EVEC(5,4) = b1*t6*t8
      EVEC(5,5) = t86
      EVEC(5,6) = 0
      EVEC(6,1) = 0
      EVEC(6,2) = 0
      EVEC(6,3) = 0
      EVEC(6,4) = t6*b2*t8
      EVEC(6,5) = 0
      EVEC(6,6) = t86

      EVAL(1) = e44
      EVAL(2) = c+cmplx(0.D0,1.D0)*d
      EVAL(3) = c+cmplx(0.D0,-1.D0)*d
      EVAL(4) = e11
      EVAL(5) = 0
      EVAL(6) = 0

c$$$  t1 = e44*e41
c$$$  t2 = c**2
c$$$  t6 = e44**2
c$$$  t9 = b1*a
c$$$  t10 = e43*c
c$$$  t12 = c*b2
c$$$  t16 = d**2
c$$$  t26 = e42*d
c$$$  t28 = e43*e44
c$$$  t30 = b2*d
c$$$  t35 = 1/e44
c$$$  t39 = cmplx(0.D0,1.D0)
c$$$  t40 = t39*d
c$$$  t42 = 1/(-c+t40+e44)
c$$$  t45 = 1/(c+t40-e44)
c$$$  t52 = -d*e43+e42*c-e42*e44
c$$$  t55 = t10-t28+t26
c$$$  t60 = t35*t42*t45
c$$$  t65 = cmplx(0.D0,-1.D0/4.D0)
c$$$  t66 = t65*a
c$$$  t67 = cmplx(0.D0,-4.D0)
c$$$  t69 = d*c
c$$$  t73 = b1*d
c$$$  t81 = cmplx(0.D0,2.D0)
c$$$  t84 = cmplx(0.D0,-1.D0)
c$$$  t86 = b1*c
c$$$  t88 = t84*b1
c$$$  t94 = t39*e42
c$$$  t95 = e43+t94
c$$$  t96 = e11**2
c$$$  t97 = c*e11
c$$$  t100 = 1/(t96-2*t97+t2+t16)
c$$$  t103 = 1/(t16+t2)
c$$$  t107 = cmplx(0.D0,-1.D0/2.D0)
c$$$  t112 = t39*c
c$$$  t116 = t95*t103*t45
c$$$  t118 = cmplx(0.D0,1.D0/2.D0)
c$$$  t125 = cmplx(0.D0,-2.D0)
c$$$  t126 = t125*b2
c$$$  t142 = -e43+t94
c$$$  t155 = t103*t42*t142
c$$$  t162 = 1/e11
c$$$  t163 = t162*b2
c$$$  INVEVEC(1,1) = -(-2*t1*t2+4*c*e41*t6+t9*t10+2*t12*e42*a-2*t1*t16
c$$$  &  -2*e41*t6*e44-2*b2*e42*a*e44+t9*t26-t9*t28-2*t30*a*e43)*t35
c$$$  &  /(e11-e44)*t42*t45/2
c$$$  INVEVEC(1,2) = t52*t42*t45
c$$$  INVEVEC(1,3) = t55*t42*t45
c$$$  INVEVEC(1,4) = 1
c$$$  INVEVEC(1,5) = t55*a*t60/2
c$$$  INVEVEC(1,6) = t52*a*t60
c$$$  INVEVEC(2,1) = t66*(t67*b2*t69-2*t12*e11-t73*e11+2*b2*t2-2*b2
c$$$  &  *t16+2*t73*c+t81*e11*t30+t84*e11*t86+t88*t16+t39*b1*t2)*t95
c$$$  &  *t100*t103*t45
c$$$  INVEVEC(2,2) = t107*t45*t95
c$$$  INVEVEC(2,3) = t45*t95/2
c$$$  INVEVEC(2,4) = 0
c$$$  INVEVEC(2,5) = t66*(t112+d)*t116
c$$$  INVEVEC(2,6) = t118*a*(-c+t40)*t116
c$$$  INVEVEC(3,1) = a*(-4*t30*c+t126*t97+t88*d*e11+t81*b2*t2+t126*t16
c$$$  &  +t81*b1*t69+2*t30*e11-t86*e11-b1*t16+b1*t2)*t142*t100*t103
c$$$  &  *t42/4
c$$$  INVEVEC(3,2) = t118*t42*t142
c$$$  INVEVEC(3,3) = t42*t142/2
c$$$  INVEVEC(3,4) = 0
c$$$  INVEVEC(3,5) = a*(c+t40)*t155/4
c$$$  INVEVEC(3,6) = a*(-d+t112)*t155/2
c$$$  INVEVEC(4,1) = t163
c$$$  INVEVEC(4,2) = 0
c$$$  INVEVEC(4,3) = 0
c$$$  INVEVEC(4,4) = 0
c$$$  INVEVEC(4,5) = 0
c$$$  INVEVEC(4,6) = 0
c$$$  INVEVEC(5,1) = -b1*t162
c$$$  INVEVEC(5,2) = 0
c$$$  INVEVEC(5,3) = 0
c$$$  INVEVEC(5,4) = 0
c$$$  INVEVEC(5,5) = 1
c$$$  INVEVEC(5,6) = 0
c$$$  INVEVEC(6,1) = -t163
c$$$  INVEVEC(6,2) = 0
c$$$  INVEVEC(6,3) = 0
c$$$  INVEVEC(6,4) = 0
c$$$  INVEVEC(6,5) = 0
c$$$  INVEVEC(6,6) = 1
      END SUBROUTINE

      SUBROUTINE SOLVEUODES(A,DT,X0,X1)
#ifdef USE_IMSL
      USE MSIMSL
#endif
      PARAMETER (N = 6)
      REAL*8 :: A(N,N),DT,X0(N),X1(N)
      INTENT(IN) :: A,DT,X0
      INTENT(OUT) :: X1
      INTEGER :: I,J
      COMPLEX(8) :: EVAL(N), EVEC(N,N), INVEVEC(N,N), C(N), CX0(N)
      CX0 = X0
      
#ifdef DEBUG_OUTPUT
#ifdef USE_IMSL
      CALL DWRRRN('X0',1,N,X0,1,0)
#endif
#endif
      
#ifdef USE_IMSL 
      CALL EVALEVEC(A,N,EVAL,EVEC)
      CALL DLSACG(N,EVEC,N,CX0,1,C)
#ifdef DEBUG_OUTPUT
      CALL DWRCRN('EVAL',1,N,EVAL,1,0)
      CALL DWRCRN('EVEC',N,N,EVEC,N,0)
#endif
#else
      CALL EVALEVEC2(A,N,EVAL,EVEC)
      CALL ZGETRS()
c$$$  DO I = 1,N
c$$$  C(I) = 0
c$$$  DO J = 1,N
c$$$  C(I) = C(I) + INVEVEC(I,J)*CX0(J)
c$$$  END DO
c$$$  END DO
#endif
      DO J = 1,N
         C(J) = C(J) * CDEXP(EVAL(J)*DT)
      END DO
      DO I = 1,N
         X1(I) = 0
         DO J = 1,N
            X1(I) = X1(I) + C(J)*EVEC(I,J) 
         END DO
      END DO
c     CALL DWRRRN('X1',1,N,X1,1,0)
      END SUBROUTINE

C=====Make step using theoretical solution of ODE system dx/dt = Ax + b
      SUBROUTINE SOLVEODES(A,B,DT,X0,X1)
      PARAMETER (N = 6)
      REAL*8 :: A(N,N),B(N),DT,X0(N),X1(N),Y(N),Z0(N)
      INTENT(IN) :: A,B,DT,X0
      INTENT(OUT) :: X1
#ifdef USE_IMSL
      CALL DLSARG(N,A,N,B,1,Y)
#else
      CALL DGETRS()
#endif

#ifdef DEBUG_OUTPUT
#ifdef USE_IMSL
      CALL DWRRRN('Y',1,N,Y,1,0)
#endif
#endif
      Z0(:) = X0(:) + Y(:)
      CALL SOLVEUODES(A,DT,Z0,X1)
      X1(:) = X1(:) - Y(:)
      END SUBROUTINE

#endif  ! if 0 -- cause use euler solver

      SUBROUTINE EULERSOLVER(A,B,DT,X0,X1)
      PARAMETER (N = 6)
      REAL*8 :: A(N,N),B(N),DT,X0(N),X1(N)
      INTENT(IN) :: A,B,DT,X0
      INTENT(OUT) :: X1
      DO I=1,N
         X1(I)=B(I)
         DO J=1,N
            X1(I) = X1(I) + X0(J)*A(I,J)
         END DO
         X1(I) = X1(I) * DT + X0(I)
      END DO
      END SUBROUTINE

      SUBROUTINE ADAMSSOLVER(F,FP,DT,X0,X1)
      PARAMETER (N = 6)
      REAL*8 :: F(N),FP(N),DT,X0(N),X1(N)
      INTENT(IN) :: F,FP,DT,X0
      INTENT(OUT) :: X1
      X1(:) = X0(:) + DT*(3/2d0*F(:)-1/2d0*FP(:))
      END SUBROUTINE

C===================================================
C=====Right-hand side step
C===================================================
      SUBROUTINE DISSIPATIONSTEP(DT,V0,V1,FP,m,N)
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 :: DT,V0(7),V1(7),X0(6),X1(6),B(6),D(6,6),F(6),FP(6)
      COMPLEX(8) :: DEVAL(6), DEVEC(6,6)
      INTENT(IN) :: DT,V0,m,N
      INTENT(OUT) :: V1
      INTENT(INOUT) :: FP
      common/bconst/cs,csv,cze,cme,ckb,ci,cih,ca2,ce0,cedm,cmdm !obschie bloki dlya rasch. pravoj chasti
      common/bconatm/cmb,cnw,cbv,cno2,czi,czb,cne0,cni0,ct0,cti !obschie bloki dlya rasch. pravoj chasti
      common/batmosf/zsi,zn,zni0,zne0,zn00,zne,zni,zn0,zno2 !obschie bloki dlya rasch. pravoj chasti
      common/bimpuls/ chi,csi,cpi,cep,z0 !obschie bloki dlya rasch. pravoj chasti
      common/bconist/cf1v,cvei,cve0,se00,vei0,ve00,chn0 !obschie bloki dlya rasch. pravoj chasti
      common/bistoch/tjx,tjz,sie,pesie,vee,pevee,tje,petje,tjeh,petjeh,
     1     tjn,petjn,tjd,petjd,tjp,petjp,vei,pnvei,pevei,ve0,peve0,fie,
     2     pefie,chn,pechn,dee,pedee !obschie bloki dlya rasch. pravoj chasti
      
#ifdef DEBUG_OUTPUT
      IF(m .eq. 1) write(*,*) 'nextstep=', N
      write(*,*) 'm=', m
#endif

      WPI = 3.14159265d0
      WZ = 1d0
      Wevlt = cze
      WM = cmb                  ! = 4.84d-26 ! average air ion mass = ??????? in SI
      WEt0 = 0.025d0
      Wmu0 = 4*WPI*1d-7

      Wni = zni                 ! ions
      Wn0 = zn0                 ! neutrals
      Wn = zn                   ! air
      Wno2 = zno2               ! oxygen

#ifdef DEBUG_OUTPUT
      write(*,*) 'ni=', Wni
      write(*,*) 'n0=', Wn0
      write(*,*) 'n=', Wn
      write(*,*) 'no2=', Wno2
#endif

      Wne = V0(1)
      Wvz = V0(2)               !! vz
      Wvx = V0(3)               !! vx
      WEt = V0(4)
      Wfi = V0(5)
      Wxi = V0(6)
      WEz = V0(7)

      WEx = (Wfi+Wxi)/2d0
      WHy = (Wfi-Wxi)/2d0/csv

#ifdef DEBUG_OUTPUT
      write(*,*) 'ne=', Wne
      write(*,*) 'vz=', Wvz
      write(*,*) 'vx=', Wvx
      write(*,*) 'Et=', WEt
      write(*,*) 'fi=', Wfi
      write(*,*) 'xi=', Wxi
      write(*,*) 'Ex=', WEx
      write(*,*) 'Hy=', WHy
      write(*,*) 'Ez=', WEz
#endif


      Wf_ = 11.67*(1d0 - DEXP(-0.0083*(-11.35d0 + WEt)))*
     1     (1d0 + (0.64d0*(-11.35d0 + WEt))/(88.65d0 + WEt))

      IF (WEt .gt. 14.86d0) THEN
         Wsigma = Wf_*4d0*0.88d-20*(13.6d0/WEt)**2d0
     1        * (WEt - 14.86d0) / 14.86d0
         if (Wsigma .lt. 0d0) then
            Wsigma = 0d0
         end if
      else
         Wsigma = 0d0
      END IF

      Wsigmae0 = 12.47*0.88d-20*(0.4d0+(0.84d0*WEt)/(0.5d0+WEt))


      Wf = 1.5d0 - 0.64d0 + 0.11*Dlog(14.86d0*3d0/2d0/WEt)/2.3026d0

#ifdef DEBUG_OUTPUT
      if (WEt.lt.0d0) then
         WRITE(5,*)'V0= ',V0
      end if
      write(*,*) 'Wj*'
#endif

      Wj0e = (5.3d7*Dsqrt(WEt)*Wsigma)                   * 1d-2
      Wjei = (8.75d-27/(2d0/3d0*WEt)**4.5d0)             * 1d-12
      Wjv = (2.7d-13/(2d0/3d0*WEt)**0.75d0)              * 1d-6
      Wjg = (1.16d-8/WEt)                                * 1d-6
      if (WEt .gt. 1e-5) then
         Wjp = ((3.8d-31/WEt)*DEXP(-0.103d0/WEt))        * 1d-12
      else
         Wjp = 0
      end if

      WL = 25.2 + Dlog(2*WEt/3/ckb)- 0.5*Dlog(Wne        * 1d-6)
#ifdef DEBUG_OUTPUT
      write(*,*) 'L<>=', WL
#endif
      
      IF (WL .LT. 15) THEN
         WL = 15
      END IF
      IF (WL .GT. 18) THEN
         WL = 18
      END IF

      Wv = 6d-8*Wn0*Dsqrt(WEt)*(0.4d0+(0.84d0*WEt)/(0.5d0+WEt))
     1     * 1d-6     
      WDelta = (1.7d-3*((1d0+0.2d0*(WEt/0.9d0)**5d0)/
     1     (1+3.7d-2*(1+0.2d0*(WEt/0.9d0)**5d0))))
     2     * 1d-2


      Wvei = (2d0*(cze**4d0)*WL*Wni*Dsqrt(3d0*WPI)*(WZ**2d0))/
     1     (((WEt*Wevlt)**1.5)*Dsqrt(cme))
      Wve0 = 8d0/3d0/Dsqrt(WPI)*Wsigmae0*
     1     (Dsqrt(4d0/3d0*WEt*Wevlt/cme 
     2     +      Wvx*Wvx + Wvz*Wvz))*Wn0
      Wve02 = 8d0/3d0/Dsqrt(WPI)*Wsigmae0*
     1     (Dsqrt(4d0/3d0*WEt*Wevlt/cme))*Wn0

      WSee_ne = (14.86d0+WEt)*Wni*Wjei
      WSee_vz = 0
      WSee_vx = 0
      WSee_Et = -Wn0*Wj0e + Wf*2d0/3d0*Wni*Wjv - Wjg*Wni - Wjp*Wno2*Wn
      WSee_1  = -14.86d0*Wn0*Wj0e

      WQe_ne = 0
      WQe_vz = (-Wvei*cme*Wvz)                           / Wevlt
      WQe_vx = (-Wvei*cme*Wvx)                           / Wevlt
      WQe_Et = -2d0*Wvei*cme/WM - Wv*WDelta
      WQe_1  = 2d0*Wvei*cme*WEt0/WM + WEt0*Wv*WDelta

      WQw_ne = 0
      WQw_vz = DSIGN(WEz,Wvz)                  ! *cze/Wevlt == 1
      WQw_vx = DSIGN(WEx,Wvx)                  ! *cze/Wevlt == 1
c      WQw_vz = -cze*WEz                        / Wevlt
c      WQw_vx = -cze*WEx                        / Wevlt
      WQw_Et = 0
      WQw_1  = 0

c      write(555,*)WQe_vz*Wvz,WQe_vx*Wvx,WQe_Et*WEt,WQe_1,WSee_ne*Wne,
c     1     WSee_Et*WEt,WSee_1

#ifdef DEBUG_OUTPUT
      Wvv = dsqrt(cze*dabs(WEx)*2d0/(Wn0*Wsigmae0*cme))
      WvEt = dsqrt(WEt*Wevlt/cme)
      write(*,*) 'vv=', Wvv
      write(*,*) 'vEt=', WvEt
      write(*,*) 'j0e=', Wj0e
      write(*,*) 'jg=', Wjg
      write(*,*) 'jv=', Wjv
      write(*,*) 'jp=', Wjp
      write(*,*) 'jei=', Wjei
      write(*,*) 'L=', WL
      write(*,*) 'f=', Wf
      write(*,*) 'sigma=', Wsigma
      write(*,*) 'sigmae0=', Wsigmae0
      write(*,*) 'f_=', Wf_
      write(*,*) '(1 - DEXP(-0.0083*(-11.35 + WEt)))=',
     1     (1 - DEXP(-0.0083*(-11.35 + WEt)))
      write(*,*) '(1 + (0.64*(-11.35 + WEt))/(88.65 + WEt))=',
     1     (1 + (0.64*(-11.35 + WEt))/(88.65 + WEt))
      write(*,*) 'vei=', Wvei
      write(*,*) 've0=', Wve0
      write(*,*) 've02=', Wve02
      write(*,*) 'v=', Wv

      write(*,*) 'See_ne=', WSee_ne
      write(*,*) 'See_vz=', WSee_vz
      write(*,*) 'See_vx=', WSee_vx
      write(*,*) 'See_Et=', WSee_Et
      write(*,*) 'See_1=', WSee_1

      write(*,*) 'Qe_ne=', WQe_ne
      write(*,*) 'Qe_vz=', WQe_vz
      write(*,*) 'Qe_vx=', WQe_vx
      write(*,*) 'Qe_Et=', WQe_Et
      write(*,*) 'Qe_1=', WQe_1

      write(*,*) 'Qw_ne=', WQw_ne
      write(*,*) 'Qw_vz=', WQw_vz
      write(*,*) 'Qw_vx=', WQw_vx
      write(*,*) 'Qw_Et=', WQw_Et
      write(*,*) 'Qw_1=', WQw_1

      WSe = (Wj0e*Wn0*Wne - Wne*Wne*Wni*Wjei)
     1     -(Wjg+Wjv)*Wni*Wne - Wne*Wjp*Wn*Wno2

      write(*,*) 'Se=', WSe
      write(*,*) 'WSee_1 + WQe_1 + WQw_1=', WSee_1 + WQe_1 + WQw_1

#endif

      D(1,1) = Wj0e*Wn0-(Wjg+Wjv)*Wni-Wjp*Wn*Wno2-Wne*Wni*Wjei
      D(1,2) = 0
      D(1,3) = 0
      D(1,4) = 0
      
      D(2,1) = 0
      D(2,2) = -Wve0-Wvei
      D(2,3) = -cze*Wmu0*WHy/cme
      D(2,4) = 0
      
      D(3,1) = 0
      D(3,2) = cze*Wmu0*WHy/cme
      D(3,3) = -Wve0-Wvei
      D(3,4) = 0

      D(4,1) = WSee_ne + WQe_ne + WQw_ne
      D(4,2) = WSee_vz + WQe_vz + WQw_vz
      D(4,3) = WSee_vx + WQe_vx + WQw_vx
      D(4,4) = WSee_Et + WQe_Et + WQw_Et
      
      D(1,5) = 0
      D(1,6) = 0
      D(2,5) = 0
      D(2,6) = -cze/cme       ! Почему тут нет деления на 2.0
      D(3,5) = -cze/cme/2d0   ! А тут есть???
      D(3,6) = 0
      D(4,5) = 0
      D(4,6) = 0

      D(5,1) = 0
      D(5,2) = 0
      D(5,3) = csv*cs*cze*Wne
      D(5,4) = 0
      D(6,1) = 0
      D(6,2) = csv*cs*cze*Wne
      D(6,3) = 0
      D(6,4) = 0

      D(5:6,5:6) = 0

      B(1) = 0
      B(2) = 0
      B(3) = 0
      B(4) = WSee_1 + WQe_1 + WQw_1
      B(5) = 0
      B(6) = 0


#ifdef DEBUG_OUTPUT
c     CALL DWRRRN('A',4,4,A,4,0)
c     write(*,*) 'X0=', X0

#ifdef USE_IMSL
      CALL DWRRRN('D',6,6,D,6,0)
      CALL DWRRRN('B',1,6,B,1,0)
#endif
#endif

c     CALL DEVCRG(6,D,6,DEVAL,DEVEC,6)
c#ifdef DEBUG_OUTPUT
c     CALL DWRCRN('DEVAL',1,6,DEVAL,1,0)
c#endif

      X0(1:4) = V0(1:4)
      X0(5) = V0(5) + V0(6)
      X0(6) = V0(7)

c     CALL SOLVEODES(D,B,DT,X0,X1)
c      CALL EULERSOLVER(D,B,DT,X0,X1)

      DO I=1,6
         F(I)=B(I)
         DO J=1,6
            F(I) = F(I) + X0(J)*D(I,J)
         END DO
      END DO
	IF (N .LE. 1) THEN
		CALL EULERSOLVER(D,B,DT,X0,X1)
	ELSE
		CALL ADAMSSOLVER(F,FP,DT,X0,X1)
	ENDIF
	FP(:)=F(:)

      V1(1:4) = X1(1:4)
      V1(5) = V0(5) + (X1(5) - X0(5))/2d0
      V1(6) = V0(6) + (X1(5) - X0(5))/2d0
      V1(7) = X1(6)

#ifdef DEBUG_OUTPUT
      write(*,*) 'delta1=', V1(1)-V0(1)
      write(*,*) 'delta2=', V1(2)-V0(2)
      write(*,*) 'delta3=', V1(3)-V0(3)
      write(*,*) 'delta4=', V1(4)-V0(4)
      write(*,*) 'delta5=', V1(5)-V0(5)
      write(*,*) 'delta6=', V1(6)-V0(6)
      write(*,*) 'delta7=', V1(7)-V0(7)
#endif

c     B5 = csv*cs*cze*X1(1)*X1(3)
c     B6 = csv*cs*cze*X1(1)*X1(3)
c     B7 = csv*cs*cze*X1(1)*X1(2)
c     if (DBG .eq. 1) then
c     write(*,*) 'B5=', B5
c     write(*,*) 'B6=', B6
c     write(*,*) 'B7=', B7
c     end if

c     V1(5) = Wfi + DT*B5
c     V1(6) = Wxi + DT*B6
c     V1(7) = WEz + DT*B7

c      if (m .eq. 1) then
c         write (666,*) 'Et(prch)=', V1(4) - V0(4)
c      end if

      END SUBROUTINE

C===================================================
C=====Atmosphere
C===================================================
      SUBROUTINE ATMOSPHERE(t,z,ne)
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 :: ne
      common/bconst/cs,csv,cze,cme,ckb,ci,cih,ca2,ce0,cedm,cmdm
      common/bconatm/cmb,cnw,cbv,cno2,czi,czb,cne0,cni0,ct0,cti
      common/batmosf/zsi,zn,zni0,zne0,zn00,zne,zni,zn0,zno2
      if(z.lt.czi)then          !Nighnyaya granica ionosfery
         zsi=0.001              !Dolya ionov v ed. obema nevozmusch.atmosfery
      else
         zsi=1.
      end if
      zn=cnw*dexp(cbv*z)        !Obschee chislo nejtralov i ionov v ed. obema nevozmusch.atmosfery
      zni0=zn*zsi               !Chislo ionov v ed. obema nevozmusch.atmosfery
      zne0=zn*zsi*czb           !Chislo elektronov v ed. obema nevozmusch.atmosfery
      zn00=zn*(1.-zsi)          !Chislo nejtralov v ed. obema nevozmusch.atmosfery
      zne=zne0+ne               !Chislo elektronov v ed. obema
      zni=ne                    !zni0+ne/czb !Chislo ionov v ed. obema
      zn0=zn00-ne/czb           !Chislo nejtralov v ed. obema
      if(zn0.le.0.)then
         zn0=0.
         zni=zn
         zne=zni*czb
      end if
      zno2=cno2*zn0             !Chislo molekul kisloroda v ed. obema nevozmusch.atmosfery    
      RETURN
      END
