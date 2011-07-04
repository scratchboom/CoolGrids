C=====  Make step using theoretical solution of ODE system dx/dt = Ax
	SUBROUTINE SOLVEUODES(A,DT,X0,X1)
	USE MSIMSL
	PARAMETER (DBG = 1)
	PARAMETER (N = 4)
	REAL*8 :: A(N,N),DT,X0(N),X1(N)
	INTENT(IN) :: A,DT,X0
	INTENT(OUT) :: X1
	INTEGER :: I,J
	COMPLEX(8) :: EVAL(N), EVEC(N,N), C(N), CX0(N)
	CX0 = X0
	
	if (DBG .eq. 1) then
		CALL DWRRRN('X0',1,N,X0,1,0)
	end if
	
	CALL DEVCRG(N,A,N,EVAL,EVEC,N)
	
	if (DBG .eq. 1) then
		CALL DWRCRN('EVAL',1,N,EVAL,1,0)
	end if
	
	CALL DLSACG(N,EVEC,N,CX0,1,C)
	DO J = 1,N
		C(J) = C(J) * CDEXP(EVAL(J)*DT)
	END DO
	DO I = 1,N
		X1(I) = 0
		DO J = 1,N
			X1(I) = X1(I) + C(J)*EVEC(I,J) 
		END DO
	END DO
c	CALL DWRRRN('X1',1,N,X1,1,0)
	END SUBROUTINE

C=====  Make step using theoretical solution of ODE system dx/dt = Ax + b
	SUBROUTINE SOLVEODES(A,B,DT,X0,X1)
	PARAMETER (DBG = 1)
	PARAMETER (N = 4)
	REAL*8 :: A(N,N),B(N),DT,X0(N),X1(N),Y(N),Z0(N)
	INTENT(IN) :: A,B,DT,X0
	INTENT(OUT) :: X1
	CALL DLSARG(N,A,N,B,1,Y)
	if (DBG .eq. 1) then
		CALL DWRRRN('Y',1,N,Y,1,0)
	end if
	Z0(:) = X0(:) + Y(:)
	CALL SOLVEUODES(A,DT,Z0,X1)
	X1(:) = X1(:) - Y(:)
	END SUBROUTINE

C===================================================
C=====  Right-hand side step
C===================================================
	SUBROUTINE DISSIPATIONSTEP(DT,V0,V1,m)
	PARAMETER (DBG = 1)
	IMPLICIT REAL*8 (a-h,o-z)
	REAL*8 :: DT,V0(7),V1(7),X0(4),X1(4),A(4,4),B(4),D(7,7)
	COMPLEX(8) :: DEVAL(7), DEVEC(7,7)
	INTENT(IN) :: DT,V0
	INTENT(OUT) :: V1
	common/bconst/cs,csv,cze,cme,ckb,ci,cih,ca2,ce0,cedm,cmdm           !obschie bloki dlya rasch. pravoj chasti
	common/bconatm/cmb,cnw,cbv,cno2,czi,czb,cne0,cni0,ct0,cti           !obschie bloki dlya rasch. pravoj chasti
	common/batmosf/zsi,zn,zni0,zne0,zn00,zne,zni,zn0,zno2               !obschie bloki dlya rasch. pravoj chasti
	common/bimpuls/	chi,csi,cpi,cep,z0                                  !obschie bloki dlya rasch. pravoj chasti
	common/bconist/cf1v,cvei,cve0,se00,vei0,ve00,chn0                   !obschie bloki dlya rasch. pravoj chasti
	common/bistoch/tjx,tjz,sie,pesie,vee,pevee,tje,petje,tjeh,petjeh,   
     *tjn,petjn,tjd,petjd,tjp,petjp,vei,pnvei,pevei,ve0,peve0,fie,pefie,
     *chn,pechn,dee,pedee                                                 !obschie bloki dlya rasch. pravoj chasti

	if (DBG .eq. 1) then
		write(*,*) 'm=', m
	end if

	WPI = 3.14159265d0
	WZ = 1d0
	Wevlt = 1.602176487d-19
	WM = cmb ! = 4.84d-26 ! average air ion mass = ??????? in SI
	WEt0 = 0.025d0
	Wmu0 = 4*WPI*1d-7

	Wni = zni ! ions
	Wn0 = zn0 ! neutrals
	Wn = zn ! air
	Wno2 = zno2 ! oxygen

	Wne = V0(1)
	Wvz = V0(2) !! vz
	Wvx = V0(3) !! vx
	WEt = V0(4)
	Wfi = V0(5)
	Wxi = V0(6)
	WEz = V0(7)

	if (DBG .eq. 1) then
		write(*,*) 'ne=', Wne
		write(*,*) 'vz=', Wvz
		write(*,*) 'vx=', Wvx
		write(*,*) 'Et=', WEt
		write(*,*) 'fi=', Wfi
		write(*,*) 'xi=', Wxi
		write(*,*) 'Ez=', WEz
	end if

	WEx = (Wfi+Wxi)/2d0
	WHy = (Wfi-Wxi)/2d0/csv

	Wf_ = 11.67*(1d0 - DEXP(-0.0083*(-11.35d0 + WEt)))*
     &      (1d0 + (0.64d0*(-11.35d0 + WEt))/(88.65d0 + WEt))

	IF (WEt .gt. 14.86d0) THEN
		Wsigma = Wf_*4d0*0.88d-20*(13.6d0/WEt)**2d0
     &             * (WEt - 14.86d0) / 14.86d0
		if (Wsigma .lt. 0d0) then
			Wsigma = 0d0
		end if
	else
		Wsigma = 0d0
	END IF

	Wsigmae0 = 12.47*0.88d-20*(0.4d0+(0.84d0*WEt)/(0.5d0+WEt))


	Wf = 1.5d0 - 0.64d0 + 0.11*Dlog(14.86d0*3d0/2d0/WEt)/2.3026d0

	if (DBG .eq. 1) then
		if (WEt.lt.0d0) then
			WRITE(5,*)'V0= ',V0
		end if
		write(*,*) 'Wj*'
	end if

	Wj0e = (5.3d7*Dsqrt(WEt)*Wsigmae0)                   * 1d-2
	Wjei = (8.75d-27/(2d0/3d0*WEt)**4.5d0)             * 1d-12
	Wjv = (2.7d-13/(2d0/3d0*WEt)**0.75d0)              * 1d-6
	Wjg = (1.16d-8/WEt)                                * 1d-6
	if (WEt .gt. 1e-5) then
		Wjp = ((3.8d-31/WEt)*DEXP(-0.103d0/WEt))        * 1d-12
	else
		Wjp = 0
	end if

	WL = 25.2 + Dlog(2*WEt/3/ckb)- 0.5*Dlog(Wne        * 1d-6)
	if (DBG .eq. 1) then
		write(*,*) 'L<>=', WL
	end if
	
	IF (WL .LT. 15) THEN
		WL = 15
	END IF
	IF (WL .GT. 18) THEN
		WL = 18
	END IF

	Wv = 6d-8*Wn0*Dsqrt(WEt)*(0.4d0+(0.84d0*WEt)/(0.5d0+WEt))
     &                                                   * 1d-6     
	WDelta = (1.7d-3*((1d0+0.2d0*(WEt/0.9d0)**5d0)/
     &         (1+3.7d-2*(1+0.2d0*(WEt/0.9d0)**5d0))))
     &                                                   * 1d-2


	Wvei = (2d0*(cze**4d0)*WL*Wni*Dsqrt(3d0*WPI)*(WZ**2d0))/
     &      (((WEt*Wevlt)**1.5)*Dsqrt(cme))
	Wve0 = 8d0/3d0/Dsqrt(WPI)*Wsigmae0*
     &       (Dsqrt(4d0/3d0*WEt*Wevlt/cme 
     & +      Wvx*Wvx + Wvz*Wvz))*Wn0
	Wve02 = 8d0/3d0/Dsqrt(WPI)*Wsigmae0*
     &       (Dsqrt(4d0/3d0*WEt*Wevlt/cme))*Wn0

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
	WQw_vz = cze*DSIGN(WEz,Wvz)                        / Wevlt
	WQw_vx = cze*DSIGN(WEx,Wvx)                        / Wevlt
	WQw_Et = 0
	WQw_1  = 0

	if (DBG .eq. 1) then
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
     & (1 - DEXP(-0.0083*(-11.35 + WEt)))
	write(*,*) '(1 + (0.64*(-11.35 + WEt))/(88.65 + WEt))=',
     & (1 + (0.64*(-11.35 + WEt))/(88.65 + WEt))
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
	end if

	A(1,1) = Wj0e*Wn0-(Wjg+Wjv)*Wni-Wjp*Wn*Wno2-Wne*Wni*Wjei
	A(1,2) = 0
	A(1,3) = 0
	A(1,4) = 0
	
	A(2,1) = 0
	A(2,2) = -Wve0-Wvei
	A(2,3) = -cze/cme*Wmu0*WHy
	A(2,4) = 0
	
	A(3,1) = 0
	A(3,2) = cze/cme*Wmu0*WHy
	A(3,3) = -Wve0-Wvei
	A(3,4) = 0

	A(4,1) = WSee_ne + WQe_ne + WQw_ne
	A(4,2) = WSee_vz + WQe_vz + WQw_vz
	A(4,3) = WSee_vx + WQe_vx + WQw_vx
	A(4,4) = WSee_Et + WQe_Et + WQw_Et
	

	B(1) = 0
	B(2) = -cze/cme*WEz
	B(3) = -cze/cme*WEx
	B(4) = WSee_1 + WQe_1 + WQw_1

	X0(1) = Wne
	X0(2) = Wvz
	X0(3) = Wvx
	X0(4) = WEt

c	D(1:4,1:4) = A(1:4,1:4)

c	D(1,5) = 0
c	D(1,6) = 0
c	D(1,7) = 0
c	D(2,5) = 0
c	D(2,6) = 0
c	D(2,7) = -cze/cme
c	D(3,5) = -cze/cme/2d0
c	D(3,6) = -cze/cme/2d0
c	D(3,7) = 0
c	D(4,5) = 0
c	D(4,6) = 0
c	D(4,7) = 0

c	D(5,1) = 0
c	D(5,2) = 0
c	D(5,3) = Wne
c	D(5,4) = 0
c	D(6,1) = 0
c	D(6,2) = 0
c	D(6,3) = Wne
c	D(6,4) = 0
c	D(7,1) = 0
c	D(7,2) = Wne
c	D(7,3) = 0
c	D(7,4) = 0

c	CALL DWRRRN('A',4,4,A,4,0)
c	CALL DWRRRN('B',1,4,B,1,0)
c	write(*,*) 'X0=', X0

c	CALL DWRRRN('D',7,7,D,7,0)
c	CALL DEVCRG(7,D,7,DEVAL,DEVEC,7)
c	CALL DWRCRN('DEVAL',1,7,DEVAL,1,0)

	CALL SOLVEODES(A,B,DT,X0,X1)

	V1(1) = X1(1)
	V1(2) = X1(2)
	V1(3) = X1(3)
	V1(4) = X1(4)

	B5 = csv*cs*cze*X1(1)*X1(3)
	B6 = csv*cs*cze*X1(1)*X1(3)
	B7 = csv*cs*cze*X1(1)*X1(2)
	if (DBG .eq. 1) then
		write(*,*) 'B5=', B5
		write(*,*) 'B6=', B6
		write(*,*) 'B7=', B7
	end if

	V1(5) = Wfi + DT*B5
	V1(6) = Wxi + DT*B6
	V1(7) = WEz + DT*B7

	END SUBROUTINE

C===================================================
C=====  Atmosphere
C===================================================
	SUBROUTINE ATMOSPHERE(t,z,ne)
	IMPLICIT REAL*8 (a-h,o-z)
	REAL*8 :: ne
	common/bconst/cs,csv,cze,cme,ckb,ci,cih,ca2,ce0,cedm,cmdm
	common/bconatm/cmb,cnw,cbv,cno2,czi,czb,cne0,cni0,ct0,cti
	common/batmosf/zsi,zn,zni0,zne0,zn00,zne,zni,zn0,zno2
	if(z.lt.czi)then	!Nighnyaya granica ionosfery
	zsi=0.001				!Dolya ionov v ed. obema nevozmusch.atmosfery
	else
	zsi=1.
	end if
	zn=cnw*dexp(cbv*z)	 !Obschee chislo nejtralov i ionov v ed. obema nevozmusch.atmosfery
	zni0=zn*zsi 	   !Chislo ionov v ed. obema nevozmusch.atmosfery
	zne0=zn*zsi*czb    !Chislo elektronov v ed. obema nevozmusch.atmosfery
	zn00=zn*(1.-zsi)   !Chislo nejtralov v ed. obema nevozmusch.atmosfery
	zne=zne0+ne	   !Chislo elektronov v ed. obema
	zni=zni0+ne/czb !Chislo ionov v ed. obema
	zn0=zn00-ne/czb !Chislo nejtralov v ed. obema
	if(zn0.le.0.)then
	zn0=0.
	zni=zn
	zne=zni*czb
	end if
	zno2=cno2*zn0	  !Chislo molekul kisloroda v ed. obema nevozmusch.atmosfery	
	RETURN
	END
