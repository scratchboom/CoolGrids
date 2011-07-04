#ifdef WIN32
#include "win32.conf"
#else
#include "linux.conf"
#endif

c#define RHS_DISABLE
c#define STEP_PAUSE
c#define FPOT_OUT
c#define AFTER_OUT

	PROGRAM MAIN
C -----------------------------------------------------------------------
c	This routine solves 1D CFD problems in the form dU/dt + dF(U)/dx = 0
C------------------------------------------------------------------------
#ifdef USE_AV
	use AVDef  ! AVDef is the AViz module file
	use AVViewer
	use DFLib
#endif
C-----------------------------------------------------------------------
	IMPLICIT REAL*8 (a-h,o-z)
	INTEGER status,hv
	CHARACTER(1) :: key,key1
	CHARACTER(5) :: var,Chm
	CHARACTER(4096) :: ArrayName
	PARAMETER (IU = 3)
	PARAMETER (IUU = 4)
	PARAMETER (IK = 7)
	LOGICAL PRINTMORE
C-----------------------------------------------------------------------------------------------------------------
	COMMON /MIN/ Rmin,Emin,Cmin,Pmin
	common/bconst/cs,csv,cze,cme,ckb,ci,cih,ca2,ce0,cedm,cmdm			!obschie bloki dlya rasch. pravoj chasti
	common/bconatm/cmb,cnw,cbv,cno2,czi,czb,cne0,cni0,ct0,cti			!obschie bloki dlya rasch. pravoj chasti
	common/batmosf/zsi,zn,zni0,zne0,zn00,zne,zni,zn0,zno2				!obschie bloki dlya rasch. pravoj chasti
	common/bimpuls/ chi,csi,cpi,cep,z0									!obschie bloki dlya rasch. pravoj chasti
	common/bconist/cf1v,cvei,cve0,se00,vei0,ve00,chn0					!obschie bloki dlya rasch. pravoj chasti
	common/bistoch/tjx,tjz,sie,pesie,vee,pevee,tje,petje,tjeh,petjeh,	
     *tjn,petjn,tjd,petjd,tjp,petjp,vei,pnvei,pevei,ve0,peve0,fie,pefie,
     *chn,pechn,dee,pedee 												!obschie bloki dlya rasch. pravoj chasti
C-----------------------------------------------------------------------------------------------------------------
	ALLOCATABLE :: W(:),W1(:),V0(:),UN1(:),UN1U(:),UN2(:),UN2U(:),
     &               UB(:),S(:),FF(:),
     &			   UP1(:),UP2(:),UP1U(:),UP2U(:),FM(:),FMU(:),B(:),
     &               CC(:),GIB(:),Chm(:)
	ALLOCATABLE :: X(:),ALM(:), RUS(:),  OM(:,:),OMO(:,:)
	ALLOCATABLE :: VF(:,:),V(:,:),V1(:,:),UN(:,:),UP(:,:),FV(:,:)
     &,vs(:,:),FP(:,:),FPtmp(:)															!massiv dlya pechaty (vs)
	ALLOCATABLE :: VM1(:),VN(:),VN1(:),yn(:),y(:),amn(:),amx(:),		
     &vm2(:),F(:),R(:),F1(:),GPR(:),FT(:),FT1(:),fmn(:),fmx(:)	        !massivy dlya pravoj chasti (yn,y,F,R,F1,G,FT,FT1) 
                                                                          ! pechaty (amn,amx). Vspomog. massiv VM2
	ALLOCATABLE :: FTY(:,:),FYYI(:,:),RF(:,:),FY(:,:),FY1(:,:)			!massivy  dlya rasch. pravoj chasti FTY,FYYI,RF,FY,FY1
	ALLOCATABLE :: MKP(:),LKP(:)										!massivy  dlya rasch. pravoj chasti MK ,LK	 
C-----------------------------------------------------------------
C	Reading of constants
C--------------------------------------------------------------------------------------------------------------------------------------------------------------------

c	call vvodconstprch(J)												!vvod const pravoj chasti iz fajla bconst
	OPEN(77,file="befor-hyp.txt",status="replace",iostat=istat)
	PRINTMORE=.FALSE.
#ifdef AFTER_OUT
	OPEN(78,file="after-hyp.txt",status="replace",iostat=istat)
	OPEN(79,file="after-rhs.txt",status="replace",iostat=istat)
#endif
#ifdef FPOT_OUT
	OPEN(80,file="fpot.txt",status="replace",iostat=istat)
#endif
	IF(istat /= 0) STOP "***Cannot open file***"
	OPEN(UNIT=42,FILE='bconst.cfg',STATUS='UNKNOWN',FORM='FORMATTED')
	REWIND 42
	OPEN(UNIT=43,FILE='parskhgs.cfg',STATUS='UNKNOWN',FORM='FORMATTED')
	REWIND 43
	OPEN(UNIT=44,FILE='bconatm.cfg',STATUS='UNKNOWN',FORM='FORMATTED')
	REWIND 44
	OPEN(UNIT=45,FILE='bimpuls.cfg',STATUS='UNKNOWN',FORM='FORMATTED')
	REWIND 45
	OPEN(UNIT=46,FILE='bconist.cfg',STATUS='UNKNOWN',FORM='FORMATTED')
	REWIND 46
CC	Vvod constant iz fajla BCONST
	WRITE(*,*)'Skorost sveta (m/s) cs=?'
	READ(42,*)cs
	WRITE(*,*)'cs=',cs
	WRITE(*,*)'Soprotivlenie vakuuma (s/m?) csv=?'
	READ(42,*)csv
	WRITE(*,*)'csv=',csv
	WRITE(*,*)'Zaryad elektrona (?) cze=?'
	READ(42,*)cze
	WRITE(*,*)'cze=',cze
	WRITE(*,*)'Massa elektrona (?) cme=?'
	READ(42,*)cme
	WRITE(*,*)'cme=',cme
	WRITE(*,*)'Konstanta Bolcmana (?) ckb=?'
	READ(42,*)ckb
	WRITE(*,*)'ckb=',ckb 
	WRITE(*,*)'Konstanta I v f-le secheniya ionizacii (eV) ci=?'
	READ(42,*)ci
	WRITE(*,*)'ci=',ci 
	WRITE(*,*)'Konstanta Ih v f-le secheniya ionizacii (eV) ci=?'
	READ(42,*)cih
	WRITE(*,*)'cih=',cih 
	WRITE(*,*)'Konstanta alfa2 v f-le secheniya ionizacii (sm2) ca2=?'
	READ(42,*)ca2
	WRITE(*,*)'ca2=',ca2 
	WRITE(*,*)'Konstanta e0 v f-le dlya Qei (eV) ce0=?'
	READ(42,*)ce0
	WRITE(*,*)'ce0=',ce0 
C	Vvod parametrov atmosfery z fajla BCONST
	WRITE(*,*)'massa iona vozdukha(g?) cmb=?'
	READ(42,*)cmb
	WRITE(*,*)'cmb=',cmb 
	WRITE(*,*)'Konzentraciya vozdukha u poverkhnosti zemli(1/sm3) cnw=?'
	READ(42,*)cnw
	WRITE(*,*)'cnw=',cnw 
	WRITE(*,*)'Pokazatel exponenty v raspr.konc.vozd.po vysote (?) cbv=?'
	READ(42,*)cbv
	WRITE(*,*)'cbv=',cbv 
	WRITE(*,*)'Dolya kisloroda (0-1) cno2=?'
	READ(42,*)cno2
	WRITE(*,*)'cno2=',cno2 
	WRITE(*,*)'Vysota nachala ionosfery (m) czi=?'
	READ(42,*)czi
	WRITE(*,*)'czi=',czi 
	WRITE(*,*)'Stepen ionizacii (1-...) czb=?'
	READ(42,*)czb
	WRITE(*,*)'czb=',czb 
	WRITE(*,*)'Nachalnaya koncentraciya elektronov (1/sm3) cne0=?'
	READ(42,*)cne0
	WRITE(*,*)'cne0=',cne0 
	WRITE(*,*)'Nachalnaya koncentraciya ionov (1/sm3) cni0=?'
	READ(42,*)cni0
	WRITE(*,*)'cni0=',cni0 
	WRITE(*,*)'Temperatura vozdukha v atmocfere (grK) ct0=?'
	READ(42,*)ct0
	WRITE(*,*)'ct0=',ct0 
	WRITE(*,*)'Temperatura vozdukha v ionosfere (grK) cti=?'
	READ(42,*)cti
	WRITE(*,*)'cti=',cti
C	Vvod parametrov impulsa iz fajla BCONST
	WRITE(*,*)'Opornaya chastota (Gc) chi=?'
	READ(42,*)chi
	WRITE(*,*)'chi=',chi
	WRITE(*,*)'Chastota sledovaniya impulsov (Gc) csi=?'
	READ(42,*)csi
	WRITE(*,*)'csi=',csi
	WRITE(*,*)'Chislo impulsov v pachke ( cpi=?'
	READ(42,*)cpi
	WRITE(*,*)'cpi=',cpi
	WRITE(*,*)'Amplituda electricheskogo polya (V/m) cep=?'
	READ(42,*)cep
	WRITE(*,*)'cep=',cep
C	Vvod parametrov skhemy i nachalnykh dannykh
	IKK=20 !!!100
	EK=0.0001
	WRITE(*,*)'Nachalnyj shag integrirowanija pravoj chasti DTPR=?'
	READ(42,*)DTPR
	WRITE(*,*)'dtpr=',DTPR
	WRITE(*,*)'Nachalnaya vysota Z0=?'
	READ(42,*)Z0
	WRITE(*,*)'Z0=',Z0
	READ(46,*)cf1v
	WRITE(*,*)'cf1v=',cf1v
	READ(46,*)cvei
	WRITE(*,*)'cvei=',cvei
	READ(46,*)cve0
	WRITE(*,*)'cve0=',cve0
	READ(42,*)cedm
	WRITE(*,*)'cedm=',cedm
	READ(42,*)cmdm
	WRITE(*,*)'cmdm=',cmdm
c	WRITE(*,*)cs,csv,cze,cme,ckb,ci,cih,ca2,ce0,cedm,cmdm
c	WRITE(*,*)chi,csi,cpi,cep,z0  
c	WRITE(*,*)cmb,cnw,cbv,cno2,czi,czb,cne0,cni0,ct0,cti
	C00=0.																 ! const skhemy pravoj chasti
	C10=0.																 ! const skhemy pravoj chasti
	A0=-1.																 ! const skhemy pravoj chasti
	B00=0. !0.5+C00+C10 													 ! const skhemy pravoj chasti
	B10=1. !0.5-C00-C10 													 ! const skhemy pravoj chasti
C--------------------------------------------------------------------------------------------------------
	OPEN(1,file="Data1D.in",status="old",iostat=istat)
	IF(istat /=0) STOP "***Cannot open Data1D.in***"
	REWIND 1
	READ(1,*)MK 	!Number of the grid steps in one dimension
	READ(1,*)G		!Addiabatic constant (= 5/3 by default for our case!)
	READ(1,*)XXK		!Region size
	READ(1,*)NL 	!Type of Boundary Conditions on the left side (NL = 1 by default for our case!)
	READ(1,*)NR 	!Type of Boundary Conditions on the right side (NR = 1 by default for our case!)
	READ(1,*)AKYR		!Currant number
	READ(1,*)BG 	!Order of Scheme accuracy (1, 2 or 3)
	READ(1,*)RID		!Parameter of the discontinuity (0.001-0.01)
	READ(1,*)Rmin	  !Minimal Electron Concentration - Rmin
	Rmin=Rmin*cme !!!!
	IF(Rmin <= 0.) STOP "Programm can't work then Rmin <= 0"
	READ(1,*)Emin	  !Minimal Electron Internal Energy - Emin
	Pmin=(G-1.0)*Emin*Rmin
	IF(Emin <= 0.) STOP "Programm can't work then Emin <= 0"
	READ(1,*)EPS	  !Parameter of the Scheme Hybridity  (0.001-0.01)
	READ(1,*)Elch	  !Electron Charge
	READ(1,*)Vlt	  !Velocity of Light
	READ(1,*)VcIm	  !Impedance of Vacuum
	write(*,*) 'VcIm=', VcIm
	READ(1,*)
	ALLOCATE(X(MK),ALM(MK),W(IUU),W1(IU),UN1(IU),UN2(IU),UN1U(IUU),
     &		   UN2U(IUU),UB(IU),UP1U(IUU),UP2U(IUU),
     &		   V0(IUU),FM(IU),FMU(IUU),B(IU),CC(IU),GIB(IU),UP1(IU),UP2(IU),
     &		   UN(MK,IU),UP(MK,IU),V(MK,IUU),V1(MK,IUU),RUS(IU),Chm(MK),
     &		   FF(IU),OM(IU,IU),OMO(IU,IU),VF(MK,IK),FV(MK,IK),S(IU),
     &		   vs(mk,ik),yn(ik),y(ik),amn(ik),amx(ik),VM1(MK),VN(MK),
     &		   F(ik),R(ik),F1(ik),GPR(ik),FT(ik),FT1(ik),
     &		   FTY(ik,ik),FYYI(ik,ik),RF(ik,ik),FY(ik,ik),FY1(ik,ik),
     &		   vm2(mk),VN1(MK),MKP(ik),LKP(ik),fmn(ik),fmx(ik),
     &           FP(mk,6),FPtmp(6),
     &		   stat=istat)										 !massivy dlya pravoj chasti (yn,y,F,R,F1,G,FT,FT1) i pechaty (vs,amn,amx). Vspomog. massiv VM2
	IF(istat /=0) STOP "***ERROR ON ALLOCATE***"
c	write(*,*)mk,g,xxk,nl,nr,akyr,bg,rid,rmin,emin,eps,elch,vlt,vcim
C -------------------------------------------------------
C	Input of the initial data
C -------------------------------------------------------
	WRITE(*,*)'variant novyj (0) ili prodolghenie starogo (1) nvar=?'
	READ(*,*)nvar
	WRITE(*,*)'nvar=',nvar
	DO J = 1, IK
		READ(1,*) (VF(M,J), M=1,MK) !Data arrays of the Concentration, VelocityX,VelocityZ,Internal Energy,Ex,Hy,Ez
		READ(1,*)
	END DO
	CLOSE(1)
C -------------- TEST initial data -----------------------
c	vf(:,1) = 100
c	vf(:,4) = 2
c	vf(100:200,1) = 200
c	vf(100:200,4) = 1
c	vf(:,3) = 1d8
c	vf(150:250,2) = 12345
C --------------- Field test -----------------------------
c	vf(:,5) = 1d7
c	vf(:,6) = 1d7
c	vf(:,7) = 0

C --------------------------------------------------------
	if(nvar.eq.1)then
		OPEN(UNIT=5,FILE='StateDump',STATUS='OLD',iostat=istat)
		IF(istat /=0) STOP "***Cannot open StateDump***"
		REWIND 5
		read(5,*)N,T,ZT
		read(5,*) X
		read(5,*) VF
		read(5,*)
		CLOSE(UNIT=5)
	end if
C -------------------------------------------------------
	CALL URSOST(Rmin,Emin,G,Cmin,Pmin,PR,PE)
	WRITE(*,*)'Ending time(sec) TK=?'
	READ(*,*)TK
	OPEN(1,file="Data1D.in",status="old",iostat=istat)
	IF(istat /=0) STOP "***Cannot open Data1D.in***"
	WRITE(*,*)' Output saving step (number steps) NPK=?'
	READ(*,*)NPK
	WRITE(*,*)' Output visualization step (number steps) NPG=?'
	READ(*,*)NPG
	WRITE(*,*)' Would you like to move the region boundaries? (y/n)?'
	READ(*,*)key1
	IF(key1 == "y")THEN
		WRITE(*,*)'Boundary flow disturbance checking parameter DELTA=?'
		WRITE(*,*)'(0.000001-0.001)'
		READ(*,*)DELTA
	END IF
	if(nvar.eq.0)then
		N=0
		T = 0.
	end if
	NP=0
	NG=0
	NN=0
	nchsh=0 														!obschee chislo shagov po vremeni
	HX = XXK/(MK-1) 												!prostranstvennyj shag
	dtm=hx/cs														!"magnitnyj" shag po vremeni
	DO M = 1,MK
		X(M) = HX*(M-1) 											 !Z-koordinata uzlov setki
	END DO
#ifdef USE_AV
	call faglStartWatch(Vs, status)
	call faglStartWatch(X, status)
	call favStartViewer(hv, status)
	call favSetArray(hv, Vs, status)
	call favSetDimScale(hv, 1, X, status)
	call favSetDimName (hv, 1, "X", status)
	call favShowWindow(hv, AV_TRUE, status)
#endif
C --------------------------------------------------------
C	Make up one time step
C --------------------------------------------------------
	write(*,*) 'Emin=',Emin
	write(*,*) 'Rmin=',Rmin
	write(*,*) 'Cmin=',Cmin
	write(*,*) 'Pmin=',Pmin


3		N=N+1
	NP=NP+1
	NG=NG+1
	ALMAX=0.
c ----------- boundary impulse -----------------------------
c	if (T .lt. 1d0/chi) then
c		VF(1,5) = cep*dsin(6.2831858*chi*T)
c		VF(1,6) = cep/csv*dsin(6.2831858*chi*T)
c	end if

c	if (T .lt. 1d0/chi) then
c		Wfi1 = 2*cep*dsin(6.2831858*chi*T)
c		Wxi1 = VF(MK-1,5) - csv*VF(MK-1,6)
c		VF(MK-1,5) = (Wfi1+Wxi1)/2d0
c		VF(MK-1,6) = -(Wfi1+Wxi1)/2d0/csv
c	end if

c --- One Impulse -------
c	if (T .lt. 1d0/chi) then
c		Wfi1 = 2*cep*dsin(6.2831858*chi*T)
c		Wxi1 = VF(1,5) - csv* VF(1,6)
c		VF(1,5) = (Wfi1+Wxi1)/2d0
c		VF(1,6) = (Wfi1+Wxi1)/2d0/csv
c	end if
c --- Infinite Series of Impulse Packs -------
	  apacknum = T * csi
	  mpacknum = apacknum
	  Tps = mpacknum/csi
	  if (T-Tps .lt. cpi/chi) then
		 Wfi1 = 2*cep*dsin(6.2831858*chi*(T-Tps))
		 Wxi1 = VF(1,5) - csv* VF(1,6)
		 VF(1,5) = (Wfi1+Wxi1)/2d0
		 VF(1,6) = (Wfi1+Wxi1)/2d0/csv
	  end if
	  

C -------------------------------------------------------
	V(:,1)=VF(:,1)*cme				  !!!!V(:,1)=VF(:,1)
	V(:,2)=VF(:,3) !+cs !-cs
	V(:,3)=(VF(:,4)*cze)/cme    !!!!!	V(:,3)=VF(:,4)
	V(:,4)=VF(:,2)
C -------------------------------------------------------    
C	novaya pechat (v otnositelnykh - meghdu min i max - velichinakh)
	do i=1,ik !1,ik
		amn(i)=1.e20
		amx(i)=-1.e20
		do m=1,mk
			if(vf(m,i).lt.amn(i))then
				amn(i)=vf(m,i)
			end if
			if(vf(m,i).gt.amx(i))then
				amx(i)=vf(m,i)
			end if
		end do
	end do
	do m=1,mk
c		do i=1,ik
c			vs(m,i)=(vf(m,i)-amn(i))/(amx(i)-amn(i)+1.e-6)*(i+3)
c		end do

		vs(m,1)=(vf(m,1))/(amx(1)+1.e-6)*(1+3)
		vs(m,2)=(vf(m,2)-amn(2))/(amx(2)-amn(2)+1.e-6)*(2+3)
		vs(m,3)=(vf(m,3)-amn(3))/(amx(3)-amn(3)+1.e-6)*(3+3)
		vs(m,4)=(vf(m,4))/(amx(4)+1.e-6)*(4+3)
		vs(m,5)=(vf(m,5)-amn(5))/(amx(5)-amn(5)+1.e-6)*(5+3)
		vs(m,6)=(vf(m,6)-amn(6))/(amx(6)-amn(6)+1.e-6)*(6+3)
		vs(m,7)=(vf(m,7)-amn(7))/(amx(7)-amn(7)+1.e-6)*(7+3)

c		vs(m,1)=(vf(m,1))/(amx(1)+1.e-6)+(0)
c		vs(m,2)=(vf(m,2)-amn(2))/(amx(2)-amn(2)+1.e-6)+(1)
c		vs(m,3)=(vf(m,3)-amn(3))/(amx(3)-amn(3)+1.e-6)+(2)
c		vs(m,4)=(vf(m,4))/(amx(4)+1.e-6)+(3)
c		vs(m,5)=(vf(m,5)-amn(5))/(amx(5)-amn(5)+1.e-6)+(4)
c		vs(m,6)=(vf(m,6)-amn(6))/(amx(6)-amn(6)+1.e-6)+(5)
c		vs(m,7)=(vf(m,7)-amn(7))/(amx(7)-amn(7)+1.e-6)+(6)

	end do
C --------------------------------------------------------
C	Output (zapominanie) gazodinamicheskikh results 
C --------------------------------------------------------
	if (.true.)then
	IF(NP >= NPK .or. T >= TK-DT)THEN
		NP=0
		PRINTMORE=.TRUE.
		if (.true.) then
			OPEN(UNIT=2,FILE='Test1D.out',STATUS='REPLACE',iostat=istat)
			IF(istat /=0) STOP "***Cannot open Test1D.out***"
			REWIND 2
			WRITE(2,*)N,T
			WRITE(2,*) X
			WRITE(2,*) V
			WRITE(2,*)
			CLOSE(UNIT=2)
			OPEN(UNIT=5,FILE='StateDump',STATUS='REPLACE',iostat=istat)
			IF(istat /=0) STOP "***Cannot open StateDump***"
			REWIND 5
			WRITE(5,*)N,T,ZT
			WRITE(5,*) X
			WRITE(5,*) VF
			WRITE(5,*)
			CLOSE(UNIT=5)
		end if

		WRITE(77,*)
		WRITE(77,*)N,T,ZT
		do m=1,mk
			WRITE(77,*)VF(m,1),VF(m,2),VF(m,3),VF(m,4),VF(m,5),VF(m,6),VF(m,7)
		end do
4	END IF
	end if
C --------------------------------------------------------
C	grafika VF(TN,Z):sinij-Ne, zelenyj-Vx , krasnyj-Vz, goluboj-Ee ,fioletovyj-Ex, gheltyj- Hy, korichnevyj-Ez
C --------------------------------------------------------
	IF (NG >= NPG) THEN
#ifdef USE_AV
		call favUpdate(hv, 0, status)
#endif
		NG=0
	END IF
C --------------------------------------------------------
C	Pre-step checks
C --------------------------------------------------------
	DO M=1,MK
		IF (VF(M,1) .LT. 0) THEN
			WRITE(*,*) 'Negative ne (m=',m,')' 
		END IF
	END DO
C --------------------------------------------------------
C	Gas-Dynamics-Step:
C	Evaluate time step
C --------------------------------------------------------
	DO 9 M = 1, MK
		R=V(M,1)
		U=V(M,2)
		E=V(M,3)
		CALL URSOST(R,E,G,C,P,PR,RE)
		ALM(M) = DABS(U)+C
		IF(ALM(M) > ALMAX) ALMAX = ALM(M)
9		CONTINUE
	IF(ALMAX <= 1.E-8) ALMAX=1.E-8
	DTG = HX *AKYR /ALMAX !gazodinamicheskij shag po vremeni
C-----------------------------------------------------------------------------------------------------------------
C	Vybor shaga iz dtg i dtpr i popravka do celogo chisla dtm
	if(dtg.gt.dtpr)then
		dt=dtpr
	else
		dt=dtg
	end if
c		if(t.lt.dt)then
c			dt=1.e-19
c		end if
	amag=dt/dtm
	mmag=amag				  !cislo "magnitnykh" shagov v gazodinamicheskom
	nchsh=nchsh+1			  !obschee chislo shagov po vremeni
	if(mmag.ge.1)then
		dt=mmag*dtm 			  !Popravka gazodinamicheskogo shaga po vremeni (do celogo chisla "magnitnykh" shagov)
	end if
	SIG=dt*AKYR*ALMAX/hx
	B0=B00*DT																 ! const skhemy pravoj chasti
	B1=B10*DT																 ! const skhemy pravoj chasti
	C0=C00*DT*DT															 ! const skhemy pravoj chasti
	C1=C10*DT*DT															 ! const skhemy pravoj chasti
C-----------------------------------------------------------------------------------------------------------------
	T = T + DT																 !tn1=tn+dt
C --------------------------------------------------------
	zt=z0+cs*t																!tekyschaya Z-koordinata
C --------------------------------------------------------
	IF(NN.EQ.NG)THEN
		NN=0
		WRITE(*,*)'N = ',N,' T = ',T,' dT = ',DT
		WRITE(*,*)'DTG= ',DTG,' DTPR= ',dtpr,' DTM= ',DTM
		WRITE(*,*)'Z= ',zt,'mmag= ',mmag
		WRITE(*,*)'zn= ',zn,'zne= ',zne,'zni= ',zni,'zn0= ',zn0,'zno2= ',
     &			zno2
		WRITE(*,*)'AMN = ',AMN													!pechat min znachenij parametrov
		WRITE(*,*)'AMX = ',AMX													!pechat max znachenij parametrov

#ifdef USE_AV
	write(ArrayName,'(A,E10.3,A,E10.3,A,E10.3,A,E10.3,A,E10.3,A,E10.3
     &,A,E10.3,A,E10.3,A,E10.3,A,E10.3,A,E10.3
     &,A,E10.3,A,E10.3,A,E10.3,A,E10.3)')
     &'ne=',AMN(1),':',AMX(1)
     &,',vx=',AMN(2),':',AMX(2)
     &,',vz=',AMN(3),':',AMX(3)
     &,',Et=',AMN(4),':',AMX(4)
     &,',Ex=',AMN(5),':',AMX(5)
     &,',Hy=',AMN(6),':',AMX(6)
     &,',Ez=',AMN(7),':',AMX(7)
#endif
#ifdef STEP_PAUSE
		read(*,*)
#endif
#ifdef USE_AV
		call favSetArrayName(hv, trim(ArrayName),status)
#endif
	END IF
C --------------------------------------------------------
C	Prediction step (Godunov scheme)
C --------------------------------------------------------

#ifdef FPOT_OUT
	WRITE(80,*)
	WRITE(80,*)N,T,ZT
#endif

	V0 = V(1,:)
	W = V(2,:)
	CALL FPOT(N,1,RID, G, V0, IUU, W, UN1U,var)
#ifdef FPOT_OUT
	WRITE(80,*)UN1U(1),UN1U(2),UN1U(3),UN1U(4)
#endif
	Chm(1)=var
C -----------------------------------------------------------
	DO 79 M = 2, MK-1
		V0 = V(M,:)
		W = V(M+1,:)
		CALL FPOT(N,M,RID, G, V0, IUU, W, UN2U,var)

#ifdef FPOT_OUT
		WRITE(80,*)UN2U(1),UN2U(2),UN2U(3),UN2U(4)
#endif

		Chm(M)=var
		CALL UPOT( V0, IUU, UP1U )
		UN(M,1:3) = UP1U
	    
		UP2U = UP1U - dt*(UN2U - UN1U)/hx
	    if(UP2U(1).lt.0.0) then
	     write(*,*) 'V0:',V0
	     write(*,*) 'UP1U:',UP1U
		 write(*,*) 'UN2U:',UN2U
		 write(*,*) 'UN1U:',UN1U
		 write(*,*) 'SIG:',SIG
		 write(*,*)'M:',M
	     stop
	    end if
		UP(M,1:3) = UP2U
		UN1U = UN2U
		IF(BG < 1.5)THEN
			CALL VPOT( UP2U, IUU, FMU )
			V1(M,:) = FMU
		END IF
79	CONTINUE

C -------------------------------------------------------
C	Impose boundary conditions
C -------------------------------------------------------
	Chm(MK)=Chm(MK-1)
	IF(BG <1.5)THEN
		UN1U = V1(2,:)
		UN2U = V1(MK-1,:)
		CALL BOUNDARY(IUU,UN1U,UN2U,UP1U,UP2U,NL,NR)
		V1(1,:) = UP1U
		V1(MK,:) = UP2U
	ELSE
C -------------------------------------------------------
		UN1 = UN(2,:)
		UN2 = UN(MK-1,:)
		CALL BOUNDARY(IU,UN1,UN2,UP1,UP2,NL,NR)
		UN(1,:) = UP1
		UN(MK,:) = UP2
C -------------------------------------------------------
		UN1 = UP(2,:)
		UN2 = UP(MK-1,:)
		CALL BOUNDARY(IU,UN1,UN2,UP1,UP2,NL,NR)
		UP(1,:) = UP1
		UP(MK,:) = UP2
	END IF
C --------------------------------------------------------
C	Correction step (2-3 order approximate scheme GKhM2li)
C --------------------------------------------------------
	IF(BG > 1.5)THEN
C---------------------------------------------------------------------
C	  Calculation of the flux in the point m =1/2
C---------------------------------------------------------------------
		UN1 = UN(1,:)
		UN2 = UN(2,:)
		UP1 = UP(1,:)
		UP2 = UP(2,:)
		B = (UN1+UN2+UP1+UP2)/4.
		CALL VPOT( B, IU, V0 )
		CALL SZN(G,V0,IU,SIG,S)
		CALL OMEGA(G,V0,OM,OMO,IU)	   !Calculate matrixes Omega and inverse to it
		CALL GIBR(0,BG,IU,EPS,UN1,UN1,UN2,UN2,UP1,UP2,OM,S,
     &			GIB,RUS,V0,FM)	  !Calculate parametres of hybridity
		CALL BCD(RUS,S,IU,GIB,B,CC) 	 !Calculate coeficients of diagonal matrixes [B] and [C]
		DO I=1,IU
			UN1(I) = B(I)*V0(I)
			UP1(I) = CC(I)*FM(I)
		END DO
		W = MATMUL(OMO,(UN1 + UP1))
C---------------------------------------------------------
C	  Start cycle by x
C---------------------------------------------------------
		DO 99 M =2,MK-1
			UB	= UN(M-1,:)
			UN1 = UN(M,:)
			UN2 = UN(M+1,:)
			FF =  UP(M-1,:)
			UP1 = UP(M,:)
			UP2 = UP(M+1,:)
			IF(M < MK-1)THEN
				FM = UN(M+2,:)
				CC = UP(M+2,:)
				MM = 1
			ELSE
				FM = UN(M+1,:)
				CC = UP(M+1,:)
				MM = 0
			END IF
			B=(UB+UN1+UN2+FM+FF+UP1+UP2+CC)/8.
			CALL VPOT( B, IU, V0 )
			CALL SZN(G,V0,IU,SIG,S)
			CALL OMEGA(G,V0,OM,OMO,IU)		 !Calculate matrixes Omega and inverse to it
			CALL GIBR(MM,BG,IU,EPS,UB,UN1,UN2,FM,UP1,UP2,OM,S,
     &			   GIB,RUS,V0,FF)	 !Calculate parametres of hybridity
			CALL BCD(RUS,S,IU,GIB,B,CC) 	 !Calculate coeficients of diagonal matrixes [B] and [C]
			DO I=1,IU
				UN1(I) = B(I)*V0(I)
				UP1(I) = CC(I)*FF(I)
			END DO
			W1 = MATMUL(OMO,(UN1 + UP1))
C--------------------------------------------------------------------
C	  Calculation of the final values of U on the layer t=tn+dt
C---------------------------------------------------------------------
			IF(Chm(M-1) == "vacum".or. Chm(M) == "vacum"
     &					.or.Chm(M+1) == "vacum")THEN
				FM = UP(M,:)
			ELSE
				FM = UP(M,:) + W1 - W
			END IF
			CALL VPOT(FM,IU,V0)
			V1(M,:) = V0
			W = W1
99			CONTINUE
C---------------------------------------------------------
C	  End cycle by z
C--------------------------------------------------------------------
C	  Impose boundary conditions
C--------------------------------------------------------------------
		UN1 = V1(2,:)
		UN2 = V1(MK-1,:)
		CALL BOUNDARY(IU,UN1,UN2,UP1,UP2,NL,NR)
		V1(1,:) = UP1
		V1(MK,:) = UP2
	END IF
	if (.true.) then ! hypoff
C ----------------------------------------------------------
	FV(:,1)=VF(:,1) ! Remember the values on the N layer
	VF(:,1)=V1(:,1)/cme    !!!!
C ----------------------------------------------------------
	FV(:,3)=VF(:,3) ! Remember the values on the N layer
	VF(:,3)=V1(:,2) !+cs
C ----------------------------------------------------------
	FV(:,4)=VF(:,4) ! Remember the values on the N layer
C	VF(:,4)=V1(:,3)
	VF(:,4)=(V1(:,3)/cze)*cme
	FV(:,2)=VF(:,2) ! Remember the values on the N layer
	VF(:,2)=V1(:,4)    !!!!
C	write (666,*) 'Et(hyp1)=', VF(1,4) - FV(1,4)
C--------------------------------------------------------------------
C	  Impose additional boundary conditions
C --------------------------------------------------------
	UN1(1) = VF(1,2)
	UN1(2) = VF(1,5)-VcIm*VF(1,6)
	UN1(3) = VF(1,7)
	
	UN2(1) = VF(MK,2)
	UN2(2) = VF(MK,5)-VcIm*VF(MK,6)
	UN2(3) = VF(MK,7)
C --------------------------------------------------------
	CALL BOUNDARY(IU,UN1,UN2,UP1,UP2,NL,NR)
C ----------------------------------------------------------
C	Culculation the Transfer Equations:
C ----------------------------------------------------------
C	FV(:,2)=VF(:,2) ! Remember the values on the N layer
C ----------------------------------------------------------
	FV(:,7)=VF(:,7) ! Remember the values on the N layer
C -----------------------------------------------------------
C	Invariant Psi = Ex-Ro*Hy on the N+1 layer:
C ----------------------------------------------------------
	VU = UP2(2)
	vl=-cs
	Cr=DABS(Vl*DT/HX)
	VM1(:)=FV(:,5)-VcIm*FV(:,6)
	VN(:)=VF(:,5)-VcIm*VF(:,6)
C -----------------------------------------------------------
	IF( mmag.ne.0)THEN
		DO M=MK,1,-1
			NB=(mmag+m-mk)/(mk-2)
			if(mmag.le.mk-m)then
				jmag=mmag+m
			else
				jmag=mmag-nb*(mk-2)-(mk-m)+2
			end if
c				write(*,*)'t=',t,'m=',m,'mma=',mmag,'mk=',mk,'nb=',nb,'jma=',jmag
			vn1(m)=vn(jmag) 					   !snoc pri Ejlerov. s/k
		end do
	else
		CALL SKHEMA(N,VU,VM1,VN,VN1,MK,Cr,Vl) !t)
	end if
	vm2(:)=vn1(:)						   !dobavlen operator pri Ejlerovoj s/k
C --------------------------------------------------------
C	  Impose additional boundary conditions
C --------------------------------------------------------
	UN1(1) = VF(1,2)
	UN1(2) = VF(1,5)+VcIm*VF(1,6)
	UN1(3) = VF(1,7)

	UN2(1) = VF(MK,2)
	UN2(2) = VF(MK,5)+VcIm*VF(MK,6)
	UN2(3) = VF(MK,7)
	CALL BOUNDARY(IU,UN1,UN2,UP1,UP2,NL,NR)
C ----------------------------------------------------------
C ----------------------------------------------------------
C	Invariant Fi = Ex+Ro*Hy on the N+1 layer:
C ----------------------------------------------------------
	VU =UP1(2)								   !dobavlen operator pri Ejlerovoj s/k
	vl=cs
	Cr=DABS(Vl*DT/HX)
	VM1(:)=FV(:,5)+VcIm*FV(:,6) 		   !dobavlen operator pri Ejlerovoj s/k
	VN(:)=VF(:,5)+VcIm*VF(:,6)
	IF( mmag.ne. 0)THEN !rascet pri DTgaz ili DTpravchast < DTmagnitnom, inache - snos
		DO M=1,mk !2,MK
			NB=(mmag-m)/(mk-2)
			if(mmag-m.le.-1)then
				jmag=m-mmag
			else
				jmag=m-mmag+(nb+1)*(mk-2)
			end if
c				write(*,*)'t=',t,'m=',m,'mma=',mmag,'mk=',mk,'nb=',nb,'jma=',jmag
			vn1(m)=vn(jmag) !m+1)	  !snoc pri Ejlerov. s/k
		end do
	else
		CALL SKHEMA(N,VU,VM1,VN,VN1,MK,Cr,Vl) !t) !dobavlen operator pri Ejlerovoj s/k
	end if
C ----------------------------------------------------------
	FV(:,5)=VF(:,5) ! Remember the values on the N layer
	FV(:,6)=VF(:,6) ! Remember the values on the N layer
C -----------------------------------------------------------
	!izmenen operator pri Ejlerovoj s/k(VN(:)+VN1(:))/2. !Electric field in X direction on the N+1 layer
	VF(:,5)=(vn1(:)+vm2(:))/2.
	!izmenen operator pri Ejlerovoj s/k(VN(:)-VN1(:))/VcIm/2.  !Magnetic field in Y direction on the N+1 layer		
	VF(:,6)=(vn1(:)-vm2(:))/VcIm/2.
	end if
C -----------------------------------------------------------
C	End calculation of the correction, nachalo ucheta pravoj chasti
C -----------------------------------------------------------
#ifdef AFTER_OUT
	if (PRINTMORE) then
		WRITE(78,*)
		WRITE(78,*)N,T,ZT
		do m=1,mk
			WRITE(78,*)VF(m,1),VF(m,2),VF(m,3),VF(m,4),VF(m,5),VF(m,6),VF(m,7)
		end do
	end if
#endif

c		WRITE(555,*)
c		WRITE(555,*)N,T,ZT

#ifndef RHS_DISABLE
	do m=1,mk
		yn(1)=VF(m,1)
		yn(2)=VF(m,3)
		yn(3)=VF(m,2)
		yn(4)=VF(m,4)
		yn(5)=VF(m,5)+VcIm*VF(m,6)
		yn(6)=VF(m,5)-VcIm*VF(m,6)
		yn(7)=VF(m,7)
		call ATMOSPHERE(t,zt,yn(1))
		FPtmp(:)=FP(m,:)
		call DISSIPATIONSTEP(dtpr,yn,y,FPtmp,m,N)
		FP(m,:)=FPtmp(:)
		VF(m,1)=y(1)														   !okonchatelnyj vector VF(tn+dt,z) s pr. chastyu F(y)
		VF(m,2)=y(3)
		VF(m,3)=y(2)
		VF(m,4)=y(4)
		VF(m,5)=(y(5)+y(6))/2d0
		VF(m,6)=(y(5)-y(6))/VcIm/2d0
		VF(m,7)=y(7)
	end do
#endif

	if (PRINTMORE) then
		PRINTMORE=.FALSE.
#ifdef AFTER_OUT
		WRITE(79,*)
		WRITE(79,*)N,T,ZT
		do m=1,mk
			WRITE(79,*)VF(m,1),VF(m,2),VF(m,3),VF(m,4),VF(m,5),VF(m,6),VF(m,7)
		end do
#endif
	end if

321	continue
cc	dtpr=dtpr*1. !1.01 !1.005 !  novyj shag dlya pravykh chastej													   !Konec vychisleniya pravykh chastej
C -----------------------------------------------------------------------
	IF(key1 == "y")THEN !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !Rasschirenie oblasti ?
C ----------------------------------------------------------------------
C	Evaluate new region dimensions
C ----------------------------------------------------------------------
		DL=0.; DR=0.
		DO I=1,IK
			Dtl=2.*DABS(VF(10,I)-FV(10,I))/(VF(10,I)+FV(10,I))
			Dtr=2.*DABS(VF(MK-10,I)-FV(MK-10,I))/(VF(MK-10,I)+FV(MK-10,I))
			IF(Dtl >= DL) DL=Dtl
			IF(Dtr >= DR) DR=Dtr
		END DO
		WRITE(*,*)'Dl= ',DL,' Dr= ',DR
C ----------------------------------------------------------------------
		IF(DL > DELTA .or. DR > DELTA)THEN
			XXK=XXK*2.
			HX=XXK/(MK-1) ! New step in X
			DO M=1,MK
				X(M)=(M-1)*HX  ! X-koordinata uzlov
			END DO
C ----------------------------------------------------------------------
#ifdef USE_AV
			call faglUpdate (X, status)
#endif
C ----------------------------------------------------------------------
			DO I=1,IK
				VM1(:)=FV(:,I)
				VN(:)=VF(:,I)
				DO M=2,MK-1
					IF(M <= MK/4)THEN
						FV(M,I)=FV(1,I)
						VF(M,I)=VF(1,I)
					ELSE
						IF(M <= 3*MK/4)THEN
							FV(M,I)=VM1(2*(M-MK/4))
							VF(M,I)=VN(2*(M-MK/4))
						ELSE
							FV(M,I)=FV(MK,I)
							VF(M,I)=VF(MK,I)
						END IF
					END IF
				END DO
			END DO
		END IF
	END IF
C ------------------------------------------------------------
C	New time step
C ------------------------------------------------------------
	IF( T < TK )THEN
		GO TO 3
	ENDIF
C -----------------------------------------------------------
#ifdef USE_AV
	call faglClose(VF, status)
	call faglClose(X, status)
	call favEndViewer(hv, status)
#endif
C -----------------------------------------------------------
	CLOSE(77)
	STOP
	END PROGRAM MAIN



C ----------------------------------------------------------
C	Calculation boundary conditions
C ----------------------------------------------------------
	SUBROUTINE BOUNDARY(IK,V2,VM1,V1,VM,NL,NR)
	IMPLICIT REAL*8 (a-h,o-z)
	DIMENSION:: V2(IK),VM1(IK),V1(IK),VM(IK)
	INTENT(IN):: IK,V2,VM1,NL,NR
	INTENT(OUT):: V1,VM
	V1 = V2 						   !transparent BC
	VM = VM1						   !transparent BC
	IF(NL == 0) V1(2)= -V2(2)		   !reflective BC
	IF(NR == 0) VM(2)= -VM1(2)		   !reflective BC
	IF(NL == 2) V1 = VM1			   !periodic BC
	IF(NR == 2) VM = V2 			   !periodic BC
	RETURN
	END SUBROUTINE BOUNDARY
C ----------------------------------------------------------
C	Calculation fuction "F(V)"
C ----------------------------------------------------------
	SUBROUTINE FPOT(N,M,RID, G, V, IU, W, F, var)
	IMPLICIT REAL*8 (a-h,o-z)
	DIMENSION V(IU), W(IU), F(IU), VP(IU)
	INTENT(IN):: N, M, RID, G, V, IU, W
	INTENT(OUT):: F
	CHARACTER(5) :: var
	COMMON /MIN/ Rmin,Emin,Cmin,Pmin
	R1=V(1)
	U1=V(2)
	E1=V(3)
	R2=W(1)
	U2=W(2)
	E2=W(3)
	IF(E1 <= Emin) E1=Emin 
	IF(R1 <= Rmin) R1=Rmin 
	IF(E2 <= Emin) E2=Emin
	IF(R2 <= Rmin) R2=Rmin 
	G1=G
	G2=G
C	CALL RIEMANN(N,M,RID,R1,E1,U1,G1,R2,E2,U2,G2,P,R,E,U,var)
	CALL RIEMANN(N,M,RID,R1,E1,U1,G1,R2,E2,U2,G2,P,R,E,U)
	F(1) = R*U
	F(2) = F(1)*U+P
	F(3) = (R*(E+U*U/2.)+P)*U
	DO I = 4,IU
		IF (U .GT. 0D0) THEN
			F(I) = R*U*V(I)
		ELSE
			F(I) = R*U*W(I)
		END IF
	END DO
	RETURN
	END SUBROUTINE FPOT
C --------------------------------------------------------
C	  Calculation fuction "U(V)"
C --------------------------------------------------------
	SUBROUTINE UPOT( V, IU, U )
	IMPLICIT REAL*8 (a-h,o-z)
	DIMENSION V( IU ), U( IU )
	INTENT(IN):: V, IU
	INTENT(OUT):: U
	  U(1)=V(1)
	  U(2)=V(1)*V(2)
	  SQRD = V(2)*V(2)
	  DO I = 4,IU
	    U(I)=V(1)*V(I)
	    SQRD=SQRD+V(I)*V(I)
	  END DO
	  U(3)=V(1)*(V(3)+SQRD/2d0)
	RETURN
	END SUBROUTINE UPOT
C --------------------------------------------------------
C	  Calculation fuction "V(U)"
C --------------------------------------------------------
	SUBROUTINE VPOT( U, IU, V )
	IMPLICIT REAL*8 (a-h,o-z)
	DIMENSION V( IU ), U( IU )
	INTENT(IN):: U, IU
	INTENT(OUT):: V
	COMMON /MIN/ Rmin,Emin,Cmin,Pmin
	IF(U(1) <= Rmin)THEN
        WRITE(*,*)'-----------',U,Rmin
	  stop
	  V(1)=Rmin; V(2)=0; V(3)=Emin
	ELSE
	  V(1)=U(1)
	  V(2)=U(2)/V(1)
	  SQRD = V(2)**2
	  DO I = 4,IU
	    V(I)=U(I)/V(1)
	    SQRD=SQRD+V(I)**2
	  END DO
	  V(3)=U(3)/V(1)-SQRD/2.
	END IF
	IF(V(3) < Emin) V(3)=Emin
	RETURN
	END SUBROUTINE VPOT
C --------------------------------------------------------------
C	  Calculation new values of  P,R,E,U  by Godunov Method
C---------------------------------------------------------------
      SUBROUTINE RIEMANN(N,M,RID,R1in,E1,U1,GN1,R2in,E2,U2,GN2,P,R,E,U)
      IMPLICIT REAL*8 (a-h,o-z)
      INTENT(IN):: N,M,RID,R1in,E1,U1,GN1,R2in,E2,U2,GN2
      INTENT(OUT):: P,R,E,U
      CHARACTER(5) :: var
      COMMON /MIN/ Rmin,Emin,Cmin,Pmin
            logical vak_case
C --------------------------------------------------------------
      Amin=Rmin*Cmin
      CALL URSOST(R1in,E1,GN1,C1,P1,PR,PE)
      CALL URSOST(R2in,E2,GN2,C2,P2,PR,PE)
     
!        Uvak=dabs(2.0*(sqrt(C1)+sqrt(C2)))/(GN1-1.0)
        Uvak=dabs(2.0*(C1+C2))/(GN1-1.0)
        DU=U2-U1
!        if(DU.gt.Uvak) then
!                vak_case=.TRUE.
!!                write(*,*)'Vacuum! Uvak=',Uvak,'DU=',DU
!        else
                vak_case=.FALSE.
!        end if
       C1=max(C1,Cmin)
       C2=max(C2,Cmin)
      		 P1=max(P1,Pmin)
      		 P2=max(P2,Pmin)
      		 R1=max(R1in,Rmin)
      		 R2=max(R2in,Rmin)
C --------------------------------------------------------------

C	Simplest Variant

C --------------------------------------------------------------
      RC=R1*C1+R2*C2
      IF(RC < Amin) RC=Amin
      Pa=(P2*R1*C1+P1*R2*C2+(U1-U2)*R1*C1*R2*C2)/RC
      IF(Pa < Pmin) Pa=Pmin
      Ua=(U1*R1*C1+U2*R2*C2+P1-P2)/RC
      D12=DSQRT((U1-U2)*(U1-U2)+(P1-P2)*(P1-P2))
      IF(D12 < 1.E-10)THEN
        var = "simpl"
        D1=U1-C1
        D2=U2+C2
        IF(Ua >= 0.)THEN
          IF(D1 > 0.)THEN
            P=P1
            U=U1
      	  R=R1
            E=E1
          ELSE
            P=Pa
            U=Ua
            R=GN1*Pa/(C1*C1)
      	  CALL VNENER(R,P,GN1,E)
      	END IF
        ELSE
          IF(D2 < 0.)THEN
            P=P2
            U=U2
      	  R=R2
            E=E2
          ELSE
            P=Pa
            U=Ua
            R=GN2*Pa/(C2*C2)
      	  CALL VNENER(R,P,GN1,E)
          END IF
        END IF
       GO TO 11
      END IF
C---------------------------------------------------------------
C     Check the Variants Dp-Dp, Sh-Dp, Dp-Sh or Sh-Sh:
C---------------------------------------------------------------
      Gd1=(GN1-1.)/2./GN1
      Gd2=(GN2-1.)/2./GN2
      DU = U1-U2
      Uvc = -2.*(C1/(GN1-1.)+C2/(GN2-1.))
      IF(DU < Uvc)THEN  ! We can not get solution!!
         var = "vacum"
         GO TO 10
      END IF
      A1=DSQRT(DABS(R1*(P2*(GN1+1.)/2.+P1*(GN1-1.)/2.)))
      A2=DSQRT(DABS(R2*(P1*(GN2+1.)/2.+P2*(GN2-1.)/2.)))
      IF(A1 < Amin) A1=Amin
      IF(A2 < Amin) A2=Amin
      Ush1 = (P2-P1)/A1
      Ush2 = (P1-P2)/A2
  
      Udp1 = -2.*C2*(1.-DABS(P1/P2)**Gd2)/(GN2-1.)
      Udp2 = -2.*C1*(1.-DABS(P2/P1)**Gd1)/(GN1-1.)
C------------------------------------------------------------------------
C     Logical Block: dp - depresision, tn - tangential, sh - shock waves.
C------------------------------------------------------------------------
      IF(P2 >= P1)THEN
         IF(DU >= Uvc .AND. DU < Udp1) var = "dp-dp"
         IF(DU == Udp1) var = "tn-dp"
         IF(DU > Udp1 .AND. DU < Ush1) var = "sh-dp"
         IF(DU == Ush1) var = "sh-tn"
         IF(DU > Ush1) var = "sh-sh"
      ELSE
         IF(DU >= Uvc .AND. DU < Udp2) var = "dp-dp"
         IF(DU == Udp2) var = "dp-tn"
         IF(DU > Udp2 .AND. DU < Ush2) var = "dp-sh"
         IF(DU == Ush2) var = "tn-sh"
         IF(DU > Ush2) var = "sh-sh"
      END IF
      IF(var == "dp-dp".or.var == "tn-dp".or.var == "dp-tn".or.
     &   var == "tn-sh".or.var == "sh-tn") GO TO 10
C-------------------------------------------------------------------------
C     Noniterative Case: we calculate P and U in case of weak discontinuity  
C-------------------------------------------------------------------------
      SELECT CASE(var)
C-----Depresision - Shock Waves-------------------------------------------
      CASE("dp-sh")
      G1=GN1
      G2=(GN2+1.)/2.
C-----Shock - Depresision Waves-------------------------------------------
      CASE("sh-dp")
      G1=(GN1+1.)/2.
      G2=GN2
C-----Shock - Shock Waves-------------------------------------------------
      CASE DEFAULT
      G1=(GN1+1.)/2.
      G2=(GN2+1.)/2.
      END SELECT
C-------------------------------------------------------------------------
      RG=DABS(R1*G1*R2*G2)**(1.5d0)
      DN=2.*(DSQRT(DABS(R1*G1))*R2*G2+DSQRT(DABS(R2*G2))*R1*G1)
      IF(DN < Amin) DN=Amin
      Pe=((R1*R2*(U1-U2)*G1*G2+
     &DSQRT(DABS(4.*(P1+P2)*RG+4.*R1*R2*G1*G2*(P2*R1*G1+P1*R2*G2)+
     &        (R1*R2*U1*G1*G2-R1*R2*U2*G1*G2)**2)))/DN)**2

      IF(Pe < Pmin) Pe=Pmin
      IF(DABS((Pe-Pa)/Pa) <= RID)THEN
        P=Pe
        SELECT CASE(var)
C-------Depresision - Shock Waves-------------------------------------------
        CASE("dp-sh")
        Pd1 = 1.-DABS(Pe/P1)**Gd1
        IF(DABS(Pd1) < 1.E-10) Pd1=1.E-10
        A1=Gd1*R1*C1*(1.-Pe/P1)/Pd1
        A2=R2*C2*DSQRT(DABS((GN2+1.)*Pe/P2/GN2/2.+Gd2))
C-------Shock - Depresision Waves-------------------------------------------
        CASE("sh-dp")
        A1=R1*C1*DSQRT(DABS((GN1+1.)*Pe/P1/GN1/2.+Gd1))
        Pd2 = 1.-DABS(Pe/P2)**Gd2
        IF(DABS(Pd2) < 1.E-10) Pd2=1.E-10
        A2=Gd2*R2*C2*(1.-Pe/P2)/Pd2
C-------Shock - Shock Waves-------------------------------------------------
        CASE DEFAULT
        A1=R1*C1*DSQRT(DABS((GN1+1.)*Pe/P1/GN1/2.+Gd1))
        A2=R2*C2*DSQRT(DABS((GN2+1.)*Pe/P2/GN2/2.+Gd2))
        END SELECT
C-------------------------------------------------------------------------
        IF(A1 < Amin) A1=Amin
        IF(A2 < Amin) A2=Amin
        U=(U1-(Pe-P1)/A1+U2+(Pe-P2)/A2)/2.
        GO TO 10
      END IF
C-------------------------------------------------------------------------
C     Iterative Block: Newton Iterations
C-------------------------------------------------------------------------
      PI=Pa
      I=0
9	CONTINUE
      IF(PI < Pmin) PI=Pmin
      IF(PI >= P1)THEN
        A1=R1*C1*DSQRT(DABS((GN1+1.)*PI/P1/GN1/2.+Gd1))
        IF(A1 < Amin) A1=Amin
        F1=(PI-P1)/A1
      ELSE
        F1=2.*C1*(DABS(PI/P1)**Gd1-1.)/(GN1-1.)
      END IF
      IF(PI >= P2)THEN
        A2=R2*C2*DSQRT(DABS((GN2+1.)*PI/P2/GN2/2.+Gd2))
        IF(A2 < Amin) A2=Amin
        F2=(PI-P2)/A2
      ELSE
        F2=2.*C2*(DABS(PI/P2)**Gd2-1.)/(GN2-1.)
      END IF
      IF(PI >= P1)THEN
        AA1=4.*GN1*R1*C1*DABS((GN1+1.)*PI/P1/GN1/2.+Gd1)**1.5d0
        IF(AA1 < Amin) AA1=Amin
        DF1=((GN1+1.)*PI/P1+3.*GN1-1.)/AA1
      ELSE
        DF1=C1*((PI/P1)**Gd1)/PI/GN1
      END IF
      IF(PI >= P2)THEN
        AA2=4.*GN2*R2*C2*DABS((GN2+1.)*PI/P2/GN2/2.+Gd2)**1.5d0
        IF(AA2 < Amin) AA2=Amin
        DF2=((GN2+1.)*PI/P2+3.*GN2-1.)/AA2
      ELSE
        DF2=C2*(DABS(PI/P2)**Gd2)/PI/GN2
      END IF
      DF=DF1+DF2
      IF(DF < 1.E-10) DF=1.E-10
      PI1=PI-(F1+F2-U1+U2)/DF
      I=I+1
      IF(DABS(PI-PI1) > 1.E-7 .AND. I < 40)THEN
        PI=PI1
        GO TO 9
      END IF
      P=PI1
      U=(U1+U2+F2-F1)/2.0d0
C-------------------------------------------------------------------------
C     Now we calculate E and R with known P and U
C-------------------------------------------------------------------------
10	SELECT CASE(var)
C-----Depresision - Depresision Waves-------------------------------------
      CASE("dp-dp")       !  This is true only, if GN1=GN2!!!
      Udv = Udp1-Uvc
      IF(DABS(Udv) < 1.E-10) Udv=1.E-10
      Pt=P1*DABS((U1-U2-Uvc)/Udv)**(1./Gd1)
      Ut=U1-2.*C1*(DABS(Pt/P1)**Gd1-1.)/(GN1-1.)
      D1=U1-C1
      C11=C1+(GN1-1.)*(U1-Ut)/2.
      IF(C11 < Cmin) C11=Cmin
      D11=Ut-C11
      D2=U2+C2
      C22=C2-(GN2-1.)*(U2-Ut)/2.
      IF(C22 < Cmin) C22=Cmin
      D22=Ut+C22
      IF(Ut >= 0.)THEN
        IF(D1 >= 0.)THEN
          P=P1
          U=U1
          R=R1
          E=E1
        ELSE
          IF(D11 > 0.)THEN
            CU1=(2.*C1+(GN1-1.)*U1)/(GN1+1.)
            IF(CU1 < Cmin) CU1=Cmin
      	  U=CU1
      	  R=R1*(CU1/C1)**(2./(GN1-1.))
            P=P1*(CU1/C1)**(2.*GN1/(GN1-1.))
            CALL VNENER(R,P,GN1,E)
      	ELSE
            P=Pt
      	  U=Ut
      	  R=GN1*Pt/(C11*C11)
            CALL VNENER(R,P,GN1,E)
          END IF
        END IF
      ELSE
        IF(D2 <= 0.)THEN
          P=P2
          U=U2
          R=R2
          E=E2
        ELSE
          IF(D22 < 0.)THEN
            CU2=(2.*C2-(GN2-1.)*U2)/(GN2+1.)
            IF(CU2 < Cmin) CU2=Cmin
      	  U=-CU2
      	  R=R2*(CU2/C2)**(2./(GN2-1.))
            P=P2*(CU2/C2)**(2.*GN2/(GN2-1.))
            CALL VNENER(R,P,GN2,E)
      	ELSE
            P=Pt
      	  U=Ut
      	  R=GN2*Pt/(C22*C22)
            CALL VNENER(R,P,GN2,E)
          END IF
        END IF
      END IF
C-----Tangential - Depresision Waves--------------------------------------
      CASE("tn-dp")
      D2=U2+C2
      C22=C2-(GN2-1.)*(U2-U1)/2.
      IF(C22 < Cmin) C22=Cmin
      D22=U1+C22
      IF(U1 >= 0.)THEN
        P=P1
        U=U1
        R=R1
        E=E1
      ELSE
        IF(D2 <= 0.)THEN
          P=P2
          U=U2
          R=R2
          E=E2
        ELSE
          IF(D22 < 0.)THEN
            CU2=(2.*C2-(GN2-1.)*U2)/(GN2+1.)
            IF(CU2 < Cmin) CU2=Cmin
      	  U=-CU2
      	  R=R2*(CU2/C2)**(2./(GN2-1.))
            P=P2*(CU2/C2)**(2.*GN2/(GN2-1.))
            CALL VNENER(R,P,GN2,E)
      	ELSE
            P=P1
      	  U=U1
      	  R=GN2*P/(C22*C22)
            CALL VNENER(R,P,GN2,E)
          END IF
        END IF
      END IF
C-----Depresision - Tangential Waves--------------------------------------
      CASE("dp-tn")
      D1=U1-C1
      C11=C1+(GN1-1.)*(U1-U2)/2.
      IF(C11 < Cmin) C11=Cmin
      D11=U2-C11
      IF(U2 >= 0.)THEN
        IF(D1 >= 0.)THEN
          P=P1
          U=U1
          R=R1
          E=E1
        ELSE
          IF(D11 > 0.)THEN
            CU1=(2.*C1+(GN1-1.)*U1)/(GN1+1.)
            IF(CU1 < Cmin) CU1=Cmin
      	  U=CU1
      	  R=R1*(CU1/C1)**(2./(GN1-1.))
            P=P1*(CU1/C1)**(2.*GN1/(GN1-1.))
            CALL VNENER(R,P,GN1,E)
      	ELSE
      	  P=P2
      	  U=U2
      	  R=GN1*P/(C11*C11)
            CALL VNENER(R,P,GN1,E)
          END IF
        END IF
      ELSE
        P=P2
        U=U2
        R=R2
        E=E2
      END IF
C-----Shock - Depresision Waves-------------------------------------------
      CASE("sh-dp")
      D1=U1-A1/R1
      D2=U2+C2
      C22=C2-(GN2-1.)*(U2-U)/2.
      IF(C22 < Cmin) C22=Cmin
      D22=U+C22
      IF(U >= 0.)THEN
        IF(D1 >= 0.)THEN
          P=P1
          U=U1
          R=R1
          E=E1
        ELSE
          RA=A1-R1*(U1-U)
          IF(RA < Amin) RA=Amin
      	R=R1*A1/RA
          CALL VNENER(R,P,GN1,E)
        END IF
      ELSE
        IF(D2 <= 0.)THEN
          P=P2
          U=U2
          R=R2
          E=E2
        ELSE
          IF(D22 < 0.)THEN
            CU2=(2.*C2-(GN2-1.)*U2)/(GN2+1.)
            IF(CU2 < Cmin) CU2=Cmin
      	  U=-CU2
      	  R=R2*(CU2/C2)**(2./(GN2-1.))
            P=P2*(CU2/C2)**(2.*GN2/(GN2-1.))
            CALL VNENER(R,P,GN2,E)
      	ELSE
      	  R=GN2*P/(C22*C22)
            CALL VNENER(R,P,GN2,E)
          END IF
        END IF
      END IF
C-----Depresision - Shock Waves-------------------------------------------
      CASE("dp-sh")
      D1=U1-C1
      C11=C1+(GN1-1.)*(U1-U)/2.
      IF(C11 < Cmin) C11=Cmin
      D11=U-C11
      D2=U2+A2/R2
      IF(U >= 0.)THEN
        IF(D1 >= 0.)THEN
          P=P1
          U=U1
          R=R1
          E=E1
        ELSE
          IF(D11 > 0.)THEN
            CU1=(2.*C1+(GN1-1.)*U1)/(GN1+1.)
            IF(CU1 < Cmin) CU1=Cmin
      	  U=CU1
      	  R=R1*(CU1/C1)**(2./(GN1-1.))
            P=P1*(CU1/C1)**(2.*GN1/(GN1-1.))
            CALL VNENER(R,P,GN1,E)
      	ELSE
      	  R=GN1*P/(C11*C11)
            CALL VNENER(R,P,GN1,E)
          END IF
        END IF
      ELSE
        IF(D2 < 0.)THEN
          P=P2
          U=U2
          R=R2
          E=E2
        ELSE
          RA=A2+R2*(U2-U)
          IF(RA < Amin) RA=Amin
      	R=R2*A2/RA
          CALL VNENER(R,P,GN2,E)
        END IF
      END IF
C-----Tangential - Shock Waves--------------------------------------------
      CASE("tn-sh")
      D2=U2+A2/R2
      IF(U1 >= 0.)THEN
        P=P1
        U=U1
        R=R1
        E=E1
      ELSE
        IF(D2 < 0.)THEN
          P=P2
          U=U2
          R=R2
          E=E2
        ELSE
          P=P1
          U=U1
          RA=A2+R2*(U2-U1)
          IF(RA < Amin) RA=Amin
      	R=R2*A2/RA
          CALL VNENER(R,P,GN2,E)
        END IF
      END IF
C-----Shock - Tangential Waves--------------------------------------------
      CASE("sh-tn")
      D1=U1-A1/R1
      IF(U2 >= 0.)THEN
        IF(D1 >= 0.)THEN
          P=P1
          U=U1
          R=R1
          E=E1
        ELSE
          P=P2
      	U=U2
      	RA=A1-R1*(U1-U2)
          IF(RA < Amin) RA=Amin
      	R=R1*A1/RA
          CALL VNENER(R,P,GN1,E)
        END IF
      ELSE
        P=P2
        U=U2
        R=R2
        E=E2
      END IF
C-----Shock - Shock Waves-------------------------------------------------
      CASE("sh-sh")
      D1=U1-A1/R1
      D2=U2+A2/R2
      IF(U >= 0.)THEN
        IF(D1 >= 0.)THEN
          P=P1
          U=U1
          R=R1
          E=E1
        ELSE
          RA=A1-R1*(U1-U)
          IF(RA < Amin) RA=Amin
      	R=R1*A1/RA
          CALL VNENER(R,P,GN1,E)
        END IF
      ELSE
        IF(D2 < 0.)THEN
          P=P2
          U=U2
          R=R2
          E=E2
        ELSE
          RA=A2+R2*(U2-U)
          IF(RA < Amin) RA=Amin
      	R=R2*A2/RA
          CALL VNENER(R,P,GN2,E)
        END IF
      END IF
C-----Vacuum Case---------------------------------------------------------
      CASE("vacum")
      Pt=Pmin
      Ut=(U1-2.*C1*(DABS(Pt/P1)**Gd1-1.)/(GN1-1.)+
     *    U2+2.*C2*(DABS(Pt/P2)**Gd2-1.)/(GN2-1.))/2.
      D1=U1-C1
      D11=U1+C1*2./(GN1-1.)
      D2=U2+C2
      D22=U2-C2*2./(GN2-1.)
      IF(Ut >= 0.)THEN
        IF(D1 >= 0.)THEN
          P=P1
          U=U1
          R=R1
          E=E1
        ELSE
          IF(D11 < Ut .and. D11 > 0.)THEN
            CU1=(2.*C1+(GN1-1.)*U1)/(GN1+1.)
            IF(CU1 < Cmin) CU1=Cmin
      	  U=CU1
      	  R=R1*(CU1/C1)**(2./(GN1-1.))
            P=P1*(CU1/C1)**(2.*GN1/(GN1-1.))
            CALL VNENER(R,P,GN1,E)
      	ELSE
            P=Pmin
      	  U=Ut
      	  R=Rmin
            E=Emin
          END IF
        END IF
      ELSE
        IF(D2 <= 0.)THEN
          P=P2
          U=U2
          R=R2
          E=E2
        ELSE
          IF(D22 > Ut .and. D22 < 0. )THEN
            CU2=(2.*C2-(GN2-1.)*U2)/(GN2+1.)
            IF(CU2 < Cmin) CU2=Cmin
      	  U=-CU2
      	  R=R2*(CU2/C2)**(2./(GN2-1.))
            P=P2*(CU2/C2)**(2.*GN2/(GN2-1.))
            CALL VNENER(R,P,GN2,E)
      	ELSE
            P=Pmin
      	  U=Ut
      	  R=Rmin
            E=Emin
          END IF
        END IF
      END IF
C-----DEFAULT Case---------------------------------------------------------
      CASE DEFAULT
      RG=DABS(R1*G1*R2*G2)**(3/2)
      DN=2.*(DSQRT(DABS(R1*G1))*R2*G2+DSQRT(DABS(R2*G2))*R1*G1)
      IF(DN < Amin) DN=Amin
      Pt=((R1*R2*(U1-U2)*G1*G2+
     &DSQRT(DABS(4.*(P1+P2)*RG+4.*R1*R2*G1*G2*(P2*R1*G1+P1*R2*G2)+
     &        (R1*R2*U1*G1*G2-R1*R2*U2*G1*G2)**2)))/DN)**2
      Ut=(U1-2.*C1*(DABS(Pt/P1)**Gd1-1.)/(GN1-1.)+
     *    U2+2.*C2*(DABS(Pt/P2)**Gd2-1.)/(GN2-1.))/2.
      D1=U1-C1
      D2=U2+C2
      IF(Ut >= 0.)THEN
        IF(D1 >= 0.)THEN
          P=P1
          U=U1
          R=R1
          E=E1
        ELSE
          P=Pt
          U=Ut
      	R=GN1*Pt/(C1*C1)
          CALL VNENER(R,P,GN1,E)
        END IF
      ELSE
        IF(D2 <= 0.)THEN
          P=P2
          U=U2
          R=R2
          E=E2
        ELSE
          P=Pt
      	U=Ut
      	R=GN2*Pt/(C2*C2)
          CALL VNENER(R,P,GN2,E)
        END IF
      END IF
      END SELECT
C-------------------------------------------------------------------------
11    CONTINUE
      IF(R < Rmin )R=Rmin
      IF(P < Pmin )P=Pmin
      IF(E < Emin )E=Emin
      if(vak_case) then

!  Если газ равномерно расширяется,
!   n~= n0*(1-DV*t/h) 
!                tau
!   rho~n
!
!  Пусть адиабатически, тогда
! 
!   p~rho**gamma
!  
!   p = p0*(r/r0)**gamma
!
!        write(*,*)'U=',U,'Uav=',(0.5*(U1+U2))
        C1=0.5
        C2=0.5
C        C1=(0.5d0+u2*tau)/(1.0d0+(u2-u1)*tau)
C        C2=(0.5d0-u1*tau)/(1.0d0+(u2-u1)*tau)
        if((C1.lt.0.0).or.(C1.gt.1.0)) then
          write(*,*)'call to riemann bug:'
          write(*,*)'C1=',c1
          write(*,*)'U1=',U1
          write(*,*)'U2=',U2
C          write(*,*)'tau=',tau
          stop
        end if
        if((C2.lt.0.0).or.(C2.gt.1.0)) then
          write(*,*)'call to riemann bug:'
          write(*,*)'C2=',c2
          write(*,*)'U1=',U1
          write(*,*)'U2=',U2
C          write(*,*)'tau=',tau
          stop
        end if
          
          
!       U=0.5*(U1+U2)
!       R=0.5*(R1+R2)!*(1.0-tau*(U2-U1)*0.5)
!       P=0.5*(P1+P2)!(1.0-tau*(U2-U1)*0.5)
       u=c1*u1+c2*u2
       r=c1*r1+c2*r2
       p=c1*p1+c2*p2
        
!       p=0.5*(p1+p2)*((r/(0.5*(r1+r2)))**(0.5*(GN1+GN2)))
!       if(dabs((p-p1)/p1).gt.0.1) then
!         write(*,*)p1,p2,p
!       end if
      			call VNENER(R,P,GN2,E)
      end if
!      if(var.eq."dp-dp") then
!        write(*,*)'dp-dp U=',U,'Uav=',(0.5*(U1+U2))
!      end if
      RETURN

      END SUBROUTINE RIEMANN

C-------------------------------------------------------------------------
C	Calculation value of E through values pressure - P and dencity - R
C-------------------------------------------------------------------------
	SUBROUTINE VNENER(R,P,G,E)
	IMPLICIT REAL*8 (a-h,o-z)
	INTENT(IN):: R, P, G
	INTENT(OUT):: E
	COMMON /MIN/ Rmin,Emin,Cmin,Pmin
	IF(R <= Rmin)THEN
	  E=P/(G-1.)/Rmin
	ELSE
	  E=P/(G-1.)/R
	END IF
	IF(E <= Emin) E=Emin
	RETURN
	END  SUBROUTINE VNENER
C---------------------------------------------------------------------------------
C	Calculation values of P,C through values dencity - R and internal energy - E
C---------------------------------------------------------------------------------
	SUBROUTINE URSOST(R,E,G,C,P,PR,PE)
	IMPLICIT REAL*8 (a-h,o-z)
	INTENT(IN):: R, E, G
	INTENT(OUT):: C, P, PR, PE
	COMMON /MIN/ Rmin,Emin,Cmin,Pmin
	CC=G*(G-1.)*E
	if(cc.lt.1.e-35)then
		write(*,*)'9sq cc=',cc, ' E=',E,' G=',G
		cc=1.e-35
	end if
	C = DSQRT(CC)
	IF(C <= Cmin) C=Cmin
	P =(G-1.)*R*E
	IF(P <= Pmin) P=Pmin
	PR=(G-1.)*E
	IF(PR < (G-1.)*Emin) PR=(G-1.)*Emin
	PE=(G-1.)*R
	IF(PE < (G-1.)*Rmin) PE=(G-1.)*Rmin
	RETURN
	END SUBROUTINE URSOST
C-------Podprogramma vychislenija matricy OM iz levykh sobstv.vektorov
C-------matricy Df/Du iskhodnoj sist. i obratnoj matricy OMO ---------
	SUBROUTINE OMEGA(G,V0,OM,OMO,IU) 
	IMPLICIT REAL*8 (a-h,o-z)
	DIMENSION V0(IU),OM(IU,IU),OMO(IU,IU)
	INTENT(IN):: G, V0, IU
	INTENT(OUT):: OM, OMO 
	COMMON /MIN/ Rmin,Emin,Cmin,Pmin
C----------------------------------------------------------------------
	R=V0(1)
	U=V0(2)
	E=V0(3)
	IF(R < Rmin )THEN
	  C=Cmin; U=0.
	ELSE
	  CALL URSOST(R,E,G,C,P,PR,PE)
	END IF
C----------------------------------------------------------------------
	GM=G-1.
	OM(1,1)=U*U/C/2.-U/GM
	OM(1,2)=1./GM-U/C
	OM(1,3)=1./C
	OM(2,1)=U*U/C/2.-C/GM
	OM(2,2)=-U/C
	OM(2,3)=1./C
	OM(3,1)=U*U/C/2.+U/GM
	OM(3,2)=-1./GM-U/C
	OM(3,3)=1./C
	OMO(1,1)=GM/C/2.
	OMO(1,2)=-GM/C
	OMO(1,3)=OMO(1,1)
	OMO(2,1)=(C+U)*GM/C/2.
	OMO(2,2)=-U*GM/C
	OMO(2,3)=(U-C)*GM/C/2.
	OMO(3,1)=C/2.+U*GM/2.+U*U*GM/C/4.
	OMO(3,2)=-U*U*GM/C/2.
	OMO(3,3)=C/2.-U*GM/2.+U*U*GM/C/4.
	RETURN
	END SUBROUTINE OMEGA
C-------Podpr. vychislenija sobstvennuh znacheni matrizi Jacobi
	SUBROUTINE SZN(G,V,IU,SIG,S)  
	IMPLICIT REAL*8 (a-h,o-z)
	DIMENSION  V(IU),S(IU),V0(IU)
	INTENT(IN)::  G, V, IU, SIG
	INTENT(OUT):: S 
	COMMON /MIN/ Rmin,Emin,Cmin,Pmin
C----------------------------------------------------------------------
	R = V(1)
	U = V(2)
	E = V(3)
	IF(R <= Rmin )THEN
	  C=Cmin; U=0.
	ELSE
	  CALL URSOST(R,E,G,C,P,PR,PE)
	END IF
	S(1)=(U+C)*SIG
	S(2)=U*SIG
	S(3)=(U-C)*SIG
	RETURN
	END SUBROUTINE SZN
C-------Podprogr. vychislenija parametrov gibridnosti G -------------
	SUBROUTINE GIBR(M,BG,IU,EPS,UN,UNM,UNP,UFM,UPM,UPP,OM,S,	  
     &				GIB,RUS,D0,C0)
	IMPLICIT REAL*8 (a-h,o-z)
	DIMENSION UN(IU),UNM(IU),UNP(IU),UFM(IU),UPM(IU),UPP(IU),
     *		  OM(IU,IU),GIB(IU),S(IU),RUS(IU),D0(IU),C0(IU),
     *		  DM1(IU),D1(IU)
	INTENT(IN)::  M,BG,IU,EPS,UN,UNM,UNP,UFM,UPM,UPP,OM,S
	INTENT(OUT):: GIB,RUS,D0,C0
	D0 = MATMUL(OM,(UNM-UNP))
	C0 = MATMUL(OM,(UPM-UPP))
	DO I=1,IU
	   GIB(I) = D0(I)*D0(I) - C0(I)*C0(I)
	END DO
	RUS=0.
	IF(BG > 2.5 .AND. M > 0)THEN
	  DM1 = MATMUL(OM,(UN-UNM))
	  D1 = MATMUL(OM,(UNP-UFM))
	  DO I=1,IU
		 SM = DABS(S(I))
		 IF(SM < EPS .OR. DABS(D0(I)) < EPS)THEN
		   RUS(I)=0.
		 ELSE
		   IF(S(I) > 0.)THEN
			 D1(I)	= D1(I)/D0(I)
			 DM1(I) = DM1(I)/D0(I)
			 W = SM*(1.+D1(I))/2.+SM*SM*(1.-D1(I))/2.+
     *			  SM*(SM*SM-1.)*(DM1(I)+D1(I)-2.)/6.
			 IF(W>=EPS.AND.W<=1.-EPS.AND.DM1(I)>=0..AND.D1(I)>=0.)
     *		 THEN
			   RUS(I)=1.
			 ELSE
			   RUS(I)=0.
			 END IF
		   ELSE
			 D1(I)	= DM1(I)/D0(I)
			 DM1(I) = D1(I)/D0(I)
			 W = SM*(1.+D1(I))/2.+SM*SM*(1.-D1(I))/2.+
     *			  SM*(SM*SM-1.)*(DM1(I)+D1(I)-2.)/6.
			 IF(W>=EPS.AND.W<=1.-EPS.AND.DM1(I)>=0..AND.D1(I)>=0.)
     *		 THEN
			   RUS(I)=1.
			 ELSE
			   RUS(I)=0.
			 END IF
		   END IF
		 END IF 	   
	  END DO
	END IF
	D0 = -D0
	C0 = -C0
	RETURN
	END SUBROUTINE GIBR
C-------Podpr. vychislenija parametrov B,C,D sglaghivajusch.operatora-
	SUBROUTINE BCD(RUS,S,IU,GIB,B,CC)
	IMPLICIT REAL*8 (a-h,o-z)
	DIMENSION  S(IU),B(IU),CC(IU),GIB(IU),RUS(IU)
	INTENT(IN)::  RUS, S, IU, GIB 
	INTENT(OUT):: B, CC
	DO I=1,IU
	   SM=DABS(S(I))
	   IF(DABS(RUS(I)) > 0.)THEN
		  B(I)=(SM-1.)*(2.*SM-1.)/6.
		  CC(I)=(SM*SM-1.)/6.
	   ELSE 
		  IF(GIB(I) > 0.)THEN
			 B(I)=(1.-SM)**2/2.
			 CC(I)=(SM-1.)/2.
		  ELSE
			 B(I)=SM*(SM-1.)/2.
			 CC(I)=0.
		  END IF
	   END IF
	END DO
	RETURN
	END SUBROUTINE BCD
C----------------------------------------------------------------------
C	Numerical Scheme
C----------------------------------------------------------------------
	SUBROUTINE SKHEMA(N,U1,UM1,UN,UN1,MK,S,AL)
	IMPLICIT REAL*8 (a-h,o-z)
	INTENT(IN):: N,MK,U1,UM1,UN,S,AL
	INTENT(OUT):: UN1
	DIMENSION UM1(MK),UN(MK),UN1(MK)
	IF( AL >= 0.)THEN
	 UN1(1)=U1
	 DO M=2,MK
	   UN1(M)=UN1(M-1)*2.*(S-1.)*S/(1.+3.*S+2.*S*S)+
     *		  UN(M-1)*2.*S/(S+1.)+UN(M)*(-2.+4./(1.+S))+
     *		  UM1(M)*(S-1.)/(1.+3.*S+2.*S*S)
	   UMN=UN1(M)
C----------------------------------------------------------------------
	   IF(UN(M-1) <= UN(M))THEN
		  IF(UN1(M) < UN(M-1) .OR. UN1(M) > UN(M)) 
     *		 UN1(M)=UN(M-1)*2.*S*S/(S+1.)+2.*UN(M)*(1.-S)+
     *				UM1(M)*(S-1.)/(1.+S)
			 IF(UN1(M) < UN(M-1) .OR. UN1(M) > UN(M))THEN
			   IF(UMN < UN(M-1))UN1(M)=UN(M-1)
			   IF(UMN > UN(M))UN1(M)=UN(M)
			 END IF 
	   ELSE
		  IF(UN1(M) < UN(M) .OR. UN1(M) > UN(M-1))
     *		 UN1(M)=UN(M-1)*2.*S*S/(S+1.)+2.*UN(M)*(1.-S)+
     *				UM1(M)*(S-1.)/(1.+S)
			IF(UN1(M) < UN(M) .OR. UN1(M) > UN(M-1))THEN
			   IF(UMN > UN(M-1))UN1(M)=UN(M-1)
			   IF(UMN < UN(M))UN1(M)=UN(M)
			END IF 
	   END IF
C----------------------------------------------------------------------
	   IF(N == 1)THEN
		  UN1(M)=UN(M)+S*(UN(M-1)-UN(M))
	   END IF
	 END DO
	ELSE
	 UN1(MK)=U1
	 DO M=MK-1,1,-1
	   UN1(M)=UN1(M+1)*2.*(S-1.)*S/(1.+3.*S+2.*S*S)+
     *		  UN(M+1)*2.*S/(S+1.)+UN(M)*(-2.+4./(1.+S))+
     *		  UM1(M)*(S-1.)/(1.+3.*S+2.*S*S)
	   UMN=UN1(M)
C----------------------------------------------------------------------
	   IF(UN(M) <= UN(M+1))THEN
		  IF(UN1(M) < UN(M) .OR. UN1(M) > UN(M+1)) 
     *		 UN1(M)=UN(M+1)*2.*S*S/(S+1.)+2.*UN(M)*(1.-S)+
     *				UM1(M)*(S-1.)/(1.+S)
			 IF(UN1(M) < UN(M) .OR. UN1(M) > UN(M+1))THEN
			   IF(UMN < UN(M))UN1(M)=UN(M)
			   IF(UMN > UN(M+1))UN1(M)=UN(M+1)
			 END IF 
	   ELSE
		  IF(UN1(M) < UN(M+1) .OR. UN1(M) > UN(M))
     *		 UN1(M)=UN(M+1)*2.*S*S/(S+1.)+2.*UN(M)*(1.-S)+
     *				UM1(M)*(S-1.)/(1.+S)
			IF(UN1(M) < UN(M+1) .OR. UN1(M) > UN(M))THEN
			   IF(UMN < UN(M+1))UN1(M)=UN(M+1)
			   IF(UMN > UN(M))UN1(M)=UN(M)
			END IF 
	   END IF
C----------------------------------------------------------------------
	   IF(N == 1)THEN
		  UN1(M)=UN(M)+S*(UN(M+1)-UN(M))
	   END IF
	 END DO
	END IF
	RETURN
	END
