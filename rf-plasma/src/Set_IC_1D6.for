      Program INIT_SET  
C-----------------------------------------------------------------------
c    This routine sets mesh and initial conditions for the 1D CFD code  
c    and writes them to an output file "Data1D.in"
C-----------------------------------------------------------------------
	IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER*4 (i-n)
      PARAMETER (IU = 7)
      ALLOCATABLE :: V(:,:) 
C-----------------------------------------------------------------------
	common/bconst/cs,csv,cze,cme,ckb,ci,cih,ca2,ce0,cedm,cmdm
	common/bconatm/cmb,cnw,cbv,cno2,czi,czb,cne0,cni0,ct0,cti
	common/batmosf/zsi,zn,zni0,zne0,zn00,zne,zni,zn0,zno2
	common/bimpuls/	chi,csi,cpi,cep,z0
	common/bconist/cf1v,cvei,cve0,se00,vei0,ve00,chn0
	common/bistoch/tjx,tjz,sie,pesie,vee,pevee,tje,petje,tjeh,petjeh,
     *tjn,petjn,tjd,petjd,tjp,petjp,vei,pnvei,pevei,ve0,peve0,fie,pefie,
     *chn,pechn,dee,pedee
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
c	WRITE(*,*)'Konec uchastka integrirovanija TK=?'
c	READ(42,*)TK
c	WRITE(*,*)'tk=',TK
	WRITE(*,*)'Nachalnyj shag integrirowanija po vremeni DT=?'
	READ(42,*)DT
	WRITE(*,*)'dt=',DT
	WRITE(*,*)'Nachalnaya vysota Z0=?'
	READ(42,*)Z0
	WRITE(*,*)'Z0=',Z0
	cf1v=csv*cs*cze
	WRITE(*,*)'cf1v=',cf1v
	WRITE(46,*)cf1v,'    cf1v'
	cvei=(((2.*sqrt(3.1415926*3./cme)*cze)*cze)*cze)*cze
c	cvei=2.*sqrt(3.1415926*3./cme)*cze
	WRITE(*,*)'cvei=',cvei
	WRITE(46,*)cvei,'    cvei'
	cve0=5.33333/sqrt(3.*3.1415926*cme)
	WRITE(*,*)'cve0=',cve0
	WRITE(46,*)cve0,'    cve0'
	WRITE(*,*)'Otnoshenie zaryada elektrona k ego masse cedm=?'
c	cedm=cze/cme
	READ(42,*)cedm
	WRITE(*,*)'cedm=',cedm
	WRITE(*,*)'Otnoshenie massy elektr. k masse iona vozdukha cmdm=?'
	READ(42,*)cmdm
c	cmdm=cme/cmb
	WRITE(*,*)'cmdm=',cmdm
	WRITE(*,*)'Sdvig nachalnogo raspredelen.(v dolyakh ot XXK ) x0s=?'
	READ(42,*)x0s
	WRITE(*,*)'x0s=',x0s
C-----------------------------------------------------------------------
      WRITE(*,*)'Number of the grid steps in one dimension MK=?'
      READ(43,*)MK
	write(*,*)'MK=',mk
      WRITE(*,*)' Adiabatic constant Ak MK=?'
      READ(43,*)Ak
	write(*,*)'Ak=',Ak
      WRITE(*,*)'Region size XXK=?'
      READ(43,*)XXK
	write(*,*)'XXK=',XXK
      WRITE(*,*)'Type of Boundary Conditions on the left side:'
      WRITE(*,*)'NL = 0 is the Reflective BC'
      WRITE(*,*)'NL = 1 is the Transparent BC'
      WRITE(*,*)'NL = 2 is the Periodic BC'
      READ(43,*)NL
	write(*,*)'NL=',NL
      WRITE(*,*)'Type of Boundary Conditions on the right side:'
      WRITE(*,*)'NR = 0 is the Reflective BC'
      WRITE(*,*)'NR = 1 is the Transparent BC'
      WRITE(*,*)'NR = 2 is the Periodic BC'
      READ(43,*)NR
	write(*,*)'NR=',NR
      WRITE(*,*)'Currant number AKYR=(0.9-0.99)'
      READ(43,*)AKYR
	write(*,*)'AKYR=',AKYR
      WRITE(*,*)'Input: BG=1, if you want 1 order of Scheme accuracy'
      WRITE(*,*)'       BG=2, if you want 2 order of Scheme accuracy'
      WRITE(*,*)'       BG=3, if you want 3 order of Scheme accuracy'
      READ(43,*)BG
	write(*,*)'BG=',BG
      WRITE(*,*)'Parameter of the Riemann Problem RID=0.001'
      READ(43,*)RID
	write(*,*)'RID=',RID
	WRITE(*,*)'Minimal Electron Concentration Rmin<=1E-10?'			
	READ(43,*)Rmin 
	write(*,*)'Rmin=',Rmin  
	WRITE(*,*)'Minimal Electron Internal Energy Emin<=1E-10?'			
	READ(43,*)Emin
	write(*,*)'Emin=',Emin     
      WRITE(*,*)'Parameter of the Scheme Hybridity EPS=0.001'
      READ(43,*)EPS
	write(*,*)'EPS=',EPS
C-----------------------------------------------------------------------
C	Все, что задается ниже, можно менять по желанию!
C	НАПРИМЕР: 
C-----------------------------------------------------------------------
      WRITE(*,*)'Electron Velocity in Z direction = ?'
      Uz=0. !1. !3.e10 !3.e8 !3.e4 !1.e8 !10. !1. !0. !cs !2. !
	write(*,*)'Uz=',Uz
      WRITE(*,*)'Internal region values:'
      WRITE(*,*)'Electron Concentration =?'
      R1=cne0          !Koncentraciya electronov  Ne
	write(*,*)'R1=',R1
      WRITE(*,*)'Electron Velocity in X direction = ?'
      Ux=0.
	write(*,*)'Ux=',Ux
      WRITE(*,*)'Electron Internal Energy =?'
      E1=ce0 !1.5*ckb*cti   !Vnutrennyaya energiya electronov  Ee
	write(*,*)'E1=',E1
      WRITE(*,*)'Electric field in X direction = ?'
      Ex=cep            !Amplituda electricheskogo polya po osi X Ex
	write(*,*)'Ex=',Ex
      WRITE(*,*)'Magnetic field in Y direction = ?'
      Hy=cep/csv       !Amplituda magnitnogo polya po osi Y  Hy
	write(*,*)'Hy=',Hy
      WRITE(*,*)'Electric field in Z direction = ?'
      Ez=0. ! 1. !1000. !0.            !Amplituda electricheskogo polya po osi Z Ez
	write(*,*)'Ez=',Ez
      WRITE(*,*)'External region values:'
      WRITE(*,*)'External region Electron Concentration =?:'
      R2=cne0 !50. !99. !50. !1. !Rmin !cne0 !         !Koncentraciya electronov  Ne
	write(*,*)'R2=',R2
      E2=R1*E1/R2
	write(*,*)'E2=',E2
C-----------------------------------------------------------------------
	ALLOCATE(V(MK,IU), stat=istat) 
	IF(istat /=0) STOP "***ERROR ON ALLOCATE***"  
C-----------------------------------------------------------------------
!	V(:,:) = 0.
	dx0=cs/chi
	dxs=cs/csi/cpi
	DO M = 1, MK
	xxx=(m-1.)/(mk-1.)*xxk
	yyy=xxx-x0s*xxk
	ns=yyy/dxs
	dxsd=yyy-ns*dxs
	dx0d=dxsd/dx0
c	write(*,*)dx0,dxs,m,mk,xxk,xxx,yyy,ns,dxsd,dx0d,cpi,x0s
	if(dx0d.ge.0..and.dx0d.le.1..and.ns.lt.cpi)then
           V(M,1) = R1
           V(M,2) = Ux
	     V(M,3) = Uz
           V(M,4) = E1
           V(M,5) = 0 !-Ex*dsin(6.2831858*dx0d)
           V(M,6) = 0 !-Hy*dsin(6.2831858*dx0d)
           V(M,7) = Ez  
c	write(*,*)dx0,dxs,m,mk,xxk,xxx,yyy,ns,dxsd,dx0d,cpi,x0s,v(m,5)
c	stop    
         ELSE
           V(M,1) = R2
           V(M,2) = Ux
	     V(M,3) = Uz
           V(M,4) = E2
           V(M,5) =0. ! ex !Ex
           V(M,6) =0. ! hy !Hy
           V(M,7) =0. ! ez !Ez
        END IF
      END DO
C-----------------------------------------------------------------------
	OPEN(1,file="Data1D.in",status="replace",iostat=istat)
	IF(istat /= 0) STOP "***Cannot open file***"
	WRITE(1,*)MK,'    Number of the grid steps'		
	WRITE(1,*)Ak,' Adiabatic constant'			     
	WRITE(1,*)XXK,'   Region size in Z direction'			
	WRITE(1,*)NL,'    BC on the left side  0 - Rf, 1 - Tr, 2- Perd'		
	WRITE(1,*)NR,'    BC on the right side 0 - Rf, 1 - Tr, 2- Perd'			
	WRITE(1,*)AKYR,'  Currant number'		
	WRITE(1,*)BG,  '  Approximation order'			
	WRITE(1,*)RID, '  Parameter of the Riemann Problem - RID'			
	WRITE(1,*)Rmin, ' Electron Concentration - Rmin'			
	WRITE(1,*)Emin, ' Minimal Electron Internal Energy - Emin'			
	WRITE(1,*)EPS, '  Parameter of the Scheme Hybridity - EPS'			
      WRITE(1,*)cze, '   Electron Charge'
      WRITE(1,*)cs, '   Velocity of Light'
      WRITE(1,*)csv, '   Impedance of Vacuum'
	WRITE(1,*)
	DO J = 1, IU
	   WRITE(1,*) (V(M,J), M =1,MK) !Data arrays of the Concentration, Velocity and Internal Energy 
	   WRITE(1,*)
	END DO
	CLOSE(1)
	STOP
      END Program INIT_SET
C=======================================================================











