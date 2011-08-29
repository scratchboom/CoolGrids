
#include "GridsCommon.hpp"

using namespace std;

int main(){

	omp_set_num_threads(4);

	VtiSaver3D vtiSaver;
	CImgSaver2D cimgSaver;


    double EPS0=8.854E-12;// Ф/м
    double MU0=1.256E-6;// Гн/м
    double C=299792458.0;// м/с    C=sqrt(1/EPS0/MU0)
    double VACUUM_WAVE_RESISTIVITY=sqrt(MU0/EPS0);// волновое сопротивление вакуума [Ом]



    double Nx=128;
    double Ny=81;
    double Nz=80;
    double Nt=2000;

    double SOURCE_IX = (int)(Nx*0.1)+0.5;
    double SOURCE_IY = (int)(Ny*0.5)+0.5;
    double SOURCE_IZ = Nz/2;

    double J_SOURCE_IX = (int)(Nx*0.5);
    double J_SOURCE_IY = (int)(Ny*0.5);
    double J_SOURCE_IZ = (int)(Nz*0.5)+0.5;


    double frec=3e8;

    double IMPULSE_TIME_WIDTH=1.0/frec;

    double E_AMPLITUDE=1.0;
    double H_AMPLITUDE=E_AMPLITUDE/VACUUM_WAVE_RESISTIVITY;

    double J_AMPLITUDE=1.0;

	cimgSaver.setValueRange(-2.0*H_AMPLITUDE , 2.0*H_AMPLITUDE);

    double Sx=10.0;
	double Sy= (Sx/Nx) * Ny;
	double Sz= (Sx/Nx) * Nz;
	double T=1.0;

	double dx=Sx/Nx;
	double dy=Sy/Ny;
	double dz=Sz/Nz;
	//double DT=T/Nt;

    double DT=0.5*dx/C;

//    cout << IMPULSE_TIME_WIDTH*C << endl;
//    cout << IMPULSE_TIME_WIDTH/DT << endl;
//    exit(0);


	TimedGrid3D<double> Ex("Ex");
	TimedGrid3D<double> Ey("Ey");
	TimedGrid3D<double> Ez("Ez");

	TimedGrid3D<double> Hx("Hx");
	TimedGrid3D<double> Hy("Hy");
	TimedGrid3D<double> Hz("Hz");


    TimedGrid3D<double> epsilon("epsilon");
	TimedGrid3D<double> mu("mu");
	TimedGrid3D<double> sigma("sigma");


	//building grids
	Ex.setRangeT(0,T);
	Ex.setRangeX(0.5*dx , Sx-0.5*dx);
	Ex.setRangeY(0,Sy);
	Ex.setRangeZ(0,Sz);

	Ex.setIndexRangeT(0,Nt);
	Ex.setIndexRangeX(0.5 , Nx-0.5);
	Ex.setIndexRangeY(0,Ny);
	Ex.setIndexRangeZ(0,Nz);

	Ex.setLayersCountToMaintain(7);
	Ex.build();


	Ey.setRangeT(0,T);
	Ey.setRangeX(0 , Sx);
	Ey.setRangeY(0.5*dy , Sy-0.5*dy);
	Ey.setRangeZ(0 , Sz);

	Ey.setIndexRangeT(0,Nt);
	Ey.setIndexRangeX(0,Nx);
	Ey.setIndexRangeY(0.5,Ny-0.5);
    Ey.setIndexRangeZ(0,Nz);

    Ey.setLayersCountToMaintain(7);
	Ey.build();


	Ez.setRangeT(0,T);
	Ez.setRangeX(0 , Sx);
	Ez.setRangeY(0 , Sy);
	Ez.setRangeZ(0.5*dz , Sz-0.5*dz);

	Ez.setIndexRangeT(0,Nt);
	Ez.setIndexRangeX(0 , Nx);
	Ez.setIndexRangeY(0 , Ny);
	Ez.setIndexRangeZ(0.5 , Nz-0.5);

	Ez.setLayersCountToMaintain(7);
	Ez.build();



	Hx.setRangeT(0,T);
	Hx.setRangeX(0 , Sx);
	Hx.setRangeY(0.5*dy , Sy-0.5*dy);
	Hx.setRangeZ(0.5*dz , Sz-0.5*dz);

	Hx.setIndexRangeT(0.5,Nt+0.5);
	Hx.setIndexRangeX(0   , Nx);
	Hx.setIndexRangeY(0.5 , Ny-0.5);
	Hx.setIndexRangeZ(0.5 , Nz-0.5);

	Hx.setLayersCountToMaintain(3);
	Hx.build();




	Hy.setRangeT(0,T);
	Hy.setRangeX(0.5*dx , Sx-0.5*dx);
	Hy.setRangeY(0      , Sy       );
	Hy.setRangeZ(0.5*dz , Sz-0.5*dz);

	Hy.setIndexRangeT(0.5,Nt+0.5);
	Hy.setIndexRangeX(0.5 , Nx-0.5);
	Hy.setIndexRangeY(0   , Ny    );
    Hy.setIndexRangeZ(0.5 , Nz-0.5);

    Hy.setLayersCountToMaintain(3);
	Hy.build();




	Hz.setRangeT(0,T);
	Hz.setRangeX(0.5*dx , Sx-0.5*dx);
	Hz.setRangeY(0.5*dy , Sy-0.5*dy);
	Hz.setRangeZ(0      , Sz);

	Hz.setIndexRangeT(0.5,Nt+0.5);
	Hz.setIndexRangeX(0.5 , Nx-0.5);
	Hz.setIndexRangeY(0.5 , Ny-0.5);
	Hz.setIndexRangeZ(0   , Nz    );

	Hz.setLayersCountToMaintain(3);
	Hz.build();



	epsilon.setRangeT(0,T);
	epsilon.setRangeX(0,Sx);
	epsilon.setRangeY(0,Sy);
	epsilon.setRangeZ(0,Sz);

	epsilon.setIndexRangeT(0,Nt);
	epsilon.setIndexRangeX(0,Nx);
	epsilon.setIndexRangeY(0,Ny);
	epsilon.setIndexRangeZ(0,Nz);

	epsilon.build();

	epsilon[0].iterateWhole(GRID3D_ITERATOR{
	    epsilon(0,ix,iy,iz) = 1.0;
	});




	mu.setRangeT(0,T);
	mu.setRangeX(-0.5*dx,Sx+0.5*dx);
	mu.setRangeY(-0.5*dy,Sy+0.5*dy);
	mu.setRangeZ(-0.5*dz,Sz+0.5*dz);

	mu.setIndexRangeT(0,Nt);
	mu.setIndexRangeX(-0.5,Nx+0.5);
	mu.setIndexRangeY(-0.5,Ny+0.5);
	mu.setIndexRangeZ(-0.5,Nz+0.5);

	mu.build();

	mu[0].iterateWhole(GRID3D_ITERATOR{
	    mu(0,ix,iy,iz)=1.0;
	});

	sigma.setRangeT(0,T);
	sigma.setRangeX(0,Sx);
	sigma.setRangeY(0,Sy);
	sigma.setRangeZ(0,Sz);

	sigma.setIndexRangeT(0,Nt);
	sigma.setIndexRangeX(0,Nx);
	sigma.setIndexRangeY(0,Ny);
	sigma.setIndexRangeZ(0,Nz);

	sigma.build();

	sigma[0].iterateWhole(GRID3D_ITERATOR{
	    sigma(0,ix,iy,iz)=0.0;
	});




	TimedGrid3D < Vector5D > U("U");

	TimedGrid3D < Vector5D > Fx("Fx");
	TimedGrid3D < Vector5D > Fy("Fy");
	TimedGrid3D < Vector5D > Fz("Fz");

	U.setRangeT(0, T);
	U.setRangeX(0.5 * dx, Sx - 0.5 * dx);
	U.setRangeY(0.5 * dy, Sy - 0.5 * dy);
	U.setRangeZ(0.5 * dz, Sz - 0.5 * dz);

	U.setIndexRangeT(0, Nt);
	U.setIndexRangeX(0.5, Nx - 0.5);
	U.setIndexRangeY(0.5, Ny - 0.5);
	U.setIndexRangeZ(0.5, Nz - 0.5);
	U.build();

	Fx.setRangeT(0, T);
	Fx.setRangeX(0.0 * dx, Sx);
	Fx.setRangeY(0.5 * dy, Sy - 0.5 * dy);
	Fx.setRangeZ(0.5 * dz, Sz - 0.5 * dz);

	Fx.setIndexRangeT(0, Nt);
	Fx.setIndexRangeX(0, Nx);
	Fx.setIndexRangeY(0.5, Ny - 0.5);
	Fx.setIndexRangeZ(0.5, Nz - 0.5);
	Fx.build();
	Fx.fill(Vector5D(5));

	Fy.setRangeT(0, T);
	Fy.setRangeX(0.5 * dx, Sx - 0.5 * dx);
	Fy.setRangeY(0.0 * dy, Sy);
	Fy.setRangeZ(0.5 * dz, Sz - 0.5 * dz);

	Fy.setIndexRangeT(0, Nt);
	Fy.setIndexRangeX(0.5, Nx - 0.5);
	Fy.setIndexRangeY(0, Ny);
	Fy.setIndexRangeZ(0.5, Nz - 0.5);
	Fy.build();
	Fy.fill(Vector5D(5));

	Fz.setRangeT(0, T);
	Fz.setRangeX(0.5 * dx, Sx - 0.5 * dx);
	Fz.setRangeY(0.5 * dy, Sy - 0.5 * dy);
	Fz.setRangeZ(0.0, Sz);

	Fz.setIndexRangeT(0, Nt);
	Fz.setIndexRangeX(0.5, Nx - 0.5);
	Fz.setIndexRangeY(0.5, Ny - 0.5);
	Fz.setIndexRangeZ(0, Nz);
	Fz.build();
	Fz.fill(Vector5D(5));




	double cs = 3.0E8;//Скороcть cвета, m/s [cs]
	double csv = 377.323;//Сопротивление вакуума, sqrt(Gauss/Farad)=Ом [csv]
	double cze = 1.602E-19;//Заряд электрона, Кл. [cze]
	double cme = 9.10938215E-31;//Масса электрона, кг ? [cme]
	double ckb = 0.861E-4;//Константа Больцмана, э­В/Кград [ckb]
	double ci = 14.86;//Константа I в ф-ле для сечения ионизации, эВ [ci]
	double cih = 13.6;//Консттанта Ih в ф-ле для сечения ионизации, эВ [cih]
	double ca2 = 3.52E-16;//Константа 4*Pi*alfa2 в ф-ле для сечения ионизации, sm2 [ca2]
	double ce0 = 0.025;//Константа e0 в ф-ле для Qei, эВ [ce0]
	double cmb = 48.43E-27;//Масса иона воздуха,  29 масс ат. водорода кг ? [cmb]
	double cnw = 1.14e23;//Концентрация воздуха у поверхности Земли, 1/м3 ? [cnw]
	double cbv = -0.348E-4;//Показатель компоненты в распр. концентрации воздуха по высоте, ? [cbv]
	double cno2 = 0.21;//Доля кислорода в атмосфере, ? [cno2]
	double czi = 900000.0;//Высота начала ионосферы, м ? [czi]   60-90км
	double czb = 1.0;//Степень ионизации [czb]
	double cne0 = 28E6;//начальная концентрация электронов, 1/см3 ? [cne0]
	double cni0 = 28E6;//начальная концентрация ионов, 1/см3 ? [cni0]
	//!!! НЕ ИСПОЛЬЗУЕТСЯ! double ct0=200000.0;// ЖЕСТКОВАТО!!! 1000. 115070.   Температура воздуха у поверхности Земли, гр.К [ct0]
	//!!!double cti=200000.0;//  1000. 115070.   наачальная температура электронов, гр.К [cti]
	double chi = 4E9;//  8 5.5e8 !!!!!0.55e9 !!!3.e9          Опорная частота (Гц) [chi]
	double csi = 4E5;//  1.e8 1.5e8 !!!1.e6    Частота следований импульсов (скважность) (Гц) [csi]
	double cpi = 1.0;//  25 50 5  !!!!!3 1 !!!50.                Число импульсов в "пачке" [cpi]
	double cep = 1E7;//  !!!5.e5 !!!1.e7 1.e4       амплитуда напряженности электрического поля (В/м) [cep]
	double DTPR = 0.500500525E-12;//    2.5e-13 1.e-13 5.e-14   1. e4  Nachalnyj shag integrirowanija po vremeni DT (cek) [DTPR]
	double Z0 = 42000.0;//  40000. !!!50000. Nachalnaya vysota   176121372031.662 [Z0]
	double cedm = 1.761E11;//                 e/m [cedm]
	double cmdm = 1.878174685112534E-005;//   m/M [cmdm]

	double Wevlt = 1.602176487E-19;// Джоулей в одном э­лектрон-Вольте

	double k_b=1.38E-23;

	double NE0=28E6;
	double ATM_T=200.0; //температура атмоcферы K


	U[0].iterateWhole(GRID3D_ITERATOR{
	         U(0,ix,iy,iz) = HdVec3D::fromDensityPressureVelocity(NE0*cme , 3.0/2.0*k_b*ATM_T*NE0 , 0 ,0 ,0);
	});




	Timer mainTimer;
		for(double it=0;it<=Nt;it++){
			mainTimer.logTime("calculating layer "+ toString(it));


			Ex[it+1];
			Ey[it+1];
			Ez[it+1];
			Hx[it+1.5];
			Hy[it+1.5];
			Hz[it+1.5];


	        cout << "source\n" << endl;
			double t = it * DT;


	        //double sourceHz = H_AMPLITUDE*sin(2.0*M_PI*t/IMPULSE_TIME_WIDTH)*gaussStep(t,IMPULSE_TIME_WIDTH,IMPULSE_TIME_WIDTH);
            double sourceHz = H_AMPLITUDE*gauss(t,2*IMPULSE_TIME_WIDTH,IMPULSE_TIME_WIDTH);

            double sourceJz = J_AMPLITUDE*gauss(t,2*IMPULSE_TIME_WIDTH,IMPULSE_TIME_WIDTH);

	        //Hz(it+0.5, SOURCE_IX , SOURCE_IY , SOURCE_IZ) = sourceHz;
	        cout << "sourceHz=" << sourceHz << endl;



	        //double KZ=1;
	        //double KY=0.5;


	        /*
	        double KZ=0;
	        double KY=0;

	        double PAD_SIZE=0;

	        for(double iz=SOURCE_IZ-PAD_SIZE;iz<=SOURCE_IZ+PAD_SIZE;iz++)
	            for(double iy=SOURCE_IY-PAD_SIZE;iy<=SOURCE_IY+PAD_SIZE;iy++){

	                double DZ= (iz-SOURCE_IZ)*dz;
	                double DY= (iy-SOURCE_IY)*dy;

	                //double sourceHz=H_AMPLITUDE*sin(2.0*M_PI*  (t/IMPULSE_TIME_WIDTH - KY*DY - KZ*DZ)) *gaussStep(t,IMPULSE_TIME_WIDTH,IMPULSE_TIME_WIDTH);


	                Hz(it+0.5, SOURCE_IX , iy , iz) = sourceHz;
	            }
	        */







	        cout << "Ex" << endl;
	        Ex[it+1].iterateInternal(0,1,1,GRID3D_ITERATOR{
	            double EPS = 0.5*(epsilon(0,ix-0.5,iy,iz)+epsilon(0,ix+0.5,iy,iz));
	            double SIGMA=0.5*(sigma(0,ix-0.5,iy,iz)+sigma(0,ix+0.5,iy,iz));

	            Ex(it + 1, ix, iy, iz) =

	                (1.0-SIGMA*DT/2.0/EPS/EPS0)/(1.0+SIGMA*DT/2.0/EPS/EPS0)*Ex(it, ix, iy, iz)
	                + DT/(1.0+SIGMA*DT/2.0/EPS/EPS0)/EPS/EPS0*(
	                (Hz(it+0.5, ix, iy+0.5, iz) - Hz(it+0.5, ix, iy - 0.5, iz)) / dy
	                -
	                (Hy(it+0.5, ix, iy, iz+0.5) - Hy(it+0.5, ix, iy, iz - 0.5)) / dz
	                );

	        });

	        cout << "Ey" << endl;
	        Ey[it+1].iterateInternal(1,0,1,GRID3D_ITERATOR{
	            double EPS = 0.5*(epsilon(0,ix,iy-0.5,iz)+epsilon(0,ix,iy+0.5,iz));
	            double SIGMA=0.5*(sigma(0,ix,iy-0.5,iz)+sigma(0,ix,iy+0.5,iz));

	            Ey(it+1, ix, iy, iz)=
	                (1.0-SIGMA*DT/2.0/EPS/EPS0)/(1.0+SIGMA*DT/2.0/EPS/EPS0)*Ey(it, ix, iy, iz)
	                + DT/(1.0+SIGMA*DT/2.0/EPS/EPS0)/EPS/EPS0 * (
	                (Hx(it+0.5, ix, iy, iz+0.5) - Hx(it+0.5, ix, iy, iz-0.5)) / dz
	                -
	                (Hz(it+0.5, ix+0.5, iy, iz) - Hz(it+0.5, ix-0.5, iy, iz)) / dx
	                );
	        });

	        cout << "Ez" << endl;
	        Ez[it+1].iterateInternal(1,1,0,GRID3D_ITERATOR{
	            double EPS = 0.5*(epsilon(0,ix,iy,iz-0.5)+epsilon(0,ix,iy,iz+0.5));

	            double SIGMA=0.5*(sigma(0,ix,iy,iz-0.5)+sigma(0,ix,iy,iz+0.5));

	            double jz = 0;

	            if( (ix==J_SOURCE_IX) &&
	                (iy==J_SOURCE_IY) &&
	                (iz==J_SOURCE_IZ) ){

	                jz = sourceJz;

	                cout << "sourceJz: " << jz;
	            }

	            Ez(it + 1, ix, iy, iz)=
	                (1.0-SIGMA*DT/2.0/EPS/EPS0)/(1.0+SIGMA*DT/2.0/EPS/EPS0)*Ez(it, ix, iy, iz)
	                + DT/(1.0+SIGMA*DT/2.0/EPS/EPS0)/EPS/EPS0 * (
	                (Hy(it+0.5, ix+0.5, iy, iz) - Hy(it+0.5, ix-0.5, iy, iz)) / dx
	                -
	                (Hx(it+0.5, ix, iy+0.5, iz) - Hx(it+0.5, ix, iy-0.5, iz)) / dy

	                - jz
	                );

	        });


	        /*Liao 3rd order absorbing boundary conditions*/
	        //cout << "Ex1" << endl;
	        Ex[it+1].iterateBorderMinY(GRID3D_ITERATOR{
	            double u1=Ex(it,ix,iy+1,iz);
	            double u2=Ex(it-1,ix,iy+2,iz);
	            double u3=Ex(it-2,ix,iy+3,iz);
	            double u4=Ex(it-3,ix,iy+4,iz);

	            Ex(it+1,ix,iy,iz)=interpLiao3(u1,u2,u3,u4);
	        });
	        //cout << "Ex2" << endl;
	        Ex[it+1].iterateBorderMaxY(GRID3D_ITERATOR{
	            double u1=Ex(it,ix,iy-1,iz);
	            double u2=Ex(it-1,ix,iy-2,iz);
	            double u3=Ex(it-2,ix,iy-3,iz);
	            double u4=Ex(it-3,ix,iy-4,iz);

	            Ex(it+1,ix,iy,iz)=interpLiao3(u1,u2,u3,u4);
	        });
	        //cout << "Ex3" << endl;
	        Ex[it+1].iterateBorderMinZ(GRID3D_ITERATOR{
	            double u1=Ex(it,ix,iy,iz+1);
	            double u2=Ex(it-1,ix,iy,iz+2);
	            double u3=Ex(it-2,ix,iy,iz+3);
	            double u4=Ex(it-3,ix,iy,iz+4);

	            Ex(it+1,ix,iy,iz)=interpLiao3(u1,u2,u3,u4);
	        });
	        //cout << "Ex4" << endl;
	        Ex[it+1].iterateBorderMaxZ(GRID3D_ITERATOR{
	            double u1=Ex(it,ix,iy,iz-1);
	            double u2=Ex(it-1,ix,iy,iz-2);
	            double u3=Ex(it-2,ix,iy,iz-3);
	            double u4=Ex(it-3,ix,iy,iz-4);

	            Ex(it+1,ix,iy,iz)=interpLiao3(u1,u2,u3,u4);
	        });


	        //cout << "Ex5" << endl;
	        Ey[it+1].iterateBorderMinX(GRID3D_ITERATOR{
	            double u1=Ey(it,ix+1,iy,iz);
	            double u2=Ey(it-1,ix+2,iy,iz);
	            double u3=Ey(it-2,ix+3,iy,iz);
	            double u4=Ey(it-3,ix+4,iy,iz);

	            Ey(it+1,ix,iy,iz)=interpLiao3(u1,u2,u3,u4);
	        });

	        //cout << "Ex6" << endl;
	        Ey[it+1].iterateBorderMaxX(GRID3D_ITERATOR{
	            double u1=Ey(it,ix-1,iy,iz);
	            double u2=Ey(it-1,ix-2,iy,iz);
	            double u3=Ey(it-2,ix-3,iy,iz);
	            double u4=Ey(it-3,ix-4,iy,iz);

	            Ey(it+1,ix,iy,iz)=interpLiao3(u1,u2,u3,u4);
	        });

	        //cout << "Ex7" << endl;
	        Ey[it+1].iterateBorderMinZ(GRID3D_ITERATOR{
	            double u1=Ey(it,ix,iy,iz+1);
	            double u2=Ey(it-1,ix,iy,iz+2);
	            double u3=Ey(it-2,ix,iy,iz+3);
	            double u4=Ey(it-3,ix,iy,iz+4);

	            Ey(it+1,ix,iy,iz)=interpLiao3(u1,u2,u3,u4);
	        });

	        //cout << "Ex8" << endl;
	        Ey[it+1].iterateBorderMaxZ(GRID3D_ITERATOR{
	            double u1=Ey(it,ix,iy,iz-1);
	            double u2=Ey(it-1,ix,iy,iz-2);
	            double u3=Ey(it-2,ix,iy,iz-3);
	            double u4=Ey(it-3,ix,iy,iz-4);

	            Ey(it+1,ix,iy,iz)=interpLiao3(u1,u2,u3,u4);
	        });


	        //cout << "Ex9" << endl;
	        Ez[it+1].iterateBorderMinX(GRID3D_ITERATOR{
	            double u1=Ez(it,ix+1,iy,iz);
	            double u2=Ez(it-1,ix+2,iy,iz);
	            double u3=Ez(it-2,ix+3,iy,iz);
	            double u4=Ez(it-3,ix+4,iy,iz);

	            Ez(it+1,ix,iy,iz)=interpLiao3(u1,u2,u3,u4);
	        });
	        //cout << "Ex10" << endl;
	        Ez[it+1].iterateBorderMaxX(GRID3D_ITERATOR{
	            double u1=Ez(it,ix-1,iy,iz);
	            double u2=Ez(it-1,ix-2,iy,iz);
	            double u3=Ez(it-2,ix-3,iy,iz);
	            double u4=Ez(it-3,ix-4,iy,iz);

	            Ez(it+1,ix,iy,iz)=interpLiao3(u1,u2,u3,u4);
	        });


	        //cout << "Ex11" << endl;
	        Ez[it+1].iterateBorderMinY(GRID3D_ITERATOR{
	            double u1=Ez(it,ix,iy+1,iz);
	            double u2=Ez(it-1,ix,iy+2,iz);
	            double u3=Ez(it-2,ix,iy+3,iz);
	            double u4=Ez(it-3,ix,iy+4,iz);

	            Ez(it+1,ix,iy,iz)=interpLiao3(u1,u2,u3,u4);
	        });
	        //cout << "Ex12" << endl;
	        Ez[it+1].iterateBorderMaxY(GRID3D_ITERATOR{
	            double u1=Ez(it,ix,iy-1,iz);
	            double u2=Ez(it-1,ix,iy-2,iz);
	            double u3=Ez(it-2,ix,iy-3,iz);
	            double u4=Ez(it-3,ix,iy-4,iz);

	            Ez(it+1,ix,iy,iz)=interpLiao3(u1,u2,u3,u4);
	        });




	        cout << "Hx" << endl;
	        Hx[it+1.5].iterateWhole(GRID3D_ITERATOR{
	            double MU = 0.5*(mu(0,ix-0.5,iy,iz)+mu(0,ix+0.5,iy,iz));

	            Hx(it+1.5, ix, iy, iz)=
	                Hx(it+0.5, ix, iy, iz) -
	                DT/MU/MU0  * (
	                (Ez(it+1, ix, iy+0.5, iz) - Ez(it+1, ix, iy-0.5, iz)) / dy
	                -
	                (Ey(it+1, ix, iy, iz+0.5) - Ey(it+1, ix, iy, iz-0.5)) / dz
	                );
	        });

	        cout << "Hy" << endl;
	        Hy[it+1.5].iterateWhole(GRID3D_ITERATOR{
	            double MU = 0.5*(mu(0,ix,iy-0.5,iz)+mu(0,ix,iy+0.5,iz));

	            Hy(it+1.5, ix, iy, iz)=
	                Hy(it+0.5, ix, iy, iz) -
	                DT/MU/MU0  * (
	                (Ex(it+1, ix, iy, iz+0.5) - Ex(it+1, ix, iy, iz-0.5)) / dz
	                -
	                (Ez(it+1, ix+0.5, iy, iz) - Ez(it+1, ix-0.5, iy, iz)) / dx
	                );
	        });


	        cout << "Hz" << endl;
	        Hz[it+1.5].iterateWhole(GRID3D_ITERATOR{
	            double MU = 0.5*(mu(0,ix,iy,iz-0.5)+mu(0,ix,iy,iz+0.5));

	            Hz(it+1.5, ix, iy, iz)=
	                Hz(it+0.5, ix, iy, iz) -
	                DT/MU/MU0  * (
	                (Ey(it+1, ix+0.5, iy, iz) - Ey(it+1, ix-0.5, iy, iz)) / dx
	                -
	                (Ex(it+1, ix, iy+0.5, iz) - Ex(it+1, ix, iy-0.5, iz)) / dy
	                );
	        });

#ifdef RHS_ENABLE


	       Fx[it].iterateInternal(1,0,0,GRID3D_ITERATOR {
					Fx(it,ix,iy,iz)=HdFlowVec3D::X::riemannFlux( U(it,ix-0.5,iy,iz) , U(it,ix+0.5,iy,iz) );
				});

		   Fy[it].iterateInternal(0,1,0,GRID3D_ITERATOR {
					Fy(it,ix,iy,iz)=HdFlowVec2D::Y::riemannFlux( U(it,ix,iy-0.5,iz) , U(it,ix,iy+0.5,iz) );
				});

		   Fz[it].iterateInternal(0,0,1,GRID3D_ITERATOR {
		   			Fz(it,ix,iy,iz)=HdFlowVec3D::Z::riemannFlux( U(it,ix,iy,iz-0.5) , U(it,ix,iy,iz+0.5) );
		   });

		//transparent BC
		Fx[it].iterateBorderMinX(GRID3D_ITERATOR {
					Fx(it,ix,iy,iz)=Fx(it,ix+1,iy,iz);
		});

		Fx[it].iterateBorderMaxX(GRID3D_ITERATOR {
					Fx(it,ix,iy,iz)=Fx(it,ix-1,iy,iz);
		});

		Fy[it].iterateBorderMinY(GRID3D_ITERATOR {
					Fy(it,ix,iy,iz)=Fy(it,ix,iy+1,iz);
		});

		Fy[it].iterateBorderMaxY(GRID3D_ITERATOR {
					Fy(it,ix,iy,iz)=Fy(it,ix,iy-1,iz);
		});

		Fz[it].iterateBorderMinZ(GRID3D_ITERATOR {
					Fz(it,ix,iy,iz)=Fz(it,ix,iy,iz+1);
		});

		Fz[it].iterateBorderMaxZ(GRID3D_ITERATOR {
					Fz(it,ix,iy,iz)=Fz(it,ix,iy,iz-1);
		});

		U.iterateWhole(GRID3D_ITERATOR{

		            //DISSIPATION STEP
					double ne = HdVec3D::density(U(it,ix,iy,iz)) / cme;

					double Wvx=HdVec3D::velocityX(U(it,ix,iy,iz));//TODO или тут импульс???
					double Wvy=HdVec3D::velocityY(U(it,ix,iy,iz));
					double Wvz=HdVec3D::velocityZ(U(it,ix,iy,iz));


					double WEt= (U(it,ix,iy,iz).c4 - U(it+1,ix,iy,iz).c0*sqr(Wvx,Wvy,Wvz)/2.0)/ne/Wevlt;//1.0;
					WEtGrid(ix,iy,iz)=WEt;

		});

#endif



	    if (isEveryNth(it, 5)) {

			cimgSaver.saveSliceZ(Hz[it + 0.5], Nz / 2, frame("Hz_", it, "png"),GRID3D_CALCULATOR{
				return Hz(it + 0.5,ix,iy,iz);
			});

			//vtiSaver.save(Hz[it+0.5],frame("Hz_",it,"vti"));


			vtiSaver.save(Ex[it],frame("Ex_",it,"vti"));
			vtiSaver.save(Ey[it],frame("Ey_",it,"vti"));
			vtiSaver.save(Ez[it],frame("Ez_",it,"vti"));

		}


		}
		mainTimer.logTime("calculation finished");



	return 0;
}
