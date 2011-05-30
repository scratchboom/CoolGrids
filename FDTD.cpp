
#include "GridsCommon.hpp"

using namespace std;

int main(){

	omp_set_num_threads(2);

	VtiSaver3D vtiSaver;
	CImgSaver2D cimgSaver;


    double EPS0=8.854E-12;// Ф/м
    double MU0=1.256E-6;// Гн/м
    double C=299792458.0;// м/с    C=sqrt(1/EPS0/MU0)
    double VACUUM_WAVE_RESISTIVITY=sqrt(MU0/EPS0);// волновое сопротивление вакуума [Ом]
    double cze=1.602E-19;//               Заряд электрона, Кл. [cze]
    double cme=9.10938215E-31;//               Масса электрона, кг ? [cme]


    double Nx=256;
    double Ny=257;
    double Nz=20;
    double Nt=2000;

    double SOURCE_IX = (int)(Nx*0.1)+0.5;
    double SOURCE_IY = (int)(Ny*0.5)+0.5;
    double SOURCE_IZ = Nz/2;

    double IMPULSE_FREQ = 3e8;
    double IMPULSE_TIME = 1.0/IMPULSE_FREQ;
    double IMPULSE_LENGTH = C*IMPULSE_TIME;

    double E_AMPLITUDE=1.0;
    double H_AMPLITUDE=E_AMPLITUDE/VACUUM_WAVE_RESISTIVITY;

	cimgSaver.setValueRange(-2.0*H_AMPLITUDE , 2.0*H_AMPLITUDE);

	/*
    double Sx=10.0;
	double Sy= (Sx/Nx) * Ny;
	double Sz= (Sx/Nx) * Nz;
	double T=1.0;

	double dx=Sx/Nx;
	double dy=Sy/Ny;
	double dz=Sz/Nz;
	*/
	//double DT=T/Nt;




	double dx,dy,dz;
	dx=dy=dz = IMPULSE_LENGTH/32.0;

	double Sx= dx*Nx;
	double Sy= dy*Ny;
	double Sz= dz*Nz;
	double T=1.0;

    double DT = 0.5*dx/C;


    std::cout << "stability parameters:" << std::endl;

    CHECK(C*DT/dx,<,1/sqrt(3));
    CHECK(IMPULSE_TIME/DT,>=,32);
    CHECK(IMPULSE_LENGTH/dx,>=,32);


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

	Ex.setLayersCountToMaintain(5);
	Ex.build();


	Ey.setRangeT(0,T);
	Ey.setRangeX(0 , Sx);
	Ey.setRangeY(0.5*dy , Sy-0.5*dy);
	Ey.setRangeZ(0 , Sz);

	Ey.setIndexRangeT(0,Nt);
	Ey.setIndexRangeX(0,Nx);
	Ey.setIndexRangeY(0.5,Ny-0.5);
    Ey.setIndexRangeZ(0,Nz);

    Ey.setLayersCountToMaintain(5);
	Ey.build();


	Ez.setRangeT(0,T);
	Ez.setRangeX(0 , Sx);
	Ez.setRangeY(0 , Sy);
	Ez.setRangeZ(0.5*dz , Sz-0.5*dz);

	Ez.setIndexRangeT(0,Nt);
	Ez.setIndexRangeX(0 , Nx);
	Ez.setIndexRangeY(0 , Ny);
	Ez.setIndexRangeZ(0.5 , Nz-0.5);

	Ez.setLayersCountToMaintain(5);
	Ez.build();



	Hx.setRangeT(0,T);
	Hx.setRangeX(0 , Sx);
	Hx.setRangeY(0.5*dy , Sy-0.5*dy);
	Hx.setRangeZ(0.5*dz , Sz-0.5*dz);

	Hx.setIndexRangeT(0.5,Nt+0.5);
	Hx.setIndexRangeX(0   , Nx);
	Hx.setIndexRangeY(0.5 , Ny-0.5);
	Hx.setIndexRangeZ(0.5 , Nz-0.5);

	Hx.setLayersCountToMaintain(2);
	Hx.build();




	Hy.setRangeT(0,T);
	Hy.setRangeX(0.5*dx , Sx-0.5*dx);
	Hy.setRangeY(0      , Sy       );
	Hy.setRangeZ(0.5*dz , Sz-0.5*dz);

	Hy.setIndexRangeT(0.5,Nt+0.5);
	Hy.setIndexRangeX(0.5 , Nx-0.5);
	Hy.setIndexRangeY(0   , Ny    );
    Hy.setIndexRangeZ(0.5 , Nz-0.5);

    Hy.setLayersCountToMaintain(2);
	Hy.build();




	Hz.setRangeT(0,T);
	Hz.setRangeX(0.5*dx , Sx-0.5*dx);
	Hz.setRangeY(0.5*dy , Sy-0.5*dy);
	Hz.setRangeZ(0      , Sz);

	Hz.setIndexRangeT(0.5,Nt+0.5);
	Hz.setIndexRangeX(0.5 , Nx-0.5);
	Hz.setIndexRangeY(0.5 , Ny-0.5);
	Hz.setIndexRangeZ(0   , Nz    );

	Hz.setLayersCountToMaintain(2);
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


    #define PRYSM_
    #ifdef PRYSM
	epsilon[0].fillAntialiased(1.5,GRID3D_CONDITION{
		return (ix+iy > 1.2*avg(Nx,Ny)) && (ix < 0.8*Nx) && (iy < 0.8*Ny);
	});

	mu[0].fillAntialiased(1.5,GRID3D_CONDITION{
		return (ix+iy > 1.2*avg(Nx,Ny)) && (ix < 0.8*Nx) && (iy < 0.8*Ny);
	});
    #endif




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

	vtiSaver.save(epsilon[0],"epsilon.vti");
	vtiSaver.save(mu[0],"mu.vti");


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


	        double sourceHz = H_AMPLITUDE*sin(2.0*M_PI*t/IMPULSE_TIME)*gaussStep(t,IMPULSE_TIME*4,IMPULSE_TIME);


	        //Hz(it+0.5, SOURCE_IX , SOURCE_IY , SOURCE_IZ) = sourceHz;
	        cout << "sourceHz=" << sourceHz << endl;


	        //double KZ=1;
	        //double KY=0.5;


	        double KZ = 0;
	        double KY = 0;

	        int PAD_SIZE_Y = 100;
	        int PAD_SIZE_Z = 5;

	        for(double iz=SOURCE_IZ-PAD_SIZE_Z;iz<=SOURCE_IZ+PAD_SIZE_Z;iz++)
	            for(double iy=SOURCE_IY-PAD_SIZE_Y;iy<=SOURCE_IY+PAD_SIZE_Y;iy++){

	                double DZ= (iz-SOURCE_IZ)*dz;
	                double DY= (iy-SOURCE_IY)*dy;

	                double DR = length(DZ,DY);

	                //double sourceHz=H_AMPLITUDE*sin(2.0*M_PI*  (t/IMPULSE_TIME_WIDTH - KY*DY - KZ*DZ)) *gaussStep(t,IMPULSE_TIME_WIDTH*4,IMPULSE_TIME_WIDTH);
	                sourceHz=H_AMPLITUDE * onePlusCosPulse(DR/(dx*length(PAD_SIZE_Y,PAD_SIZE_Z)) *M_PI) *sin(2.0*M_PI*  (t/IMPULSE_TIME - KY*DY - KZ*DZ)) *gaussStep(t,IMPULSE_TIME*1.5,IMPULSE_TIME);
	                //sourceHz=H_AMPLITUDE *sin(2.0*M_PI*  (t/IMPULSE_TIME - KY*DY - KZ*DZ)) *gaussStep(t,IMPULSE_TIME*1.5,IMPULSE_TIME);

	                //sourceHz=H_AMPLITUDE * cosPulse(DR/(dx*length(PAD_SIZE_Y,PAD_SIZE_Z)) * 0.5*M_PI) *sin(2.0*M_PI*  (t/IMPULSE_TIME - KY*DY - KZ*DZ)) *gaussStep(t,IMPULSE_TIME*1.5,IMPULSE_TIME);


	                Hz(it+0.5, SOURCE_IX , iy , iz) = sourceHz;
	            }







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

	            Ez(it + 1, ix, iy, iz)=
	                (1.0-SIGMA*DT/2.0/EPS/EPS0)/(1.0+SIGMA*DT/2.0/EPS/EPS0)*Ez(it, ix, iy, iz)
	                + DT/(1.0+SIGMA*DT/2.0/EPS/EPS0)/EPS/EPS0 * (
	                (Hy(it+0.5, ix+0.5, iy, iz) - Hy(it+0.5, ix-0.5, iy, iz)) / dx
	                -
	                (Hx(it+0.5, ix, iy+0.5, iz) - Hx(it+0.5, ix, iy-0.5, iz)) / dy
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


	    if (isEveryNth(it, 5)) {
			cimgSaver.saveSliceZ(Hz[it + 0.5], Nz / 2, frame("Hz_", it, "png"),GRID3D_CALCULATOR{
				return Hz(it + 0.5,ix,iy,iz);
			});
		}

	    if (isEveryNth(it, 50)) {
			vtiSaver.save(Hz[it+0.5],frame("Hz_",it,"vti"));
		}


		}
		mainTimer.logTime("calculation finished");



	return 0;
}
