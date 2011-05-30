
#include "GridsCommon.hpp"

using namespace std;

int main(){

	DBGVAL(omp_get_num_procs());
	DBGVAL(omp_get_num_threads());

	omp_set_num_threads(2);

	VtiSaver3D vtiSaver;

    double Nx=128;
    double Ny=81;
    double Nz=80;
    double Nt=20;

    double Sx=10.0;
	double Sy= (Sx/Nx) * Ny;
	double Sz= (Sx/Nx) * Nz;
	double T=1.0;


    double dx = Sx / Nx;
	double dy = Sy / Ny;
	double dz = Sz / Nz;




	TimedGrid3D<double> Ex("Ex");
	TimedGrid3D<double> Hx("Hx");

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

	Ex[0].iterateWhole(GRID3D_ITERATOR{
		    Ex(0,ix,iy,iz) = 0;
	});

	Ex[1].iterateWhole(GRID3D_ITERATOR{
			    Ex(1,ix,iy,iz) = 1;
	});

	Hx[0.5].iterateWhole(GRID3D_ITERATOR{
		    Hx(0.5,ix,iy,iz) = 0;
	});

	Hx[1.5].iterateWhole(GRID3D_ITERATOR{
		Hx(1.5,ix,iy,iz) = 1;
	});





	Timer mainTimer;
		for(double it=1;it<=Nt;it++){
			mainTimer.logTime("calculating layer "+ toString(it));


			Ex[it+1];
			Hx[it+1.5];

			Ex[it+1].iterateWhole(GRID3D_ITERATOR{
				    Ex(it+1,ix,iy,iz) = Ex(it,ix,iy,iz)+Ex(it-1,ix,iy,iz);
			});

			Hx[it+1.5].iterateWhole(GRID3D_ITERATOR{
				    Hx(it+1.5,ix,iy,iz) = Hx(it+0.5,ix,iy,iz)+Hx(it-0.5,ix,iy,iz);
			});

			vtiSaver.save(Ex[it],frame("Ex_",it,"vti"));
			vtiSaver.save(Hx[it+0.5],frame("Hx_",it,"vti"));

		}
		mainTimer.logTime("calculation finished");



	return 0;
}
