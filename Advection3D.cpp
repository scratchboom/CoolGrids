
#include "GridsCommon.hpp"

int main(){

	VtiSaver3D vtiSaver;

	const double DT=0.1;

	double hx,hy,hz;
	hx=hy=hz=0.1;

	const double Ax=0.5;
	const double Ay=0.3;
	const double Az=0.0;

	const int Nx=256;
	const int Ny=256;
	const int Nz=8;


	const int SAVE_EVERY_N=1;

    TimedGrid3D<double> U("Advection Value");

	U.setRangeT(0,1);U.setIndexRangeT(0,100);

	U.setRangeX(0,Nx*hx);U.setIndexRangeX(0,Nx);
	U.setRangeY(0,Ny*hy);U.setIndexRangeY(0,Ny);
	U.setRangeZ(0,Nz*hz);U.setIndexRangeZ(0,Nz);

	U.build();

	double it=U.minIndexT;
	U[it].iterateWhole(GRID3D_ITERATOR{
		U(it,ix,iy,iz)=gauss( length((Nx/2-ix)*hx,(Ny/2-iy)*hy,(Nz/2-iz)*hz) ,0,2);
	});


	for(double it=U.minIndexT;it<U.maxIndexT;it++){

		U[it].iterateBorder(GRID3D_ITERATOR{
				U(it,ix,iy,iz)=0.0;
		});


		U[it].iterateInternal(GRID3D_ITERATOR{
			U(it+1,ix,iy,iz) = U(it,ix,iy,iz)
					           + Ax*DT/hx*( U(it,ix-1,iy,iz)-U(it,ix,iy,iz) )
		                       + Ay*DT/hy*( U(it,ix,iy-1,iz)-U(it,ix,iy,iz) )
		                       + Az*DT/hz*( U(it,ix,iy,iz-1)-U(it,ix,iy,iz) );
		});

		if((int)it%SAVE_EVERY_N==0){
			vtiSaver.save(U[it],"U_"+toZPadString(it,5)+".vti");
		}

	}


	return 0;
}
