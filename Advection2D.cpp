
#include "GridsCommon.hpp"

int main(){



	BmpSaver bmpSaver;
	bmpSaver.setValueRange(0.0,1.1);

	//угловая скорость вращения вихря
	double W0=1;

	//координаты центра вихря
	double X0=3;
	double Y0=3;


	TimedGrid2D<double> U("Advection Value");

	U.setRangeT(0,1);U.setIndexRangeT(0,1000);

	U.setRangeX(0,10);U.setIndexRangeX(0,200);
	U.setRangeY(0,10);U.setIndexRangeY(0,200);

	U.build();

	double hx=U[0].dx;
    double hy=U[0].dy;
    double DT=hx/20.0;

	double it=U.minIndexT;
	U[it].iterateWhole(GRID2D_ITERATOR{
		U(it,ix,iy)=gauss( length(hx*ix-5,hy*iy-5),0,1.7);
	});




	Timer mainTimer;
	for(double it=U.minIndexT;it<U.maxIndexT;it++){
		Accumulator<double> lambdaAcc;
		mainTimer.logTime(toString(it));

		U[it].iterateBorder(GRID2D_ITERATOR{
				U(it,ix,iy)=0.0;
		});

		U[it].iterateInternal(GRID2D_ITERATOR{

			double x=ix*hx;
			double y=iy*hy;

			double Vx = -W0*(y-Y0);
			double Vy = W0*(x-X0);

			double lamXp=max(Vx,0);
			double lamXm=min(Vx,0);

			double lamYp=max(Vy,0);
			double lamYm=min(Vy,0);

			lambdaAcc.add(Vx);
			lambdaAcc.add(Vy);


			U(it+1,ix,iy) = U(it,ix,iy)
					        + lamXp*DT/hx*( U(it,ix-1,iy)-U(it,ix,iy) )
					        + lamXm*DT/hx*( U(it,ix,iy)-U(it,ix+1,iy) )

		                    + lamYp*DT/hy*( U(it,ix,iy-1)-U(it,ix,iy)
		                    + lamYm*DT/hy*( U(it,ix,iy)-U(it,ix,iy+1) ));
		});

		DBGVAL(lambdaAcc.getMaxAbs());

		if(isEveryNth(it,10)){
			bmpSaver.save(U[it],"./plot/vortex_"+toZPadString(it,5)+".bmp");
		}
	}


	return 0;
}
