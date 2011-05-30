
#include "GridsCommon.hpp"

int main(){

	GnuPlotSaver1D gnuPlotSaver;
	gnuPlotSaver.setLineColor("#FF00FF")
			    .setValueRange(0.0,1.1)
			    .setSaveToPNG(true);

	const double DT=0.1;
	const double h=0.1;
	const double A=0.5;

	const double LEFT_VALUE=1.0;
	const double RIGHT_VALUE=0.1;

	const int SAVE_EVERY_N=1;

	DBGLNVAL(A*DT/h);

	TimedGrid1D<double> U("Advection Value");

	U.setRangeT(0,1);U.setIndexRangeT(0,100);
	U.setRangeX(0,10);U.setIndexRangeX(0,200);

	U.build();

	//initial conditions
	double it=U.minIndexT;

	U[it].iterateWhole(GRID1D_ITERATOR{
		if(ix<U.maxIndexX/4)U(it,ix)=LEFT_VALUE;
		else U(it,ix)=RIGHT_VALUE;
	});

	for(double it=U.minIndexT;it<U.maxIndexT;it++){

		U(it,U.minIndexX)=LEFT_VALUE;
		U(it,U.maxIndexX)=RIGHT_VALUE;

		U[it].iterateInternal(GRID1D_ITERATOR{
			U(it+1,ix) = U(it,ix) + A*DT/h*( U(it,ix-1)-U(it,ix) );
		});

		if((int)it%SAVE_EVERY_N==0){
			gnuPlotSaver.save(U[it],"./plot/U_"+toZPadString(it,5)+".plot");
		}
	}

	return 0;
}
