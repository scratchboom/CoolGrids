#include "GridsCommon.hpp"

int main(){

	GnuPlotSaver1D gnuPlotSaver;
	gnuPlotSaver.setLineColor("#FF0000")
			    .setValueRange(0.0,1.2)
			    .setSaveToPNG(true);


	double N=512;
    double NT=10000;

	double S=1.0;
	double dx=S/N;

	double SGM=0.5;
	double MAX_LAMBDA=5.0;
	double dt=dx/MAX_LAMBDA * SGM;


	//Lax-Friedrichs method
	TimedGrid1D<Vec3D> U("Hydrodynamics vector");
	TimedGrid1D<Vec3D> F("Flow vector");

	U.setRangeX(0.5*dx,N*dx-0.5*dx);
	U.setIndexRangeX(0.5,N-0.5);
	U.build();

	F.setRangeX(0,dx*N);
	F.setIndexRangeX(0,N);
	F.build();


	U[0].clear();
	F[0].clear();

	U[0].iterateWhole(GRID1D_ITERATOR{
		if(ix<N/2) U(0,ix) = HydroDynVec::from_Density_Pressure_Velocity(1,1,0);
		else U(0,ix) = HydroDynVec::from_Density_Pressure_Velocity(0.125,0.1,0);
	});

	for(double it=0;it<NT;it++){
		DBGVAL(it);

		U(it,U.minIndexX)=U(it,U.minIndexX+1);
	    U(it,U.maxIndexX)=U(it,U.maxIndexX-1);

		F[it].iterateInternal(GRID1D_ITERATOR{

			double v1=HydroDynVec::velocityX( U(it,ix-0.5) );
			double rho1=HydroDynVec::density( U(it,ix-0.5) );
			double p1=HydroDynVec::pressure( U(it,ix-0.5) );

			double c1=soundVelocity(GAMMA,p1,rho1);
			double lmb11=v1;
			double lmb12=v1+c1;
			double lmb13=v1-c1;


			double v2=HydroDynVec::velocityX( U(it,ix+0.5) );
			double rho2=HydroDynVec::density( U(it,ix+0.5) );
			double p2=HydroDynVec::pressure( U(it,ix+0.5) );

			double c2=soundVelocity(GAMMA,p2,rho2);
			double lmb21=v2;
			double lmb22=v2+c2;
			double lmb23=v2-c2;

			double maxAbsLambda=maxAbs(lmb11,lmb12,lmb13,
			                           lmb21,lmb22,lmb23);

			//DBGVAL(toString(U(it,ix-0.5)));
			//DBGVAL(toString(U(it,ix+0.5)));

			Vec3D F1=HydroDynVec::toFlow( U(it,ix-0.5) );
			Vec3D F2=HydroDynVec::toFlow( U(it,ix+0.5) );

			//DBGVAL(toString(F1));
			//DBGVAL(toString(F2));

			F(it,ix)  =  (F1+F2)/2.0 - 0.5*maxAbsLambda*( U(it,ix+0.5) - U(it,ix-0.5));

			//DBGVAL(toString(F(it,ix)));
		});


		U[it+1].iterateWhole(GRID1D_ITERATOR{
			U(it+1,ix)=U(it,ix) - dt/dx*(F(it,ix+0.5)-F(it,ix-0.5));
			//DBGVAL( toString(U(it+1,ix)) );
		});


		if(isEveryNth(it,10)){

			gnuPlotSaver.save(U[it],frame("density_",it,"plot"),GRID1D_CALCULATOR{
				return HydroDynVec::density(U(it,ix));
			});

//			gnuPlotSaver.save(U[it],frame("pressure_",it,"plot"),GRID1D_CALCULATOR{
//				return HydroDynVec::pressure(U(it,ix));
//			});
//
//			gnuPlotSaver.save(U[it],frame("velocityX_",it,"plot"),GRID1D_CALCULATOR{
//				return HydroDynVec::velocityX(U(it,ix));
//			});
//
//			gnuPlotSaver.save(U[it],frame("Energy_",it,"plot"),GRID1D_CALCULATOR{
//				return HydroDynVec::fullEnergyPerVolumeUnit(U(it,ix));
//			});
//
//			gnuPlotSaver.save(F[it],frame("Flow_",it,"plot"),GRID1D_CALCULATOR{
//				return F(it,ix)[0];
//			});

		}

	}
	return 0;
}
