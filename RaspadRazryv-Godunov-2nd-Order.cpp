#include "GridsCommon.hpp"


int main(){

	GnuPlotSaver1D gnuPlotSaver;
	gnuPlotSaver.setLineColor("#FF0000")
			    .setValueRange(0.0,12.0)
			    .setSaveToPNG(true);


	double N=512;
    double NT=10000;

	double S=1.0;
	double dx=S/N;

	double SGM=0.5;
	double MAX_LAMBDA=5.0;
	double dt=dx/MAX_LAMBDA * SGM;

	soundVelocity(GAMMA,10,10);

	//Lax-Friedrichs method
	TimedGrid1D<Vec3D> U("Hydrodynamics vector");
	TimedGrid1D<Vec3D> Q("Slope of Hydrodynamics vector");
	TimedGrid1D<Vec3D> F("Flow vector");

	U.setRangeX(0.5*dx,N*dx-0.5*dx);
	U.setIndexRangeX(0.5,N-0.5);
	U.build();

	Q.setRangeX(0.5*dx,N*dx-0.5*dx);
    Q.setIndexRangeX(0.5,N-0.5);
	Q.build();

	F.setRangeX(0,dx*N);
	F.setIndexRangeX(0,N);
	F.build();


	U[0].clear();
	F[0].clear();

	U[0].iterateWhole(GRID1D_ITERATOR{
		if(ix<N/2) U(0,ix) = HydroDynVec::from_Density_Pressure_Velocity(10,10,0);
		else U(0,ix) = HydroDynVec::from_Density_Pressure_Velocity(1,1,0);
	});

	Q[0].iterateInternal(GRID1D_ITERATOR{
		Q(0,ix)[0] = minMod( (U(0,ix+1)[0]-U(0,ix)[0])/dx , (U(0,ix)[0]-U(0,ix-1)[0])/dx );
		Q(0,ix)[1] = minMod( (U(0,ix+1)[1]-U(0,ix)[1])/dx , (U(0,ix)[1]-U(0,ix-1)[1])/dx );
		Q(0,ix)[2] = minMod( (U(0,ix+1)[2]-U(0,ix)[2])/dx , (U(0,ix)[2]-U(0,ix-1)[2])/dx );
	});

	for(double it=0;it<NT;it++){
		DBGVAL(it);

		Q[it].iterateInternal(GRID1D_ITERATOR{
				Q(it,ix)[0] = minMod( (U(it,ix+1)[0]-U(it,ix)[0])/dx , (U(it,ix)[0]-U(it,ix-1)[0])/dx );
				Q(it,ix)[1] = minMod( (U(it,ix+1)[1]-U(it,ix)[1])/dx , (U(it,ix)[1]-U(it,ix-1)[1])/dx );
				Q(it,ix)[2] = minMod( (U(it,ix+1)[2]-U(it,ix)[2])/dx , (U(it,ix)[2]-U(it,ix-1)[2])/dx );
		});

		Q(it,Q.minIndexX)[0] = 0;
	    Q(it,Q.minIndexX)[1] = 0;
		Q(it,Q.minIndexX)[2] = 0;

		Q(it,Q.maxIndexX)[0] = 0;
	    Q(it,Q.maxIndexX)[1] = 0;
		Q(it,Q.maxIndexX)[2] = 0;


		F[it].iterateInternal(GRID1D_ITERATOR{

			double pressure,density,energy,velocity;

			double rho1=HydroDynVec::density(U(it,ix-0.5));
			double p1=HydroDynVec::pressure(U(it,ix-0.5));
			double s1=soundVelocity(GAMMA,p1,rho1);

			double rho2=HydroDynVec::density(U(it,ix+0.5));
			double p2=HydroDynVec::pressure(U(it,ix+0.5));
			double s2=soundVelocity(GAMMA,p2,rho2);

			Real4 result=SolveRiemannProblem(HydroDynVec::density(U(it,ix-0.5)),
							HydroDynVec::internalEnergyPerMassUnit(U(it,ix-0.5)),
							HydroDynVec::velocityX(U(it,ix-0.5)),
							GAMMA,s1,
							HydroDynVec::density(U(it,ix+0.5)),
							HydroDynVec::internalEnergyPerMassUnit(U(it,ix+0.5)),
							HydroDynVec::velocityX(U(it,ix+0.5)),
							GAMMA,s2);
            pressure=Pressure(result);
            density=Density(result);
            velocity=Velocity(result);
            energy=Energy(result);


			F(it,ix)=FlowVec::from_Density_Pressure_Velocity(density,pressure,velocity);
			//DBGVAL(toString(F(it,ix)));
		});

		F(it,F.minIndexX)=F(it,F.minIndexX+1);
	    F(it,F.maxIndexX)=F(it,F.maxIndexX-1);

		U[it+1].iterateWhole(GRID1D_ITERATOR{

			Vec3D U_ = U(it,ix) - dt/dx*(HydroDynVec::toFlow(U(it,ix)+Q(it,ix)*dx/2.0) - HydroDynVec::toFlow(U(it,ix)-Q(it,ix)*dx/2.0));
			Vec3D U__ = 0.5*(U(it,ix) + U_);

			Fp

			U(it+1,ix)=U(it,ix) - dt/dx*(F(it,ix+0.5)-F(it,ix-0.5));
			//DBGVAL( toString(U(it+1,ix)) );
		});


		if(isEveryNth(it,50)){

			gnuPlotSaver.save(U[it],frame("density_",it,"plot"),GRID1D_CALCULATOR{
				return HydroDynVec::density(U(it,ix));
			});

			gnuPlotSaver.save(U[it],frame("pressure_",it,"plot"),GRID1D_CALCULATOR{
				return HydroDynVec::pressure(U(it,ix));
			});

			gnuPlotSaver.save(U[it],frame("velocityX_",it,"plot"),GRID1D_CALCULATOR{
				return HydroDynVec::velocityX(U(it,ix));
			});

			gnuPlotSaver.save(U[it],frame("Energy_",it,"plot"),GRID1D_CALCULATOR{
				return HydroDynVec::fullEnergyPerVolumeUnit(U(it,ix));
			});

			gnuPlotSaver.save(F[it],frame("Flow_",it,"plot"),GRID1D_CALCULATOR{
				return F(it,ix)[0];
			});

		}

	}
	return 0;
}

