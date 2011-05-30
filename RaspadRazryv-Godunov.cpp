#include "GridsCommon.hpp"

int main(){

	GnuPlotSaver1D gnuPlotSaver;
	gnuPlotSaver.setLineColor("#FF0000")
			    .setValueRange(0.0,5.0)
			    .setSaveToPNG(true);


	double N=512;
    double NT=50000;

	double S=1.0;
	double dx=S/N;

	double SGM=0.5;
	double MAX_LAMBDA=5.0;
	double dt=dx/MAX_LAMBDA * SGM;

	soundVelocity(GAMMA,10,10);

	TimedGrid1D<Vector3D> U("Hydrodynamics vector");
	TimedGrid1D<Vector3D> F("Flow vector");

	U.setRangeX(0.5*dx,N*dx-0.5*dx);
	U.setIndexRangeX(0.5,N-0.5);
	U.build();

	F.setRangeX(0,dx*N);
	F.setIndexRangeX(0,N);
	F.build();


	U.fill(Vector3D(3));
	F.fill(Vector3D(3));

	U[0].iterateWhole(GRID1D_ITERATOR{
		if(ix<N/2) U(0,ix) = HdVec1D::fromDensityPressureVelocity(2,2,0);
		else U(0,ix) = HdVec1D::fromDensityPressureVelocity(1,1,0);
	});

	for(double it=0;it<NT;it++){
		DBGVAL(it);

		F[it].iterateInternal(GRID1D_ITERATOR{

			double pressure,density,energy,velocity;

			double rho1=HdVec1D::density(U(it,ix-0.5));
			double p1=HdVec1D::pressure(U(it,ix-0.5));
			double s1=soundVelocity(GAMMA,p1,rho1);

			double rho2=HdVec1D::density(U(it,ix+0.5));
			double p2=HdVec1D::pressure(U(it,ix+0.5));
			double s2=soundVelocity(GAMMA,p2,rho2);

			Real4 result=SolveRiemannProblem(HdVec1D::density(U(it,ix-0.5)),
					HdVec1D::internalEnergyPerMassUnit(U(it,ix-0.5)),
					HdVec1D::velocityX(U(it,ix-0.5)),
							GAMMA,s1,
							HdVec1D::density(U(it,ix+0.5)),
							HdVec1D::internalEnergyPerMassUnit(U(it,ix+0.5)),
							HdVec1D::velocityX(U(it,ix+0.5)),
							GAMMA,s2);
            pressure=Pressure(result);
            density=Density(result);
            velocity=Velocity(result);
            energy=Energy(result);


			F(it,ix)=HdFlowVec1D::X::fromDensityPressureVelocity(density,pressure,velocity);

		});

		//transparent
		/*
		F(it,F.minIndexX)=F(it,F.minIndexX+1);
	    F(it,F.maxIndexX)=F(it,F.maxIndexX-1);
	    */

		//reflexive
		{
			Vector3D u = U(it,  0.5);

			double rho = HdVec1D::density(u);
			double p = HdVec1D::pressure(u);
			double v = HdVec1D::velocityX(u);

			Vector3D uw = HdVec1D::fromDensityPressureVelocity(rho, p, -v);

			F(it, 0) = HdFlowVec1D::X::riemannFlux(uw, u);
		}
		{
			Vector3D u = U(it, N - 0.5);

			double rho = HdVec1D::density(u);
			double p = HdVec1D::pressure(u);
			double v = HdVec1D::velocityX(u);

			Vector3D uw = HdVec1D::fromDensityPressureVelocity(rho, p, -v);

			F(it, N) = HdFlowVec1D::X::riemannFlux(u, uw);
		}

		U[it+1].iterateWhole(GRID1D_ITERATOR {
			U(it+1,ix)=U(it,ix) - dt/dx*(F(it,ix+0.5)-F(it,ix-0.5));
		});

		double mass=0;
		U[it+1].iterateWhole(GRID1D_ITERATOR {
            #pragma omp critical
            {
                mass += (U(it,ix)[0]*dx);
            }
		});

		DBGVAL(mass);






		if(isEveryNth(it,50)){

			gnuPlotSaver.save(U[it],frame("density_",it,"plot"),GRID1D_CALCULATOR{
				return HdVec1D::density(U(it,ix));
			});

			/*
			gnuPlotSaver.save(U[it],frame("pressure_",it,"plot"),GRID1D_CALCULATOR{
				return HdVec1D::pressure(U(it,ix));
			});

			gnuPlotSaver.save(U[it],frame("velocityX_",it,"plot"),GRID1D_CALCULATOR{
				return HdVec1D::velocityX(U(it,ix));
			});

			gnuPlotSaver.save(U[it],frame("Energy_",it,"plot"),GRID1D_CALCULATOR{
				return HdVec1D::fullEnergyPerVolumeUnit(U(it,ix));
			});

			gnuPlotSaver.save(F[it],frame("Flow_",it,"plot"),GRID1D_CALCULATOR{
				return F(it,ix)[0];
			})*/;

		}

	}
	return 0;
}

