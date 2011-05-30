
#include "GridsCommon.hpp"

int main(){

	GnuPlotSaver2D gnuPlotSaver;
	gnuPlotSaver.setLineColor("#FF0000")
			    .setValueRange(0.0,12.0)
			    .setSaveToPNG(true);

	CImgSaver2D cimgSaver;
	cimgSaver.setValueRange(0.8,2.4);

	double Nx=512;
	double Ny=512;

    double NT=20000;

	double Sx=1.0;

	double dx,dy;
	dx=dy=Sx/Nx;

	double SGM=0.5;
	double MAX_LAMBDA=5.0;
	double dt=dx/MAX_LAMBDA * SGM;

	TimedGrid2D<Vector4D> U("Hydrodynamics vector");
	TimedGrid2D<Vector4D> Fx("FlowX vector");
	TimedGrid2D<Vector4D> Fy("FlowY vector");

	U.setRangeX(0.5*dx,Nx*dx-0.5*dx);
	U.setIndexRangeX(0.5,Nx-0.5);

	U.setRangeY(0.5*dy,Ny*dy-0.5*dy);
	U.setIndexRangeY(0.5,Ny-0.5);
	U.build();

	Fx.setRangeX(0,dx*Nx);
	Fx.setIndexRangeX(0,Nx);
	Fx.setRangeY(0.5*dy,Ny*dy-0.5*dy);
	Fx.setIndexRangeY(0.5,Ny-0.5);
	Fx.build();

	Fy.setRangeX(0.5*dx,Nx*dx-0.5*dx);
    Fy.setIndexRangeX(0.5,Nx-0.5);
    Fy.setRangeY(0,dy*Ny);
    Fy.setIndexRangeY(0,Ny);
	Fy.build();


	U.fill(Vector4D(4));
	Fx.fill(Vector4D(4));
	Fy.fill(Vector4D(4));

	U[0].iterateWhole(GRID2D_ITERATOR{

		/*
		if(ix>Nx/2) U(0,ix,iy) = HdVec2D::fromDensityPressureVelocity(2,2,0,0);
		else U(0,ix,iy) = HdVec2D::fromDensityPressureVelocity(1,1,0,0);
		*/

		/*
		if(iy>Ny/2) U(0,ix,iy) = HdVec2D::fromDensityPressureVelocity(2,2,0,0);
	    else U(0,ix,iy) = HdVec2D::fromDensityPressureVelocity(1,1,0,0);
        */

		//corner
		if(ix+iy<Nx/2) U(0,ix,iy) = HdVec2D::fromDensityPressureVelocity(2,2,0,0);
		else U(0,ix,iy) = HdVec2D::fromDensityPressureVelocity(1,1,0,0);

	});


	for(double it=0;it<NT;it++){

		DBGVAL(it);

		Fx[it].iterateInternal(1,0,GRID2D_ITERATOR{
			Fx(it,ix,iy)=HdFlowVec2D::X::riemannFlux( U(it,ix-0.5,iy) , U(it,ix+0.5,iy) );
		});

		Fy[it].iterateInternal(0,1,GRID2D_ITERATOR{
			Fy(it,ix,iy)=HdFlowVec2D::Y::riemannFlux( U(it,ix,iy-0.5) , U(it,ix,iy+0.5) );
		});



		//transparent BC
		/*
		Fx[it].iterateBorderMinX(GRID2D_ITERATOR{
			Fx(it,ix,iy)=Fx(it,ix+1,iy);
		});


		Fx[it].iterateBorderMaxX(GRID2D_ITERATOR {
			Fx(it,ix,iy)=Fx(it,ix-1,iy);
		});


		Fy[it].iterateBorderMinY(GRID2D_ITERATOR{
			Fy(it,ix,iy)=Fy(it,ix,iy+1);
		});

		Fy[it].iterateBorderMaxY(GRID2D_ITERATOR {
			Fy(it,ix,iy)=Fy(it,ix,iy-1);
		});
		*/



		//reflexive BC

		Fx[it].iterateBorderMinX(GRID2D_ITERATOR{
			Vector4D u=U(it,ix+0.5,iy);

			double rho=HdVec2D::density(u);
			double p=HdVec2D::pressure(u);
			double v=HdVec2D::velocityX(u);

			Vector4D uw=HdVec2D::fromDensityPressureVelocity(rho,p,-v,0);

			Fx(it,ix,iy)=HdFlowVec2D::X::riemannFlux( uw , u );
		});


		Fx[it].iterateBorderMaxX(GRID2D_ITERATOR {
			Vector4D u=U(it,ix-0.5,iy);

			double rho=HdVec2D::density(u);
			double p=HdVec2D::pressure(u);
			double v=HdVec2D::velocityX(u);

			Vector4D uw=HdVec2D::fromDensityPressureVelocity(rho,p,-v,0);

			Fx(it,ix,iy)=HdFlowVec2D::X::riemannFlux( u , uw );
		});


		Fy[it].iterateBorderMinY(GRID2D_ITERATOR{
			Vector4D u=U(it,ix,iy+0.5);

			double rho=HdVec2D::density(u);
			double p=HdVec2D::pressure(u);
			double v=HdVec2D::velocityY(u);

			Vector4D uw=HdVec2D::fromDensityPressureVelocity(rho,p,0,-v);

			Fy(it,ix,iy)=HdFlowVec2D::Y::riemannFlux( uw , u );
		});

		Fy[it].iterateBorderMaxY(GRID2D_ITERATOR {
			Vector4D u=U(it,ix,iy-0.5);

			double rho=HdVec2D::density(u);
			double p=HdVec2D::pressure(u);
			double v=HdVec2D::velocityY(u);

			Vector4D uw=HdVec2D::fromDensityPressureVelocity(rho,p,0,-v);

			Fy(it,ix,iy)=HdFlowVec2D::Y::riemannFlux( u , uw );
		});





		U[it+1].iterateWhole(GRID2D_ITERATOR{
			U(it+1,ix,iy)=U(it,ix,iy) - dt/dx*(Fx(it,ix+0.5,iy)-Fx(it,ix-0.5,iy))
		                              - dt/dy*(Fy(it,ix,iy+0.5)-Fy(it,ix,iy-0.5));
		});

		DBGVAL(U(it,100.5,0.5));
		DBGVAL(Fy(it,100.5,0.0));
		DBGVAL(Fy(it,100.5,1.0));


		if(isEveryNth(it,100)){

			/*
			gnuPlotSaver.save(U[it],frame("density_",it,"plot"),GRID2D_CALCULATOR{
				return HdVec2D::density(U(it,ix,iy));
			});

			gnuPlotSaver.save(U[it],frame("pressure_",it,"plot"),GRID2D_CALCULATOR{
				return HdVec2D::pressure(U(it,ix,iy));
			});

			gnuPlotSaver.save(U[it],frame("velocityX_",it,"plot"),GRID2D_CALCULATOR{
				return HdVec2D::velocityX(U(it,ix,iy));
			});

			gnuPlotSaver.save(U[it],frame("Energy_",it,"plot"),GRID2D_CALCULATOR{
				return HdVec2D::fullEnergyPerVolumeUnit(U(it,ix,iy));
			});

			gnuPlotSaver.save(Fx[it],frame("FlowX_",it,"plot"),GRID2D_CALCULATOR{
				return Fx(it,ix,iy)[0];
			});*/


			cimgSaver.save(U[it],frame("density_",it,"png"),GRID2D_CALCULATOR{
				return HdVec2D::density(U(it,ix,iy));
			});


/*
			cimgSaver.save(U[it],frame("pressure_",it,"png"),GRID2D_CALCULATOR{
							return HdVec2D::pressure(U(it,ix,iy));
						});

			cimgSaver.save(U[it],frame("velocityX_",it,"png"),GRID2D_CALCULATOR{
							return HdVec2D::velocityX(U(it,ix,iy));
						});

			cimgSaver.save(U[it],frame("velocityY_",it,"png"),GRID2D_CALCULATOR{
							return HdVec2D::velocityY(U(it,ix,iy));
						});

			cimgSaver.save(U[it],frame("Energy_",it,"png"),GRID2D_CALCULATOR{
							return HdVec2D::fullEnergyPerVolumeUnit(U(it,ix,iy));
						});
*/

		}

	}
	return 0;
}

