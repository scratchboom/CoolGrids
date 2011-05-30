#include "GridsCommon.hpp"

int main(){

	//auto pwcf=createPointwiseChargeField();

	double e=0.00001;
	double C=2.5;

	double t=0;
	double DT=0.1;

	omp_set_num_threads(2);

	CImgSaver2D cimgSaver;
	cimgSaver.setValueRange(-0.5,0.5);

	double Sx=40;
	double Sy=40;

	double Nx=512;
	double Ny=(int) Nx*Sy/Sx;

	Grid2D<double> Ex("Ex");

	Ex.setRangeX(-0.5*Sx,0.5*Sx);
	Ex.setRangeY(-0.5*Sy,0.5*Sy);

	Ex.setIndexRangeX(0,Nx);
	Ex.setIndexRangeY(0,Ny);

	Ex.build();


	Grid2D<double> Hz("Hz");

	Hz.setRangeX(-0.5*Sx,0.5*Sx);
	Hz.setRangeY(-0.5*Sy,0.5*Sy);

	Hz.setIndexRangeX(0,Nx);
	Hz.setIndexRangeY(0,Ny);

	Hz.build();

	double hx=Ex.dx;
	double hy=Ex.dy;

	auto field1=createPointwiseChargeField(1000,[&](double t){return 10.0*sin(1.0/20.0*t)*gaussStep(t,3.5,2.0);},
     			                             [&](double t){return 10.5*cos(1.5/20.0*t)*gaussStep(t,3.5,2.0);},
 			                                 [&](double t){return 0;});

	auto field=superposition(field1,createPointwiseChargeField(-1000,[&](double t){return -10.0*sin(1.0/20.0*t)*gaussStep(t,3.5,2.0);},
							                               [&](double t){return -10.0*cos(1.5/20.0*t)*gaussStep(t,3.5,2.0);},
						                                   [&](double t){return 0;})
						                                   );



	Timer mainTimer;
	for(double it=0;it<1000;it++){
		double t=it*DT;
		mainTimer.logTime(toString(it));

		Accumulator<double> acc;


		Ex.iterateWhole(GRID2D_ITERATOR{
			double x=Ex.minX+ix*hx;
			double y=Ex.minY+iy*hy;
			double z=0;

			EH eh=field->getField(t,DT,x,y,z);

			Ex(ix,iy) = eh.Ex;
			Hz(ix,iy) = eh.Hz;

			acc.add(eh.Hz);

		});

		DBGVAL(acc.getMaxValue());
		DBGVAL(acc.getMinValue());

		if(isEveryNth(it,1)){


		    cimgSaver.save(Ex,frame("./Ex_",it,"png"));
		    //cimgSaver.save(Hz,frame("./Hz_",it,"png"));
		}
	}

	return 0;
}

