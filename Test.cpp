#include "GridsCommon.hpp"

int main() {

	/*
	 GnuPlotSaver2D gnuPlotSaver;
	 gnuPlotSaver.setLineColor("#FF0000")
	 .setValueRange(0.0,12.0)
	 .setSaveToPNG(true);

	 TimedGrid3D<double> U("Values vector");

	 const int N=32;

	 omp_set_num_threads(2);


	 U.setRangeX(0,N);
	 U.setRangeY(0,N);
	 U.setRangeZ(0,N);
	 U.setRangeT(0,1);

	 U.setIndexRangeX(0,N);
	 U.setIndexRangeY(0,N);
	 U.setIndexRangeZ(0,N);
	 U.setIndexRangeT(0,10);

	 U.build();


	 U[0].iterateWhole(GRID3D_ITERATOR{
	 DBGVAL(ix);
	 DBGVAL(iy);
	 DBGVAL(ix);
	 U(0,ix,iy,iz)=ix+2*iy+3*iz;
	 });

	 for(double it=0;it<U.maxIndexT;it++){

	 U[it].iterateWhole(GRID3D_ITERATOR{
	 double IX=((int)ix+1)%N;
	 if(ix<0)ix+=N;
	 //DBGVAL(IX);
	 U(it+1,ix,iy,iz)=U(it,IX,iy,iz);
	 });

	 gnuPlotSaver.saveSliceZ(U[it],N/2,frame("Hz_",it,"plot"));
	 }
	 */

	VtiSaver3D vtiSaver;
	CImgSaver2D cimgSaver;

	cimgSaver.setValueRange(4.5,5);


	const int N=32;

	Grid3D<double> grid;

	grid.setRangeX(0, N);
	grid.setRangeY(0, N);
	grid.setRangeZ(0, N);

	grid.setIndexRangeX(0, N);
	grid.setIndexRangeY(0, N);
	grid.setIndexRangeZ(0, N);

	grid.build();

	grid.fillAntialiased(5,GRID3D_CONDITION{
		return length(ix,iy,iz) < N;
	});

	vtiSaver.save(grid,"gridAntialiased.vti");
	cimgSaver.saveSliceZ(grid, N / 2, "gridAntialiased.png",GRID3D_CALCULATOR{
					return grid(ix,iy,iz);
				});

	return 0;
}
