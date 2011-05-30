#include "GridsCommon.hpp"


int main(){
	Grid1D<double> grid;

	grid.setRangeX(0,1);
	grid.setIndexRangeX(0,100);

	grid.build();

	grid.iterateWhole(GRID1D_ITERATOR{
		grid(ix)=sin(10.0*ix/grid.maxIndexX);
	});


	CsvSaver csv;
	csv.save(grid,"grid.csv");

	return 0;
}
