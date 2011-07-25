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

	CsvTableSaver csvts;

	for(double x=-10.0;x<=10.0;x+=0.5){
		csvts.addRow("x",x);
		csvts.addRow("x^2",x*x);
		csvts.addRow("x^3",x*x*x);
	}

	csvts.save("csvts.csv");

	return 0;
}
