#include <iostream>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;

int main() {

	matrix<double,column_major> colMajorMatrix (3, 3);
	matrix<double,row_major> rowMajorMatrix (3, 3);

	for (int i = 0; i < 3 ; i++)
	    for (int j = 0; j < 3 ; j++){
	    	colMajorMatrix (i, j) = 3 * i + j;
	    	rowMajorMatrix (i, j) = 3 * i + j;
	    }

	std::cout << colMajorMatrix << std::endl;
	std::cout << rowMajorMatrix << std::endl;

	std::cout << "col-major: ";
	for(int i=0;i<3*3;i++)std::cout << colMajorMatrix.data().begin()[i] <<"   ";
	std::cout << std::endl;
	std::cout << "row-major: ";
	for(int i=0;i<3*3;i++)std::cout << rowMajorMatrix.data().begin()[i] <<"   ";
	std::cout << std::endl;


	double *data=new double[22];
	matrix<double,row_major> wrapedMatrix (3, 3,data[2]);/*
	for (int i = 0; i < 3 ; i++)
		    for (int j = 0; j < 3 ; j++){
		    	wrapedMatrix (i, j) = 3 * i + j;
		    }*/

	for(int i=0;i<22;i++)std::cout << data[i] <<"   ";
	std::cout << std::endl;

	std::cout << wrapedMatrix << std::endl;



}

