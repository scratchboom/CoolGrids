#include <boost/numeric/ublas/Matrix.hpp>
#include "expokit.h"

template<typename T>
boost::numeric::ublas::matrix<T, column_major> exp(const boost::numeric::ublas::matrix<T, column_major> &A) {
	int ideg = 6;
	int m = 2;
	double H[] = { 1.2, 3.0, 5.6, 4 };

	int ldh = m;

	int lwsp = 4 * m * m + ideg + 1;
	double* wsp = new double[lwsp];

	int* ipiv = new int[m];

	double t = 1;
	int iexph;
	int ns;
	int iflag;

	std::cout << "dgexpv_" << std::endl;
	dgpadm_(&ideg, &m, &t, H, &ldh, wsp, &lwsp, ipiv, &iexph, &ns, &iflag);

	std::cout << wsp[iexph - 1] << std::endl;
	std::cout << wsp[iexph + 0] << std::endl;
	std::cout << wsp[iexph + 1] << std::endl;
	std::cout << wsp[iexph + 2] << std::endl;
	std::cout << "expokit library test" << std::endl;
	std::cout << "iflag: " << iflag << std::endl;

}
