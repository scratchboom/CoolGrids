#include "GridsCommon.hpp"

using namespace std;

template<typename T>
class MpiGrid3D: Grid3D<T> {

private:
	Grid3D<T> fullGridBounds;
	Grid3D<T> subGridBounds;

	int partsX, partsY, partsZ;

public:

	MpiGrid3D() {
		Grid3D<T>();
		partsX = partsY = partsZ = 1;
	}

	void setParts(int partsX, int partsY, int partsZ) {
		this->partsX = partsX;
		this->partsY = partsY;
		this->partsZ = partsZ;
	}

	void build() {
		std::cout << "MpiGrid3D::build" << std::endl;
		buildBounds();

		if (nodesCountX % partsX == 0) {
			std::cout << nodesCountX << " = " << (nodesCountX / partsX) << "*"
					<< partsX << std::endl;
		} else {
			DBGVAL(nodesCountX % partsX);
			DBGVAL(nodesCountX / partsX);
			DBGVAL(nodesCountX / partsX +1);
		}

	}

};

double minFunc1(int Nx,int Ny,int Nz,int nx,int ny,int nz){
	return sqr((double) nx / ny - (double) Nx / Ny)
			+ sqr((double) ny / nx - (double) Ny / Nx)
			+ sqr((double) nx / nz - (double) Nx / Nz)
			+ sqr((double) nz / nx - (double) Nz / Nx)
			+ sqr((double) ny / nz - (double) Ny / Nz)
			+ sqr((double) nz / ny - (double) Nz / Ny);
}

double minFunc2(int Nx,int Ny,int Nz,int nx,int ny,int nz){
	return (nx-1)*Ny*Nz + (ny-1)*Nz*Nx +  (nz-1)*Nx*Ny;
}

void calcDims(int Nx, int Ny, int Nz, int n, int& nx, int& ny, int& nz) {
	double bestI = 1E50;

	int nxBest = 0;
	int nyBest = 0;
	int nzBest = 0;

	for (int nx = 1; nx <= n; nx++)
		for (int ny = 1; ny <= n; ny++) {

			int nz = (n / nx) / ny;

			if (nx * ny * nz != n)
				continue;


			double I = minFunc2(Nx,Ny,Nz,nx,ny,nz);

			if (I < bestI) {
				bestI = I;
				nxBest = nx;
				nyBest = ny;
				nzBest = nz;
			}
		}

	if (bestI < 1E50) {
		nx = nxBest;
		ny = nyBest;
		nz = nzBest;

		cout << "dimns calculated" << endl;
		cout << "nx: " << nxBest << endl;
		cout << "ny: " << nyBest << endl;
		cout << "nz: " << nzBest << endl;

		cout << "min proportion: " << minFunc1(Nx,Ny,Nz,nx,ny,nz) << endl;
		cout << "min area: " << minFunc2(Nx,Ny,Nz,nx,ny,nz) << endl;

	} else {
		cout << "dims calculation failed" << endl;
		throw "dims calculation failed";
	}
}

int main(int argc, char* argv[]) {

	cout << "hello there" << endl;

	int Nx = 100;
	int Ny = 100;
	int Nz = 200;

	int n = 1000;

	int nx,ny,nz;

	calcDims(Nx,Ny,Nz,n,nx,ny,nz);











	for(int i=0;i<argc;i++)cout << argv[i] << endl;

	 MpiEnvironment env(argc, argv);
	 MpiCommunicator world;

	 if(world.isMainProcess())cout << "numtasks: " << world.getSize() << endl;
	 cout << "task id: " << world.getRank() << endl;

	 if(world.isMainProcess()){
	 MpiGrid3D<double> g;
	 g.setIndexRangeX(0,2);
	 g.build();

	 }

	return 0;
}
















