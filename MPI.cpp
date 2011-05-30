#include <iostream>
#include "mpi.h"

using namespace std;

int main(int argc, char* argv[]) {

	for(int i=0;i<argc;i++)cout << argv[i] << endl;

	int	numtasks, taskid;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

	cout << "numtasks: " << numtasks << endl;
	cout << "task id: " << taskid << endl;

	MPI_Finalize();

	return 0;
}
