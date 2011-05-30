#pragma once

class MpiEnvironment{
public:
	MpiEnvironment(int argc, char* argv[]){
		MPI_Init(&argc, &argv);
	}

	~MpiEnvironment(){
		MPI_Finalize();
	}
};

class MpiCommunicator{

private:
    MPI_Comm comm;

public:

    MpiCommunicator(){
    	comm = MPI_COMM_WORLD;
    }

    int isMainProcess(){
    	return getRank()==0;
    }

    int getRank(){
    	int rank;
    	MPI_Comm_rank(comm, &rank);
    	return rank;
    }

    int getSize(){
    	int size;
    	MPI_Comm_size(comm, &size);
    	return size;
    }

};
