#pragma once

#include "GridsCommon.hpp"

#define GRID3D_ITERATOR [&](double ix,double iy,double iz)
#define GRID3D_CALCULATOR [&](double ix,double iy,double iz)->double
#define GRID3D_CONDITION [&](double ix,double iy,double iz)->bool

template <typename T>
class Grid3D{

private:
	T* data;

	int getIndex(double ix,double iy,double iz){
		int ixIndex = (int)(ix - minIndexX);
		int iyIndex = (int)(iy - minIndexY);
		int izIndex = (int)(iz - minIndexZ);

		return ixIndex * nodesCountY*nodesCountZ +
		       iyIndex * nodesCountZ +
			   izIndex;
    }

    void checkIndices(double ix,double iy,double iz){

    	std::string errorMessage;
    	if(ix < minIndexX || ix > maxIndexX || !isInteger(ix - minIndexX))
    		errorMessage += "wrong index ix:" + toString(ix) + "\n";

    	if(iy < minIndexY || iy > maxIndexY || !isInteger(iy - minIndexY))
    	    	        errorMessage += "wrong index iy:" + toString(iy) + "\n";

    	if(iz < minIndexZ || iz > maxIndexZ || !isInteger(iz - minIndexZ))
    	    	    	        errorMessage += "wrong index iz:" + toString(iz) + "\n";

    	if(!errorMessage.empty()){
    	    std::cout << errorMessage ;
    	    std::cout << getInfo();
    	    exit(1);
    	}
    }

public:

    int nodesCountX,nodesCountY,nodesCountZ;

    double minX,maxX;
    double minY,maxY;
    double minZ,maxZ;

    double minIndexX,maxIndexX;
    double minIndexY,maxIndexY;
    double minIndexZ,maxIndexZ;

    double dx,dy,dz;

    bool periodicalBoundaryX;

    std::string description;

	Grid3D(const std::string& description){
		this->description=description;
		data=NULL;
	}

	Grid3D(){
		Grid3D("unnamed Grid3D");
	}

	~Grid3D(){
		//std::cout << "destroying "+getInfo();
		delete[] data;
		data=NULL;
	}

	void setRangeX(double minX,double maxX){
		this->minX=minX;
		this->maxX=maxX;
	}

	void setRangeY(double minY,double maxY){
		this->minY=minY;
		this->maxY=maxY;
	}

	void setRangeZ(double minZ,double maxZ){
		this->minZ=minZ;
		this->maxZ=maxZ;
	}

	void setIndexRangeX(double minIndexX,double maxIndexX){
		this->minIndexX=minIndexX;
		this->maxIndexX=maxIndexX;
	}

	void setIndexRangeY(double minIndexY,double maxIndexY){
		this->minIndexY=minIndexY;
		this->maxIndexY=maxIndexY;
	}

	void setIndexRangeZ(double minIndexZ,double maxIndexZ){
		this->minIndexZ=minIndexZ;
		this->maxIndexZ=maxIndexZ;
	}

	void setPeriodicalBoundaryX(bool isPeriodical){
		periodicalBoundaryX=true;
	}

	void buildBounds(){
		std::cout << "building bounds"+getInfo();
		nodesCountX=maxIndexX-minIndexX+1;
		nodesCountY=maxIndexY-minIndexY+1;
		nodesCountZ=maxIndexZ-minIndexZ+1;

		dx=(maxX-minX)/(nodesCountX-1);
		dy=(maxY-minY)/(nodesCountY-1);
		dz=(maxZ-minZ)/(nodesCountZ-1);
	}

	void build(){
		//DBGLNVAL(this);
		buildBounds();

		data = new T[nodesCountX*nodesCountY*nodesCountZ];
		size_t dataSize=sizeof(T)*nodesCountX*nodesCountY*nodesCountZ;
		trackMemory(dataSize);
	}


	void clear(){
		for(int ix=0;ix<nodesCountX*nodesCountY*nodesCountZ;ix++)data[ix]=T();//TODO: rename index
	}

	void fill(const T& value) {
		for (int ix = 0; ix < nodesCountX*nodesCountY*nodesCountZ; ix++) data[ix] = value;//TODO: rename index
	}

	void setDescription(const std::string& description){
		this->description=description;
	}


	std::string getInfo(){
		return "Grid3D  " + description + "\n    (" + toString(minX) + ")..(" + toString(maxX) + ")\n    <"+ toString(minIndexX) + ">..<" + toString(maxIndexX)+">\n"
				                               +"(" + toString(minY) + ")..(" + toString(maxY) + ")\n    <"+ toString(minIndexY) + ">..<" + toString(maxIndexY)+">\n"
				                               +"(" + toString(minZ) + ")..(" + toString(maxZ) + ")\n    <"+ toString(minIndexZ) + ">..<" + toString(maxIndexZ)+">\n";
	}

	T& getAt(double ix,double iy,double iz){
        #ifdef GRIDS_CHECK_INDICES
	        checkIndices(ix,iy,iz);
        #endif

        return data[getIndex(ix,iy,iz)];
	}

	T& operator()(double ix,double iy,double iz){
		return getAt(ix,iy,iz);
	}


	template<typename CL>
    void iterateWhole(const CL& closure){

        #pragma omp parallel for
        for(int nix=0;nix<nodesCountX;nix++){
        	for(int niy=0;niy<nodesCountY;niy++){
        		for(int niz=0;niz<nodesCountZ;niz++){
        		    double ix=minIndexX+nix;
        		    double iy=minIndexY+niy;
        		    double iz=minIndexZ+niz;
        		    closure(ix,iy,iz);
        		}
        	}
        }
	}

	template<typename CL>
	void fill(const T& value,const CL& condition){
		iterateWhole(GRID3D_ITERATOR{
			if(condition(ix,iy,iz)) getAt(ix,iy,iz) = value;
		});
	}

	template<typename CL>
	void fillAntialiased(const T& value,const CL& condition){
		double ANTIALIASING_SAMPLES_COUNT = 20.0;

		iterateWhole(GRID3D_ITERATOR{
			T sum = T();

			for(int i=0;i<ANTIALIASING_SAMPLES_COUNT;i++){
				double jittX=rnd()-0.5;
				double jittY=rnd()-0.5;
				double jittZ=rnd()-0.5;

				if(condition(ix+jittX,iy+jittY,iz+jittZ)) sum += value;
				else sum += getAt(ix,iy,iz);
			}

			getAt(ix,iy,iz) = sum/ANTIALIASING_SAMPLES_COUNT;
		});
	}

	template<typename CL>
	void iterateInternal(const CL& closure){

        #pragma omp parallel for
		for(int nix=0+1;nix<nodesCountX-1;nix++){
		    for(int niy=0+1;niy<nodesCountY-1;niy++){
		      	for(int niz=0+1;niz<nodesCountZ-1;niz++){
		            double ix=minIndexX+nix;
		            double iy=minIndexY+niy;
		            double iz=minIndexZ+niz;
		            closure(ix,iy,iz);
		        }
		    }
		}
	}

	template<typename CL>
	void iterateInternal(int BORDER_NX,int BORDER_NY,int BORDER_NZ,const CL& closure){

        #pragma omp parallel for
		for(int nix=BORDER_NX;nix<nodesCountX-BORDER_NX;nix++)
			for(int niy=BORDER_NY;niy<nodesCountY-BORDER_NY;niy++)
				for(int niz=BORDER_NZ;niz<nodesCountZ-BORDER_NZ;niz++){
					double ix=minIndexX+nix;
					double iy=minIndexY+niy;
					double iz=minIndexZ+niz;

					closure(ix,iy,iz);
			}
	}

	template<typename CL>
	void iterateBorder(const CL& closure){

        #pragma omp parallel for
		for(int nix=0;nix<nodesCountX;nix++)
			for(int niy=0;niy<nodesCountY;niy++){
				double ix=minIndexX+nix;
				double iy=minIndexY+niy;

				closure(ix,iy,minIndexZ);
				closure(ix,iy,maxIndexZ);
		    }


        #pragma omp parallel for
		for(int nix=0;nix<nodesCountX;nix++)
			for(int niz=0;niz<nodesCountZ;niz++){
				double ix=minIndexX+nix;
				double iz=minIndexZ+niz;

				closure(ix,minIndexY,iz);
				closure(ix,maxIndexY,iz);
			}

        #pragma omp parallel for
        for(int niy=0;niy<nodesCountY;niy++)
	        for(int niz=0;niz<nodesCountZ;niz++){
	        	double iy=minIndexY+niy;
	        	double iz=minIndexZ+niz;

	        	closure(minIndexX,iy,iz);
	        	closure(maxIndexX,iy,iz);
	        }
	}


	template<typename CL>
	void iterateBorderMinX(int BORDER_WIDTH,const CL& closure){

        #pragma omp parallel for
	    for(int nix=0;nix<BORDER_WIDTH;nix++)
	        for(int niy=0;niy<nodesCountY;niy++)
			   for(int niz=0;niz<nodesCountZ;niz++){

				double ix=minIndexX+nix;
				double iy=minIndexY+niy;
				double iz=minIndexZ+niz;

				closure(ix,iy,iz);
			}
	}


	template<typename CL>
	void iterateBorderMaxX(int BORDER_WIDTH,const CL& closure){

        #pragma omp parallel for
	    for(int nix=nodesCountX-BORDER_WIDTH;nix<nodesCountX;nix++)
	        for(int niy=0;niy<nodesCountY;niy++)
			   for(int niz=0;niz<nodesCountZ;niz++){

				double ix=minIndexX+nix;
				double iy=minIndexY+niy;
				double iz=minIndexZ+niz;

				closure(ix,iy,iz);
			}
	}


	template<typename CL>
	void iterateBorderMinY(int BORDER_WIDTH,const CL& closure){

        #pragma omp parallel for
	    for(int nix=0;nix<nodesCountX;nix++)
	        for(int niy=0;niy<BORDER_WIDTH;niy++)
			   for(int niz=0;niz<nodesCountZ;niz++){

				double ix=minIndexX+nix;
				double iy=minIndexY+niy;
				double iz=minIndexZ+niz;

				closure(ix,iy,iz);
			}
	}

	template<typename CL>
	void iterateBorderMaxY(int BORDER_WIDTH,const CL& closure){

        #pragma omp parallel for
	    for(int nix=0;nix<nodesCountX;nix++)
	        for(int niy=nodesCountY-BORDER_WIDTH;niy<nodesCountY;niy++)
			   for(int niz=0;niz<nodesCountZ;niz++){

				double ix=minIndexX+nix;
				double iy=minIndexY+niy;
				double iz=minIndexZ+niz;

				closure(ix,iy,iz);
			}
	}


	template<typename CL>
	void iterateBorderMinZ(int BORDER_WIDTH,const CL& closure){

        #pragma omp parallel for
	    for(int nix=0;nix<nodesCountX;nix++)
	        for(int niy=0;niy<nodesCountY;niy++)
			   for(int niz=0;niz<BORDER_WIDTH;niz++){

				double ix=minIndexX+nix;
				double iy=minIndexY+niy;
				double iz=minIndexZ+niz;

				closure(ix,iy,iz);
			}
	}


	template<typename CL>
	void iterateBorderMaxZ(int BORDER_WIDTH,const CL& closure){

        #pragma omp parallel for
	    for(int nix=0;nix<nodesCountX;nix++)
	        for(int niy=0;niy<nodesCountY;niy++)
			   for(int niz=nodesCountZ-BORDER_WIDTH;niz<nodesCountZ;niz++){

				double ix=minIndexX+nix;
				double iy=minIndexY+niy;
				double iz=minIndexZ+niz;

				closure(ix,iy,iz);
			}
	}


	template<typename CL>
	void iterateBorderMinX(const CL& closure){
	    iterateBorderMinX(1,closure);
	}

	template<typename CL>
	void iterateBorderMaxX(const CL& closure){
	    iterateBorderMaxX(1,closure);
	}

	template<typename CL>
	void iterateBorderMinY(const CL& closure){
	    iterateBorderMinY(1,closure);
	}


	template<typename CL>
	void iterateBorderMaxY(const CL& closure){
	    iterateBorderMaxY(1,closure);
	}

	template<typename CL>
	void iterateBorderMinZ(const CL& closure){
	    iterateBorderMinZ(1,closure);
	}

	template<typename CL>
	void iterateBorderMaxZ(const CL& closure){
	    iterateBorderMaxZ(1,closure);
	}


};


class VtkSaver3D{

private:
	std::string lineColor;
	double minValue,maxValue;
	bool autoScale;
public:

	VtkSaver3D(){
	}

	VtkSaver3D& setLineColor(const std::string& lineColor){
		this->lineColor=lineColor;
		return *this;
	}

	VtkSaver3D& setValueRange(double minValue,double maxValue){
		autoScale=false;
		this->minValue=minValue;
		this->maxValue=maxValue;
		return *this;
    }

	VtkSaver3D& setAutoScale(double minValue,double maxValue){
		autoScale=true;
		return *this;
	}

	template <typename T>
	void save(Grid3D<T>& grid,const std::string& filename){
		std::ofstream out(filename.c_str());

		Timer saveTimer;

		out << "# vtk DataFile Version 3.0\n";
		out << "Insert Any Information Here.\n";
		out << "ASCII\n";
		out << "DATASET STRUCTURED_GRID\n";
		out << "DIMENSIONS "<< (int)grid.nodesCountX <<" "<< (int)grid.nodesCountY  <<" "<< (int)grid.nodesCountZ  <<"\n";
		out << "POINTS"<<" "<< grid.nodesCountX*grid.nodesCountY*grid.nodesCountZ <<" double\n";

		for (double iz = grid.minIndexZ; iz <= grid.maxIndexZ; iz++)
			for (double iy = grid.minIndexY; iy <= grid.maxIndexY; iy++)
				for (double ix = grid.minIndexX; ix <= grid.maxIndexX; ix++) {
					out << grid.minX+(ix-grid.minIndexX)*grid.dx << " " << grid.minY+(iy-grid.minIndexY)*grid.dy << " " << grid.minZ+(iz-grid.minIndexZ)*grid.dz <<"\n";
				}

		out << "POINT_DATA " << grid.nodesCountX*grid.nodesCountY*grid.nodesCountZ <<"\n";
		out << "SCALARS scalar_name float 1\n";
		out << "LOOKUP_TABLE default\n";

		for (double iz = grid.minIndexZ; iz <= grid.maxIndexZ; iz++)
			for (double iy = grid.minIndexY; iy <= grid.maxIndexY; iy++)
				for (double ix = grid.minIndexX; ix <= grid.maxIndexX; ix++) {
					out << grid.getAt(ix, iy, iz)<<"\n";
				}

		out.close();

		saveTimer.logTime(filename+" saved");
	}

};

