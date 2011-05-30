#pragma once

#include "GridsCommon.hpp"

#define GRID2D_ITERATOR [&](double ix,double iy)
#define GRID2D_CALCULATOR [&](double ix,double iy)->double

template <typename T>
class Grid2D{

private:
	T* data;

	int getIndex(double ix,double iy){
		int ixIndex = (int)(ix - minIndexX);
		int iyIndex = (int)(iy - minIndexY);

		return ixIndex * nodesCountY + iyIndex;
    }

    void checkIndices(double ix,double iy){

    	std::string errorMessage;
    	if(ix < minIndexX || ix > maxIndexX || !isInteger(ix - minIndexX))
    	        errorMessage += "wrong index ix:" + toString(ix) + "\n";

    	if(iy < minIndexY || iy > maxIndexY || !isInteger(iy - minIndexY))
    	    	        errorMessage += "wrong index iy:" + toString(iy) + "\n";

    	if(!errorMessage.empty()){
    	    std::cout << errorMessage ;
    	    std::cout << "Grid info: " << getInfo();
    	    exit(1);
    	}
    }

public:

    int nodesCountX,nodesCountY;

    double minX,maxX;
    double minY,maxY;

    double minIndexX,maxIndexX;
    double minIndexY,maxIndexY;

    double dx,dy;

    bool periodicalBoundaryX;

    std::string description;

	Grid2D(const std::string& description){
		this->description=description;
		data=NULL;
	}

	Grid2D(){
		Grid2D("unnamed Grid2D");
	}

	~Grid2D(){
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

	void setIndexRangeX(double minIndexX,double maxIndexX){
		this->minIndexX=minIndexX;
		this->maxIndexX=maxIndexX;
	}

	void setIndexRangeY(double minIndexY,double maxIndexY){
		this->minIndexY=minIndexY;
		this->maxIndexY=maxIndexY;
	}


	void setPeriodicalBoundaryX(bool isPeriodical){
		periodicalBoundaryX=true;
	}

	void build(){
		std::cout << "building "+getInfo();
		nodesCountX=maxIndexX-minIndexX+1;
		nodesCountY=maxIndexY-minIndexY+1;

		dx=(maxX-minX)/(nodesCountX-1);
		dy=(maxY-minY)/(nodesCountY-1);

		data = new T[nodesCountX*nodesCountY];
		size_t dataSize=sizeof(T)*nodesCountX*nodesCountY;
		trackMemory(dataSize);
	}


	void clear(){
		for(int ix=0;ix<nodesCountX*nodesCountY;ix++)data[ix]=T();//TODO: rename index
	}

	void fill(const T& value){
		for(int ix=0;ix<nodesCountX*nodesCountY;ix++)data[ix]=value;//TODO: rename index
	}

	void setDescription(const std::string& description){
		this->description=description;
	}


	std::string getInfo(){
		return "Grid2D  " + description + "\n    (" + toString(minX) + ")..(" + toString(maxX) + ")\n    <"+ toString(minIndexX) + ">..<" + toString(maxIndexX)+">\n";
	}

	T& getAt(double ix,double iy){
        #ifdef GRIDS_CHECK_INDICES
	        checkIndices(ix,iy);
        #endif

        return data[getIndex(ix,iy)];
	}

	T& operator()(double ix,double iy){
		return getAt(ix,iy);
	}


	template<typename CL>
    void iterateWhole(const CL& closure){

        #pragma omp parallel for schedule(dynamic)
        for(int nix=0;nix<nodesCountX;nix++){
        	for(int niy=0;niy<nodesCountY;niy++){
        		double ix=minIndexX+nix;
        		double iy=minIndexY+niy;
        		closure(ix,iy);
        	}
        }
	}

	template<typename CL>
	void iterateInternal(const CL& closure){

        #pragma omp parallel for schedule(dynamic)
		for(int nix=0+1;nix<nodesCountX-1;nix++){
			for(int niy=0+1;niy<nodesCountY-1;niy++){
				double ix=minIndexX+nix;
				double iy=minIndexY+niy;
			    closure(ix,iy);
			}
		}
	}

	template<typename CL>
	void iterateInternal(int BORDER_NX, int BORDER_NY, const CL& closure) {

        #pragma omp parallel for schedule(dynamic)
		for (int nix = BORDER_NX; nix < nodesCountX - BORDER_NX; nix++)
			for (int niy = BORDER_NY; niy < nodesCountY - BORDER_NY; niy++) {

				double ix = minIndexX + nix;
				double iy = minIndexY + niy;

				closure(ix, iy);
			}
	}

	template<typename CL>
	void iterateBorder(const CL& closure){

        #pragma omp parallel for schedule(dynamic)
		for(int nix=0;nix<nodesCountX;nix++){
			double ix=minIndexX+nix;
		    closure(ix,minIndexY);
			closure(ix,maxIndexY);
		}

        #pragma omp parallel for schedule(dynamic)
		for(int niy=0;niy<nodesCountY;niy++){
			double iy=minIndexY+niy;
			closure(minIndexX,iy);
			closure(maxIndexX,iy);
		}
	}


	template<typename CL>
	void iterateBorderMinX(int BORDER_WIDTH,const CL& closure){

        #pragma omp parallel for schedule(dynamic)
	    for(int nix=0;nix<BORDER_WIDTH;nix++)
	        for(int niy=0;niy<nodesCountY;niy++){

				double ix=minIndexX+nix;
				double iy=minIndexY+niy;

				closure(ix,iy);
			}
	}



	template<typename CL>
	void iterateBorderMaxX(int BORDER_WIDTH,const CL& closure){

        #pragma omp parallel for schedule(dynamic)
	    for(int nix=nodesCountX-BORDER_WIDTH;nix<nodesCountX;nix++)
	        for(int niy=0;niy<nodesCountY;niy++){

				double ix=minIndexX+nix;
				double iy=minIndexY+niy;

				closure(ix,iy);
			}
	}


	template<typename CL>
	void iterateBorderMinY(int BORDER_WIDTH,const CL& closure){

        #pragma omp parallel for schedule(dynamic)
	    for(int nix=0;nix<nodesCountX;nix++)
	        for(int niy=0;niy<BORDER_WIDTH;niy++){

				double ix=minIndexX+nix;
				double iy=minIndexY+niy;

				closure(ix,iy);
			}
	}



	template<typename CL>
	void iterateBorderMaxY(int BORDER_WIDTH,const CL& closure){

        #pragma omp parallel for schedule(dynamic)
	    for(int nix=0;nix<nodesCountX;nix++)
	        for(int niy=nodesCountY-BORDER_WIDTH;niy<nodesCountY;niy++){

				double ix=minIndexX+nix;
				double iy=minIndexY+niy;

				closure(ix,iy);
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


};








