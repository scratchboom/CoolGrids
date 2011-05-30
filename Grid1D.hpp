#pragma once

#include "GridsCommon.hpp"

#define GRID1D_ITERATOR [&](double ix)
#define GRID1D_CALCULATOR [&](double ix)->double

template <typename T>
class Grid1D{

private:
	T* data;

	int getIndex(double ix){
		return (int)(ix - minIndexX);
	}

    void checkIndices(double ix){

    	std::string errorMessage;
    	if(ix < minIndexX || ix > maxIndexX || !isInteger(ix - minIndexX)){
    	    errorMessage += "wrong index ix:" + toString(ix) + "\n";
    	    errorMessage += "sould be [" + toString(minIndexX) + ".." + toString(maxIndexX) + "]\n";
    	}


    	if(!errorMessage.empty()){
    	    std::cout << errorMessage ;
    	    std::cout << "Grid info: " << getInfo();
    	    exit(1);
    	}
    }

public:

    int nodesCountX;

    double minX,maxX;
    double minIndexX,maxIndexX;
    double dx;

    bool periodicalBoundaryX;

    std::string description;

	Grid1D(const std::string& description){
		this->description=description;
		data=NULL;
	}

	Grid1D(){
		Grid1D("unnamed Grid1D");
	}

	~Grid1D(){
		//std::cout << "destroying "+getInfo();
		delete[] data;
		data=NULL;
	}

	void setRangeX(double minX,double maxX){
		this->minX=minX;
		this->maxX=maxX;
	}

	void setIndexRangeX(double minIndexX,double maxIndexX){
		this->minIndexX=minIndexX;
		this->maxIndexX=maxIndexX;
	}

	void setPeriodicalBoundaryX(bool isPeriodical){
		periodicalBoundaryX=true;
	}

	void build(){
		std::cout << "building "+getInfo();
		nodesCountX=maxIndexX-minIndexX+1;
		dx=(maxX-minX)/(nodesCountX-1);

		data = new T[nodesCountX];
		size_t dataSize=sizeof(T)*nodesCountX;
		trackMemory(dataSize);
	}


	void clear(){
		for(int ix=0;ix<nodesCountX;ix++)data[ix]=T();
	}

	void fill(const T& value) {
		for (int ix = 0; ix < nodesCountX ; ix++)	data[ix] = value;//TODO: rename index
	}

	void setDescription(const std::string& description){
		this->description=description;
	}


	std::string getInfo(){
		return "Grid1D  " + description + "\n    (" + toString(minX) + ")..(" + toString(maxX) + ")\n    <"+ toString(minIndexX) + ">..<" + toString(maxIndexX)+">\n";
	}

	T& getAt(double ix){
        #ifdef GRIDS_CHECK_INDICES
	        checkIndices(ix);
        #endif

        return data[getIndex(ix)];
	}

	T& operator()(double ix){
		return getAt(ix);
	}


	template<typename CL>
    void iterateWhole(const CL& closure){

        for(int nix=0;nix<nodesCountX;nix++){
	        double ix=minIndexX+nix;
	        closure(ix);
        }
	}

	template<typename CL>
	void iterateInternal(const CL& closure){

		for(int nix=0+1;nix<nodesCountX-1;nix++){
			double ix=minIndexX+nix;
			closure(ix);
		}

	}

	template<typename CL>
	void iterateInternal(int BORDER_NX,const CL& closure){

		for(int nix=0+BORDER_NX;nix<nodesCountX-BORDER_NX;nix++){
			double ix=minIndexX+nix;
			closure(ix);
		}
	}

	template<typename CL>
	void iterateBorderMinX(const CL& closure){

		for(int nix=0;nix<1;nix++){
			double ix=minIndexX+nix;
			closure(ix);
		}
	}

	template<typename CL>
	void iterateBorderMaxX(const CL& closure) {

		for (int nix = nodesCountX-1; nix < nodesCountX; nix++) {
			double ix = minIndexX + nix;
			closure(ix);
		}
	}

};







