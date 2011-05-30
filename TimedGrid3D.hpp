#pragma once
#include "GridsCommon.hpp"

template <typename T>
class TimedGrid3D{

private:
	Grid3D<T>* layers;
	int layersCountToMaintain;
	double largestAccessableLayer;

	void switchToNextLayer(){
        #pragma omp atomic
		largestAccessableLayer++;
	}

public:

	double minT,maxT;
	double minIndexT,maxIndexT;

	double minX,maxX;
	double minY,maxY;
	double minZ,maxZ;

	double minIndexX,maxIndexX;
	double minIndexY,maxIndexY;
	double minIndexZ,maxIndexZ;

	std::string description;

	TimedGrid3D(const std::string& description){
		this->description=description;

		largestAccessableLayer=0;
		layersCountToMaintain=2;
		layers=NULL;
	}

	void setLayersCountToMaintain(int layersCountToMaintain){
		this->layersCountToMaintain=layersCountToMaintain;
	}

	void build(){
		largestAccessableLayer=minIndexT;
		layers = new Grid3D<T>[layersCountToMaintain];
		for(int i=0;i<layersCountToMaintain;i++)layers[i].setDescription(description);
		for(int i=0;i<layersCountToMaintain;i++){
			layers[i].setRangeX(minX,maxX);
			layers[i].setRangeY(minY,maxX);
			layers[i].setRangeZ(minZ,maxZ);

			layers[i].setIndexRangeX(minIndexX,maxIndexX);
			layers[i].setIndexRangeY(minIndexY,maxIndexY);
			layers[i].setIndexRangeZ(minIndexZ,maxIndexZ);
		}
		for(int i=0;i<layersCountToMaintain;i++)layers[i].build();
	}

	~TimedGrid3D(){
		delete[] layers;
	}

	void setRangeT(double minT,double maxT){
		this->minT = minT;
		this->maxT = maxT;
	}

	void setIndexRangeT(double minIndexT,double maxIndexT){
		this->minIndexT = minIndexT;
		this->maxIndexT = maxIndexT;
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

	void clear(){
		for(int i=0;i<layersCountToMaintain;i++)layers[i].clear();
	}

	void fill(const T& value) {
		for (int i = 0; i < layersCountToMaintain; i++)	layers[i].fill(value);
	}

	std::string getInfo(){
		return "TimedGrid3D  " + description + "\n    "+"    ["+toString(minIndexT)+".."+toString(maxIndexT)+"]\n    (" + toString(minX) + ")..(" + toString(maxX) + ")\n    <"+ toString(minIndexX) + ">..<" + toString(maxIndexX)+">\n"
				                                                                                                   +"(" + toString(minY) + ")..(" + toString(maxY) + ")\n    <"+ toString(minIndexY) + ">..<" + toString(maxIndexY)+">\n"
				                                                                                                   +"(" + toString(minZ) + ")..(" + toString(maxZ) + ")\n    <"+ toString(minIndexZ) + ">..<" + toString(maxIndexZ)+">\n";
	}

	Grid3D<T>& getAt(double it){
		if(it>largestAccessableLayer)switchToNextLayer();

        #ifdef GRIDS_CHECK_INDICES
		if(!inRangeInclusive(it,largestAccessableLayer-layersCountToMaintain+1,largestAccessableLayer)){
			DBGVAL(this);
			DBGVAL(largestAccessableLayer-layersCountToMaintain+1);
			DBGVAL(largestAccessableLayer);
			std::cout << "ERROR: cannot access this layer " << it << std::endl;
			std::cout << getInfo();
			exit(1);
		}
        #endif

		int index=((int)(it-minIndexT)%layersCountToMaintain);
		if(index<0)index=layersCountToMaintain+index;
		return layers[index];

	}

	Grid3D<T>& operator[](double it){
		return getAt(it);
	}

	T& operator()(double it,double ix,double iy,double iz){
		return getAt(it).getAt(ix,iy,iz);
	}

};



