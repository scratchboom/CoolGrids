#pragma once

template <typename T>
class TimedGrid1D{

private:
	Grid1D<T>* layers;
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
	double minIndexX,maxIndexX;

	std::string description;

	TimedGrid1D(const std::string& description){
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
		layers = new Grid1D<T>[layersCountToMaintain];
		for(int i=0;i<layersCountToMaintain;i++)layers[i].setDescription(description);
		for(int i=0;i<layersCountToMaintain;i++){
			layers[i].setRangeX(minX,maxX);
			layers[i].setIndexRangeX(minIndexX,maxIndexX);
		}
		for(int i=0;i<layersCountToMaintain;i++)layers[i].build();
	}

	~TimedGrid1D(){
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

	void setIndexRangeX(double minIndexX,double maxIndexX){
		this->minIndexX=minIndexX;
		this->maxIndexX=maxIndexX;
	}

	void fill(const T& value) {
		for (int i = 0; i < layersCountToMaintain; i++)	layers[i].fill(value);
	}

	Grid1D<T>& getAt(double it){
		if(it>largestAccessableLayer)switchToNextLayer();

        #ifdef GRIDS_CHECK_INDICES
		if(!inRangeInclusive(it,largestAccessableLayer-layersCountToMaintain+1,largestAccessableLayer)){
			std::cout << "ERROR: cannot access this layer " << it << std::endl;
			exit(1);
		}
        #endif

		int index=((int)(it-minIndexT)%layersCountToMaintain);
		if(index<0)index=layersCountToMaintain+index;
		return layers[index];

	}

	Grid1D<T>& operator[](double it){
		return getAt(it);
	}

	T& operator()(double it,double ix){
//		DBGLNVAL(it);
//		DBGLNVAL(ix);
		return getAt(it).getAt(ix);
	}

};



