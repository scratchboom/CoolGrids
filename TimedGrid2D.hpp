#pragma once

template <typename T>
class TimedGrid2D{

private:
	Grid2D<T>* layers;
	int layersCountToMaintain;
	double largestAccessableLayer;

	void switchToNextLayer(){
		largestAccessableLayer++;
	}

public:

	double minT,maxT;
	double minIndexT,maxIndexT;

	double minX,maxX;
	double minY,maxY;

	double minIndexX,maxIndexX;
	double minIndexY,maxIndexY;

	std::string description;

	TimedGrid2D(const std::string& description){
		this->description=description;

		largestAccessableLayer=0;
		layersCountToMaintain=2;
		layers=NULL;
	}

	void build(){
		largestAccessableLayer=minIndexT;
		layers = new Grid2D<T>[layersCountToMaintain];
		for(int i=0;i<layersCountToMaintain;i++)layers[i].setDescription(description);
		for(int i=0;i<layersCountToMaintain;i++){
			layers[i].setRangeX(minX,maxX);
			layers[i].setRangeY(minY,maxX);

			layers[i].setIndexRangeX(minIndexX,maxIndexX);
			layers[i].setIndexRangeY(minIndexY,maxIndexY);
		}
		for(int i=0;i<layersCountToMaintain;i++)layers[i].build();
	}

	~TimedGrid2D(){
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

	void setIndexRangeX(double minIndexX,double maxIndexX){
		this->minIndexX=minIndexX;
		this->maxIndexX=maxIndexX;
	}

	void setIndexRangeY(double minIndexY,double maxIndexY){
		this->minIndexY=minIndexY;
		this->maxIndexY=maxIndexY;
	}

	void fill(const T& value){
		for(int i=0;i<layersCountToMaintain;i++)layers[i].fill(value);
	}

	std::string getInfo(){
		return "TimedGrid3D  " + description + "\n    "+"    ["+toString(minIndexT)+".."+toString(maxIndexT)+"]\n    (" + toString(minX) + ")..(" + toString(maxX) + ")\n    <"+ toString(minIndexX) + ">..<" + toString(maxIndexX)+">\n"
				                                                                                                   +"(" + toString(minY) + ")..(" + toString(maxY) + ")\n    <"+ toString(minIndexY) + ">..<" + toString(maxIndexY)+">\n";
	}

	Grid2D<T>& getAt(double it){
		if(it>largestAccessableLayer)switchToNextLayer();

		if((largestAccessableLayer-it+1)>layersCountToMaintain || (largestAccessableLayer-it<0)){
			std::cout << "ERROR: cannot access this layer " << it << std::endl;
			exit(1);
		}

		return layers[(int)(it-minIndexT)%layersCountToMaintain];

	}

	Grid2D<T>& operator[](double it){
		return getAt(it);
	}

	T& operator()(double it,double ix,double iy){
		return getAt(it).getAt(ix,iy);
	}

};



