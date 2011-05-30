#pragma once

class BmpSaver{
private:

	double minValue,maxValue;

public:

	BmpSaver(){
		minValue=0.0;
		maxValue=1.0;
	}

    void setValueRange(double minValue,double maxValue){
    	this->minValue=minValue;
    	this->maxValue=maxValue;
    }

	template<typename T>
	void saveSliceZ(Grid3D<T>& grid,double iz,const std::string& fileName){
		BmpImage img(grid.nodesCountX,grid.nodesCountY);

		double amplitude=maxValue-minValue;

		for (double iy = grid.minIndexY; iy <= grid.maxIndexY; iy++)
			for (double ix = grid.minIndexX; ix <= grid.maxIndexX; ix++) {

				double value=grid.getAt(ix,iy,iz);

				if(value>maxValue || value<minValue)value=rnd();
				else value=(value-minValue)/amplitude;

				img.setPixel((int)(ix-grid.minIndexX),(int)(iy-grid.minIndexY),value,value,value);
			}

		img.save(fileName);
	}


	template<typename T>
	void save(Grid2D<T>& grid,const std::string& fileName){
		BmpImage img(grid.nodesCountX,grid.nodesCountY);

		double amplitude=maxValue-minValue;

		for (double iy = grid.minIndexY; iy <= grid.maxIndexY; iy++)
			for (double ix = grid.minIndexX; ix <= grid.maxIndexX; ix++) {

				double value=grid.getAt(ix,iy);

				if(value>maxValue || value<minValue)value=rnd();
				else value=(value-minValue)/amplitude;

				RGB rgb=redBlueColorMap.getColor(value);

				img.setPixel((int)(ix-grid.minIndexX),(int)(iy-grid.minIndexY),rgb.R,rgb.G,rgb.B);
			}

		img.save(fileName);
	}
};
