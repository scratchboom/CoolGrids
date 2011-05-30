
#pragma once

class CImgSaver2D{

private:
	double minValue,maxValue;
	//bool autoScale;
	//bool saveToPNG;
public:

	CImgSaver2D(){
	}


	CImgSaver2D& setValueRange(double minValue,double maxValue){
		//autoScale=false;
		this->minValue=minValue;
		this->maxValue=maxValue;
		return *this;
    }

	template <typename T,typename CL>
	void save(Grid2D<T>& grid,const std::string& filename,const CL& closure){

		Timer saveTimer;
		cimg_library::CImg<double> img(grid.nodesCountX, grid.nodesCountY, 1, 3);

		double amplitude=maxValue-minValue;

		for (double iy = grid.minIndexY; iy <= grid.maxIndexY; iy++)
			for (double ix = grid.minIndexX; ix <= grid.maxIndexX; ix++) {

				double value=closure(ix,iy);

				if(value>maxValue || value<minValue)value=rnd();
                else value=(value-minValue)/amplitude;

				img((int) (ix - grid.minIndexX), (int) (iy - grid.minIndexY), 0) = value * 360.0;
				img((int) (ix - grid.minIndexX), (int) (iy - grid.minIndexY), 1) = 1.0;
				img((int) (ix - grid.minIndexX), (int) (iy - grid.minIndexY), 2) = 1.0;
			}
		img.HSVtoRGB().save_png(filename.c_str());

		saveTimer.logTime(filename + " saved");
	}

	template <typename T>
	void save(Grid2D<T>& grid,const std::string& filename){
		save(grid,filename,GRID2D_CALCULATOR{
			return grid(ix,iy);
		});
	}


	template<typename T,typename CL>
	void saveSliceZ(Grid3D<T>& grid, double iz, const std::string& filename,const CL& closure) {
		Timer saveTimer;
		cimg_library::CImg<double> img(grid.nodesCountX, grid.nodesCountY, 1, 3);

		double amplitude = maxValue - minValue;

		for (double iy = grid.minIndexY; iy <= grid.maxIndexY; iy++)
			for (double ix = grid.minIndexX; ix <= grid.maxIndexX; ix++) {

				double value = closure(ix, iy, iz);

				if (value > maxValue || value < minValue) value = rnd();
				else value = (value - minValue) / amplitude;

				img((int) (ix - grid.minIndexX), (int) (iy - grid.minIndexY), 0) = value * 360.0;
				img((int) (ix - grid.minIndexX), (int) (iy - grid.minIndexY), 1) = 1.0;
				img((int) (ix - grid.minIndexX), (int) (iy - grid.minIndexY), 2) = 1.0;
			}

		img.HSVtoRGB().save_png(filename.c_str());
		saveTimer.logTime(filename + " saved");
	}

	template<typename T>
	void saveSliceZ(Grid3D<T>& grid, double iz, const std::string& filename) {
		saveSliceZ(grid,double iz,const std::string& filename,GRID3D_CALCULATOR{
			return grid(ix,iy,iz);
		});
	}


};

