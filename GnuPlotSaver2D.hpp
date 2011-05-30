#pragma once

class GnuPlotSaver2D{

private:
	std::string lineColor;
	double minValue,maxValue;
	bool autoScale;
	bool saveToPNG;
public:

	GnuPlotSaver2D(){
	}

	GnuPlotSaver2D& setLineColor(const std::string& lineColor){
		this->lineColor=lineColor;
		return *this;
	}

	GnuPlotSaver2D& setValueRange(double minValue,double maxValue){
		autoScale=false;
		this->minValue=minValue;
		this->maxValue=maxValue;
		return *this;
    }

	GnuPlotSaver2D& setAutoScale(){
		autoScale=true;
		return *this;
	}


	GnuPlotSaver2D& setSaveToPNG(bool saveToPNG){
		this->saveToPNG=saveToPNG;
		return *this;
	}

	template <typename T>
	void save(Grid2D<T>& grid,const std::string& filename){

		std::ofstream out(filename.c_str());

		if (saveToPNG) {
			out << "set term png" << "\n";
			out << "set output \"" << filename << ".png\"\n";
		}

		out << "set multiplot" << std::endl;
		out << "set title '" << grid.description << "'" << std::endl;

		out << "plot '-' matrix with image" << std::endl;

		for(double ix=grid.minIndexX;ix<grid.maxIndexX;ix++){
			for(double iy=grid.minIndexY;iy<grid.maxIndexY;iy++){
		    	out  << grid.getAt(ix,iy) << " ";
		    }
		    out << std::endl;
		}


	    out << "e" << std::endl;
	    out << "e" << std::endl;

		out.close();

		if (saveToPNG) {
			std::string gnuPlotCall("gnuplot ");
			gnuPlotCall.append(filename);
			system(gnuPlotCall.c_str());
		}
	}

	template<typename T,typename CL>
	void save(Grid2D<T>& grid, const std::string& filename,const CL& closure) {

		std::ofstream out(filename.c_str());

		if (saveToPNG) {
			out << "set term png" << "\n";
			out << "set output \"" << filename << ".png\"\n";
		}

		out << "set multiplot" << std::endl;
		out << "set title '" << grid.description << "'" << std::endl;

		out << "plot '-' matrix with image" << std::endl;

		for (double ix = grid.minIndexX; ix < grid.maxIndexX; ix++) {
			for (double iy = grid.minIndexY; iy < grid.maxIndexY; iy++) {
				out << closure(ix, iy) << " ";
			}
			out << std::endl;
		}

		out << "e" << std::endl;
		out << "e" << std::endl;

		out.close();

		if (saveToPNG) {
			std::string gnuPlotCall("gnuplot ");
			gnuPlotCall.append(filename);
			system(gnuPlotCall.c_str());
		}
	}

	template <typename T>
	void saveSliceZ(Grid3D<T>& grid,double iz,const std::string& filename){

		std::ofstream out(filename.c_str());

		if(saveToPNG){
			out << "set term png" << "\n";
			out << "set output \"" << filename << ".png\"\n";
		}

		out << "set multiplot" << std::endl;
		out << "set title '" << grid.description << "'" << std::endl;

		if(!autoScale) out << "set cbrange [" << minValue << ":" << maxValue << "]" << std::endl;

		out << "plot '-' matrix with image" << std::endl;

		for (double iy = grid.minIndexY; iy <= grid.maxIndexY; iy++){
			for (double ix = grid.minIndexX; ix <= grid.maxIndexX; ix++) {
				double value = grid.getAt(ix, iy, iz);
				out  << value << " ";
			}
		    out << std::endl;
		}

	    out << "e" << std::endl;
	    out << "e" << std::endl;

		out.close();

		if (saveToPNG) {
			std::string gnuPlotCall("gnuplot ");
			gnuPlotCall.append(filename);
			system(gnuPlotCall.c_str());
		}
	}


};
