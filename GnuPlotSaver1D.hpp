#pragma once

class GnuPlotSaver1D{

private:
	std::string lineColor;
	double minValue,maxValue;
	int xRes,yRes;
	bool autoScale;
	bool saveToPNG;
public:

	GnuPlotSaver1D(){
		this->lineColor="#FF0000";
		xRes = 640;
		yRes = 480;

		saveToPNG=true;
	}

    //sould be in format #FF0000
	GnuPlotSaver1D& setLineColor(const std::string& lineColor){
		this->lineColor=lineColor;
		return *this;
	}

	GnuPlotSaver1D& setValueRange(double minValue,double maxValue){
		autoScale=false;
		this->minValue=minValue;
		this->maxValue=maxValue;
		return *this;
    }

	GnuPlotSaver1D& setAutoScale(){
		autoScale=true;
		return *this;
	}

	GnuPlotSaver1D& setSaveToPNG(bool saveToPNG){
		this->saveToPNG=saveToPNG;
		return *this;
	}

	GnuPlotSaver1D& setSize(int xRes,int yRes){
		this->xRes=xRes;
		this->yRes=yRes;
		return *this;
	}

	template <typename T,typename CL>
	void save(Grid1D<T>& grid,const std::string& filename,const CL& closure){

		std::ofstream out(filename.c_str());

		out << "set term pngcairo size " << xRes << "," << yRes << "\n";
		out << "set output \"" << filename << ".png\"\n";
		//out << "set multiplot" << "\n";

		out << "set title '" << grid.description << "'" << "\n";
		out << "set xrange [" << grid.minX << ":" << grid.maxX << "]" << "\n";
		if(!autoScale)out << "set yrange [" << minValue << ":" << maxValue << "]" << "\n";

		out << "plot '-' with lines  linecolor rgb \"" << lineColor << "\"" << "\n";


		grid.iterateWhole(GRID1D_ITERATOR{
		   	double x = grid.minX + (ix - grid.minIndexX)*grid.dx;
		   	out << x << " " << closure(ix) << "\n";
		});

	    out << "e" << "\n";

		out.close();

		if(saveToPNG){
			std::string gnuPlotCall("gnuplot ");
			gnuPlotCall.append(filename);
			system(gnuPlotCall.c_str());
		}

	}



	template <typename T>
	void save(Grid1D<T>& grid,const std::string& filename){
		save(grid,filename,GRID1D_CALCULATOR{
			return grid.getAt(ix);
		});
	}

};
