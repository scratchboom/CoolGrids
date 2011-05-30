#pragma once

class CsvSaver{

public:
	template <typename T>
	void save(Grid1D<T>& grid,const std::string& filename){

		std::ofstream out(filename.c_str());
		out.imbue(std::locale::classic());

		out << "x,y" << "\n";
        grid.iterateWhole(GRID1D_ITERATOR{
		   	double x = grid.minX + (ix - grid.minIndexX)*grid.dx;
		   	out << x << "," << grid.getAt(ix) << "\n";
		});

        out.close();
	}

};
