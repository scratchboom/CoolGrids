#pragma once

class VtiSaver3D{

private:
	double minValue,maxValue;
	bool autoScale;
public:

	VtiSaver3D(){
	}



	VtiSaver3D& setValueRange(double minValue,double maxValue){
		autoScale=false;
		this->minValue=minValue;
		this->maxValue=maxValue;
		return *this;
    }

	VtiSaver3D& setAutoScale(double minValue,double maxValue){
		autoScale=true;
		return *this;
	}

	template <typename T>
	void save(Grid3D<T>& grid,const std::string& filename){

		Timer saveTimer;

		vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();

		imageData->SetDimensions(grid.nodesCountX, grid.nodesCountY,grid.nodesCountZ);
		imageData->SetNumberOfScalarComponents(1);
		imageData->SetScalarTypeToDouble();

		int* dims = imageData->GetDimensions();
		std::cout << "Dims: " << " x: " << dims[0] << " y: " << dims[1] << " z: " << dims[2] << std::endl;
		std::cout << "nodesCount: " << " x: " << grid.nodesCountX << " y: " << grid.nodesCountY << " z: " << grid.nodesCountZ << std::endl;

		std::cout << "Number of points: " << imageData->GetNumberOfPoints()	<< std::endl;
		std::cout << "Number of cells: " << imageData->GetNumberOfCells()	<< std::endl;
		std::cout << "Number of nodes: " << grid.nodesCountX*grid.nodesCountY*grid.nodesCountZ	<< std::endl;


		int numNodesIterated=0;
		for (int iz = 0; iz < grid.nodesCountZ; iz++)
			for (int iy = 0; iy < grid.nodesCountY; iy++)
				for (int ix = 0; ix < grid.nodesCountX; ix++) {
					double* pixel =	static_cast<double*> (imageData->GetScalarPointer(ix, iy, iz));
					pixel[0] = grid.getAt(grid.minIndexX+ix,grid.minIndexY+iy,grid.minIndexZ+iz);
					numNodesIterated++;
				}

		std::cout << "Number of nodes inserted: " << numNodesIterated << std::endl;


		vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();

		writer->SetDataModeToAppended();
		writer->SetCompressorTypeToZLib();
		writer->SetFileName(filename.c_str());
		writer->SetInput(imageData);
		writer->Write();

		saveTimer.logTime(filename+" saved");
	}


	template <typename T,typename CL>
	void save(Grid3D<T>& grid,const std::string& filename,const CL& closure){

		Timer saveTimer;

		vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();

		imageData->SetDimensions(grid.nodesCountX, grid.nodesCountY,grid.nodesCountZ);
		imageData->SetNumberOfScalarComponents(1);
		imageData->SetScalarTypeToDouble();

		int* dims = imageData->GetDimensions();
		std::cout << "Dims: " << " x: " << dims[0] << " y: " << dims[1] << " z: " << dims[2] << std::endl;
		std::cout << "nodesCount: " << " x: " << grid.nodesCountX << " y: " << grid.nodesCountY << " z: " << grid.nodesCountZ << std::endl;

		std::cout << "Number of points: " << imageData->GetNumberOfPoints()	<< std::endl;
		std::cout << "Number of cells: " << imageData->GetNumberOfCells()	<< std::endl;
		std::cout << "Number of nodes: " << grid.nodesCountX*grid.nodesCountY*grid.nodesCountZ	<< std::endl;


		int numNodesIterated=0;
		for (int iz = 0; iz < grid.nodesCountZ; iz++)
			for (int iy = 0; iy < grid.nodesCountY; iy++)
				for (int ix = 0; ix < grid.nodesCountX; ix++) {
					double* pixel =	static_cast<double*> (imageData->GetScalarPointer(ix, iy, iz));
					pixel[0] = closure(grid.minIndexX+ix,grid.minIndexY+iy,grid.minIndexZ+iz);
					numNodesIterated++;
				}

		std::cout << "Number of nodes inserted: " << numNodesIterated << std::endl;


		vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();

		writer->SetDataModeToAppended();
		writer->SetCompressorTypeToZLib();
		writer->SetFileName(filename.c_str());
		writer->SetInput(imageData);
		writer->Write();

		saveTimer.logTime(filename+" saved");
	}

};
