
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>

int main(int, char *[]) {

	// Create an image data
	vtkSmartPointer <vtkImageData> imageData	= vtkSmartPointer<vtkImageData>::New();

	// Specify the size of the image data
	imageData->SetDimensions(5, 10, 20);
	imageData->SetNumberOfScalarComponents(1);
	imageData->SetScalarTypeToDouble();

	int* dims = imageData->GetDimensions();
	// int dims[3]; // can't do this

	std::cout << "Dims: " << " x: " << dims[0] << " y: " << dims[1] << " z: " << dims[2] << std::endl;

	std::cout << "Number of points: " << imageData->GetNumberOfPoints()	<< std::endl;
	std::cout << "Number of cells: " << imageData->GetNumberOfCells()	<< std::endl;

	//fill every entry of the image data with "2.0"
	for (int z = 0; z < dims[2]; z++)
		for (int y = 0; y < dims[1]; y++)
			for (int x = 0; x < dims[0]; x++){
				double* pixel =	static_cast<double*> (imageData->GetScalarPointer(x, y, z));
				pixel[0] = 2.0;
	}

	//retrieve the entries from the image data and print them to the screen
	for (int z = 0; z < dims[2]; z++) {
		for (int y = 0; y < dims[1]; y++) {
			for (int x = 0; x < dims[0]; x++) {
				double* pixel =	static_cast<double*> (imageData->GetScalarPointer(x, y, z));
				std::cout << pixel[0] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}



	vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();

	writer->SetDataModeToAscii();
	writer->SetCompressorTypeToNone();
	writer->SetFileName("image.vti");
	writer->SetInput(imageData);
	writer->Write();

	return 0;
}
