#include "GridsCommon.hpp"

int main(){

	BmpImage bmp(7,5);
	bmp.setPixel(0,0   ,1.0,0.0,0.0);
	bmp.setPixel(1,1   ,0.0,1.0,0.0);
	bmp.setPixel(2,2   ,0.0,0.0,1.0);

	bmp.save("7x5.bmp");

	return 0;
}
