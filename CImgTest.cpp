
#include "GridsCommon.hpp"

using namespace cimg_library;

int main(){

	CImg<unsigned char> img(127,257,1,3);
	img(10,10,0)=200;
	img(20,20,1)=200;
	img(30,30,2)=200;

	img.save_bmp("bmp.bmp");
	img.save_png("png.png");
	img.save_jpeg("jpg.jpg",80);

	return 0;
}
