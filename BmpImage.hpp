#pragma once

const int BMP_CHANNELS_PER_PIXEL=4;

class BmpImage{

	struct BMP_HEADER{
		uint32_t fileSize;
		uint32_t reserved;
		uint32_t bmpOffset;
        uint32_t infoSize;
        uint32_t width;
        uint32_t height;
	    uint16_t planesNum;
	    uint16_t bits;
	    uint32_t compression;
	    uint32_t imageSize;
	    uint32_t XPelsPerMeter;
	    uint32_t YPelsPerMeter;
	    uint32_t colorsUsed;
	    uint32_t colorsImportant;
	};

private:
	unsigned char* data;
	int width,height;
	int widthPad;//for 4 bytes allign


public:

	void setPixel(int x,int y,double R,double G,double B){
		int offset=BMP_CHANNELS_PER_PIXEL*(y*(width+widthPad)+x);
		data[offset+0]=B*255.0;
		data[offset+1]=G*255.0;
		data[offset+2]=R*255.0;
	}

	BmpImage(int width,int height){
		this->width=width;
		this->height=height;
		widthPad = 0;//(((BMP_CHANNELS_PER_PIXEL*width)/4)+1)*4 - BMP_CHANNELS_PER_PIXEL*width;
		DBGVAL(widthPad);
		data=new unsigned char[BMP_CHANNELS_PER_PIXEL*(width+widthPad)*height];
	}


	~BmpImage(){
		delete data;
	}

	void save(const std::string& filename){
		std::ofstream out(filename.c_str());

		Timer saveTimer;

		BMP_HEADER hdr;

		hdr.fileSize=2 //signature
		             +sizeof(BMP_HEADER)
		             +sizeof(char)*BMP_CHANNELS_PER_PIXEL*(width+widthPad)*height;

		hdr.reserved=0;
		hdr.bmpOffset=54;
		hdr.infoSize=40;
		hdr.width=width;
		hdr.height=height;
		hdr.planesNum=1;
		hdr.bits=32;
		hdr.compression=0;
		hdr.imageSize=sizeof(char)*BMP_CHANNELS_PER_PIXEL*(width+widthPad)*height;
		hdr.XPelsPerMeter=1000;
		hdr.YPelsPerMeter=1000;
		hdr.colorsUsed=0;
		hdr.colorsImportant=0;



		out.write("BM",2);//signature

		out.write((char*)&hdr,sizeof(hdr));
		out.write((char*)data,sizeof(char)*BMP_CHANNELS_PER_PIXEL*(width+widthPad)*height);
		out.close();

		saveTimer.logTime(filename+" saved");
	}

};
