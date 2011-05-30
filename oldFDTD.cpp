#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cstdarg>
#include <string>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <map>
#include <algorithm>
#include <omp.h>

#ifdef _WIN32
#define WIN_PLATFORM
#else
#define UNIX_PLATFORM
#endif

#ifdef WIN_PLATFORM
#include <windows.h>
#define M_PI 3.14159265358979323846
#define collapse(N)
#endif

#ifdef UNIX_PLATFORM
#include <sys/time.h>
#define collapse(N) /*collapse(N)*/
#endif

#include "CImg.h"

#include "URiemann.h"

using namespace std;
using namespace cimg_library;


#pragma region platformDependent

long getTickCount(){
#ifdef WIN_PLATFORM
	return GetTickCount();
#endif

#ifdef UNIX_PLATFORM
	struct timeval tv;
    gettimeofday(&tv,NULL);
    return (tv.tv_sec*1000+tv.tv_usec/1000);
#endif
}


void pressAnyKey(){
	system("pause");
}


string trim(const string& str)
{
  string::size_type pos1 = str.find_first_not_of(" \t\r\n");
  string::size_type pos2 = str.find_last_not_of(" \t\r\n");
  return str.substr(pos1 == string::npos ? 0 : pos1,pos2 == string::npos ? str.length() - 1 : pos2 - pos1 + 1);
}

/*
string toLowerCase(string& str){
    string res=str;
    transform(res.begin(),res.end(),res.begin(),tolower);
    return res;
}

string toUpperCase(string& str){
    string res=str;
    transform(res.begin(),res.end(),res.begin(),toupper);
    return res;
}
*/


int strToInt(const string& s){
    return atoi(s.c_str());
}

double strToDouble(const string& s){
    return atof(s.c_str());
}


bool strToBool(const string& s){
    if(s.compare("true")==0)return true;
    if(s.compare("false")==0)return false;

    throw invalid_argument("Invalid argument to convert into boolean: "+s);
}



#pragma endregion

#pragma region utils

double sqr(double x){
	return x*x;
}

double sqr(double x,double y){
	return x*x+y*y;
}

double sqr(double x,double y,double z){
	return x*x+y*y+z*z;
}

int factorial(int n){
    int res=1;
    for(int i=2;i<=n;res*=n);
    return res;
}

int combination(int n,int k){
    return factorial(n)/factorial(k)/factorial(n-k);
}


const double SQRT_LN_2=0.8325546111576977563531646;
double gauss(double x,double x0,double w){
    return exp(-sqr((x-x0)*2.0*SQRT_LN_2/w));
}

const double gaussStep(double x,double x0,double w){
    if(x<x0)return gauss(x,x0,w);
    else return 1.0;
}

double linearStep(double x, double x1, double y1, double x2, double y2) {
        if (x < x1) return y1;
        else if (x > x2) return y2;
        else return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
}



/*aaSize antialias size (best is grid cell size)*/
double lens(double x, double y, double z,
                                     double x1, double y1, double z1,
                                     double x2, double y2, double z2,
                                     double r,
                                     double ior,
                                     double aaSize
                                     ){

        double r1=sqrt(sqr(x - x1) + sqr(y - y1) + sqr(z - z1));
        double r2=sqrt(sqr(x - x2) + sqr(y - y2) + sqr(z - z2));

        if(r1>r2) return linearStep(r1,r-aaSize,ior,r+aaSize,1);
        else return linearStep(r2,r-aaSize,ior,r+aaSize,1);
}

double rnd(){
    return rand()/(double(RAND_MAX)+1);
}

double _max(double d1,double d2,double d3){
    return max(d1,max(d2,d3));
}

double _max(double d1,double d2,double d3,double d4,double d5,double d6){
    return max(_max(d1,d2,d3),_max(d4,d5,d6));
}

double maxMod(double d1,double d2,double d3){
    return _max(abs(d1),abs(d2),abs(d3));
}

double maxMod(double d1,double d2,double d3,double d4,double d5,double d6){
    return _max(abs(d1),abs(d2),abs(d3),abs(d4),abs(d5),abs(d6));
}

template<typename T>
bool isInteger(T x){
	return x==floor(x);
}

void debug(string s){

}

template<typename T>
string toString(T v){
	stringstream s;
	s << v;
	return s.str();
}

template<typename T>
string toZPadString(T v,int w){
	stringstream s;
	s << setw(w) <<setfill('0') << v;
	return s.str();
}



#pragma endregion

#pragma region memory
long bytesToKb(int bytes){
	return bytes/1024;
}
long bytesToMb(int bytes){
	return bytes/1024/1024;
}
long bytesToGb(int bytes){
	return bytes/1024/1024/1024;
}

long totalBytesAllocated=0;
#pragma endregion

#pragma region classes


class Configuration{
private:
    map<string,string> paramsMap;

public:
    Configuration(const string& cfgFileName);
    string getString(const string& key);
    int getInt(const string& key);
    double getDouble(const string& key);
    bool getBoolean(const string& key);
};

Configuration::Configuration(const string& cfgFileName){

    ifstream inpStream(cfgFileName.c_str());
    if(inpStream.fail()){
        cout << "Error opening configuration file "<<  cfgFileName << "\nExiting program..." << endl;
        pressAnyKey();
        exit(1);
    }


    const int BUFF_SIZE=1024;
    char buff[BUFF_SIZE];
    string line;


    while(!inpStream.eof()){

        inpStream.getline(buff,BUFF_SIZE);
        line=buff;

        int eqIndex=line.find_first_of("=");
        int commentIndex=line.find_first_of(";");


        if(eqIndex==string::npos)continue;

        if(commentIndex!=string::npos && commentIndex<eqIndex)continue;

        int lastIndex = commentIndex!=string::npos  ?  commentIndex :  line.length();
        //cout << "lastIndex" << lastIndex << endl;

        string key=trim(line.substr(0,eqIndex));
        string value=trim(line.substr(eqIndex+1,lastIndex-eqIndex-1));

        paramsMap[key]=value;
    }

    inpStream.close();

}

string Configuration::getString(const string& key){
    auto pos=paramsMap.find(key);

    if(pos==paramsMap.end()){
        cout << "Error: there is no key '" << key <<"' in configuration file\nExiting program..." << endl;
        pressAnyKey();
        exit(1);
    }
    else return paramsMap[key];
}

int Configuration::getInt(const string& key){
    return strToInt(getString(key));
}

double Configuration::getDouble(const string& key){
    return strToDouble(getString(key));
}


bool Configuration::getBoolean(const string& key){
    bool res;
    try{
        res=strToBool(getString(key));
    }catch(invalid_argument& ia){
        cout << ia.what() << endl;
        cout << "Exiting program..." << endl;
        pressAnyKey();
        exit(1);
    }
    return res;
}



#pragma region class Timer

class Timer{

private:

    unsigned long startTick;

public:

    Timer();
    double getTime();
    void reset();
	void logTime(const string& message);
};

Timer::Timer(){
    reset();
}

void Timer::reset(){
    startTick=getTickCount();
}

double Timer::getTime(){
    return (getTickCount()-startTick)/1000.0;
}

void Timer::logTime(const string& message){
	cout<<getTime() <<": "<<message << endl;
}
#pragma endregion

#pragma region class Grid3D
template <typename T>
class Grid3D{

private:
	T* data;
    string name;

public:

	Grid3D(const string& _name);
	~Grid3D();

	int nodesCountX;
	int nodesCountY;
	int nodesCountZ;

	int cellsCountX;
	int cellsCountY;
	int cellsCountZ;

	double minIndexX,maxIndexX;
	double minIndexY,maxIndexY;
	double minIndexZ,maxIndexZ;

    double minX,maxX;
	double minY,maxY;
	double minZ,maxZ;

	double dx,dy,dz;

	void setExtentX(double _minX,double _maxX);
	void setExtentY(double _minY,double _maxY);
	void setExtentZ(double _minZ,double _maxZ);


	void setExtentIndexX(double _minIndexX,double _maxIndexX);
	void setExtentIndexY(double _minIndexY,double _maxIndexY);
	void setExtentIndexZ(double _minIndexZ,double _maxIndexZ);


	Grid3D<T> clone();
	Grid3D<T> wrap(double n);
	void build();


    void saveVTK(const string& filename);

    template<typename CL>
    void saveVTK(const string& filename,const CL& closure);

	void saveZSlice(double iz,const string& filename);


	void checkIndices(double ix,double iy,double iz);
	T& operator()(double ix,double iy,double iz);

    string getDescription();

};

template<typename T>
Grid3D<T>::Grid3D(const string& _name){
    name=_name;
	data=NULL;
};

template<typename T>
void Grid3D<T>::setExtentX(double _minX,double _maxX){
	minX=_minX;
	maxX=_maxX;
};

template<typename T>
void Grid3D<T>::setExtentY(double _minY,double _maxY){
	minY=_minY;
	maxY=_maxY;
};

template<typename T>
void Grid3D<T>::setExtentZ(double _minZ,double _maxZ){
	minZ=_minZ;
	maxZ=_maxZ;
};



template<typename T>
void Grid3D<T>::setExtentIndexX(double _minIndexX,double _maxIndexX){
	minIndexX=_minIndexX;
	maxIndexX=_maxIndexX;
};

template<typename T>
void Grid3D<T>::setExtentIndexY(double _minIndexY,double _maxIndexY){
	minIndexY=_minIndexY;
	maxIndexY=_maxIndexY;
};

template<typename T>
void Grid3D<T>::setExtentIndexZ(double _minIndexZ,double _maxIndexZ){
	minIndexZ=_minIndexZ;
	maxIndexZ=_maxIndexZ;
};

template<typename T>
void Grid3D<T>::build(){

	//cout << "building grid "<< this <<endl;

	nodesCountX=(int)(maxIndexX-minIndexX+1);
	nodesCountY=(int)(maxIndexY-minIndexY+1);
	nodesCountZ=(int)(maxIndexZ-minIndexZ+1);

	cellsCountX=(int)(maxIndexX-minIndexX);
	cellsCountY=(int)(maxIndexY-minIndexY);
	cellsCountZ=(int)(maxIndexZ-minIndexZ);

	dx=(maxX-minX)/cellsCountX;
	dy=(maxY-minY)/cellsCountY;
	dz=(maxZ-minZ)/cellsCountZ;

	data=new T[nodesCountX*nodesCountY*nodesCountZ];
	//cout << "buildt grid data"<< data <<endl;




	for(double ix=minIndexX;ix<=maxIndexX;ix++)
		for(double iy=minIndexY;iy<=maxIndexY;iy++)
			for(double iz=minIndexZ;iz<=maxIndexZ;iz++)
				operator()(ix,iy,iz)=T();//TODO was initialization with 0


	long dataSize=sizeof(T)*nodesCountX*nodesCountY*nodesCountZ;
	totalBytesAllocated+=dataSize;
	//cout << " " << dataSize << " bytes allocated    "<< bytesToKb(dataSize) << "Kb    "<< bytesToMb(dataSize) << "Mb    "<< bytesToGb(dataSize) << "Gb        " << endl;
	cout << " " << totalBytesAllocated << " total bytes allocated    "<< bytesToKb(totalBytesAllocated) << "Kb    "<< bytesToMb(totalBytesAllocated) << "Mb    "<< bytesToGb(totalBytesAllocated) << "Gb        " << endl;
};



template<typename T>
Grid3D<T> Grid3D<T>::clone(){
	Grid3D<T> res("cloned "+name);

	res.setExtentX(minX,maxX);
	res.setExtentY(minY,maxY);
	res.setExtentZ(minZ,maxZ);

	res.setExtentIndexX(minIndexX,maxIndexX);
	res.setExtentIndexY(minIndexY,maxIndexY);
	res.setExtentIndexZ(minIndexZ,maxIndexZ);

	return res;
};

template<typename T>
Grid3D<T> Grid3D<T>::wrap(double n){
	Grid3D<T> res("wraped "+name);

	res.setExtentX(minX-n*dx,maxX+n*dx);
	res.setExtentY(minY-n*dy,maxY+n*dy);
	res.setExtentZ(minZ-n*dz,maxZ+n*dz);

	res.setExtentIndexX(minIndexX-n,maxIndexX+n);
	res.setExtentIndexY(minIndexY-n,maxIndexY+n);
	res.setExtentIndexZ(minIndexZ-n,maxIndexZ+n);

	return res;
};

template<typename T>
void Grid3D<T>::saveVTK(const string& filename){

	ofstream out(filename.c_str());

	Timer saveTimer;

	out << "# vtk DataFile Version 3.0\n";
	out << "Insert Any Information Here.\n";
	out << "ASCII\n";
	out << "DATASET STRUCTURED_GRID\n";
	out << "DIMENSIONS "<< (int)nodesCountX <<" "<< (int)nodesCountY  <<" "<< (int)nodesCountZ  <<"\n";
	out << "POINTS"<<" "<< nodesCountX*nodesCountY*nodesCountZ <<" double\n";

	for (double iz = minIndexZ; iz <= maxIndexZ; iz++)
		for (double iy = minIndexY; iy <= maxIndexY; iy++)
			for (double ix = minIndexX; ix <= maxIndexX; ix++) {
				out << minX+(ix-minIndexX)*dx << " " << minY+(iy-minIndexY)*dy << " " << minZ+(iz-minIndexZ)*dz <<"\n";
			}

	out << "POINT_DATA " << nodesCountX*nodesCountY*nodesCountZ <<"\n";
	out << "SCALARS scalar_name float 1\n";
	out << "LOOKUP_TABLE default\n";

	for (double iz = minIndexZ; iz <= maxIndexZ; iz++)
		for (double iy = minIndexY; iy <= maxIndexY; iy++)
			for (double ix = minIndexX; ix <= maxIndexX; ix++) {
				out << operator()(ix, iy, iz)<<"\n";
			}

	saveTimer.logTime(filename+" saved");
}

template<typename T>
template<typename CL>
void Grid3D<T>::saveVTK(const string& filename,const CL& closure){
    ofstream out(filename.c_str());

	Timer saveTimer;

	out << "# vtk DataFile Version 3.0\n";
	out << "Insert Any Information Here.\n";
	out << "ASCII\n";
	out << "DATASET STRUCTURED_GRID\n";
	out << "DIMENSIONS "<< (int)nodesCountX <<" "<< (int)nodesCountY  <<" "<< (int)nodesCountZ  <<"\n";
	out << "POINTS"<<" "<< nodesCountX*nodesCountY*nodesCountZ <<" double\n";

	for (double iz = minIndexZ; iz <= maxIndexZ; iz++)
		for (double iy = minIndexY; iy <= maxIndexY; iy++)
			for (double ix = minIndexX; ix <= maxIndexX; ix++) {
				out << minX+(ix-minIndexX)*dx << " " << minY+(iy-minIndexY)*dy << " " << minZ+(iz-minIndexZ)*dz <<"\n";
			}

	out << "POINT_DATA " << nodesCountX*nodesCountY*nodesCountZ <<"\n";
	out << "SCALARS scalar_name float 1\n";
	out << "LOOKUP_TABLE default\n";

	for (double iz = minIndexZ; iz <= maxIndexZ; iz++)
		for (double iy = minIndexY; iy <= maxIndexY; iy++)
			for (double ix = minIndexX; ix <= maxIndexX; ix++) {
				out << closure(ix, iy, iz)<<"\n";
			}

	saveTimer.logTime(filename+" saved");
}

template<typename T>
void Grid3D<T>::saveZSlice(double iz,const string& filename){

	Timer saveTimer;
	CImg<double> img(nodesCountX,nodesCountY,1,3);

	for (double iy = minIndexY; iy <= maxIndexY; iy++)
		for (double ix = minIndexX; ix <= maxIndexX; ix++) {
					img((int)(ix-minIndexX),(int)(iy-minIndexY),0)= ((operator()(ix,iy,iz)+1.0))/2.0*360.0;
					img((int)(ix-minIndexX),(int)(iy-minIndexY),1)= 1.0;
					img((int)(ix-minIndexX),(int)(iy-minIndexY),2)= 1.0;

		}
		img.HSVtoRGB().save_bmp(filename.c_str());

	saveTimer.logTime(filename+" saved");
}

template<typename T>
void Grid3D<T>::checkIndices(double ix,double iy,double iz){

	string errorMessage;


	if(ix<minIndexX || ix>maxIndexX || !isInteger(ix - minIndexX)) errorMessage+= "wrong index ix:"+toString(ix)+"\n";
	if(iy<minIndexY || iy>maxIndexY || !isInteger(iy - minIndexY)) errorMessage+= "wrong index iy:"+toString(iy)+"\n";
	if(iz<minIndexZ || iz>maxIndexZ || !isInteger(iz - minIndexZ)) errorMessage+= "wrong index iz:"+toString(iz)+"\n";;

	if(!errorMessage.empty()){
		cout << errorMessage << endl;

        cout << "Grid description: " << getDescription();
		pressAnyKey();
		exit(1);
	}
}

template<typename T>
T& Grid3D<T>::operator()(double ix,double iy,double iz){

	checkIndices(ix,iy,iz);

	int ixIndex = (int)(ix - minIndexX);
	int iyIndex = (int)(iy - minIndexY);
	int izIndex = (int)(iz - minIndexZ);

	return data[ixIndex * nodesCountY*nodesCountZ +
		        iyIndex * nodesCountZ +
				izIndex];
}

template<typename T>
string Grid3D<T>::getDescription(){
    stringstream ss;
    ss << "Grid3D "<< name
       <<"    <"
       << minX << ","
       << minY << ","
       << minZ << ">..<"
       << maxX << ","
       << maxY << ","
       << maxZ << ">    ["
       << minIndexX << ","
       << minIndexY << ","
       << minIndexZ << "]..["
       << maxIndexX << ","
       << maxIndexY << ","
       << maxIndexZ << "]";
    return ss.str();
}



template<typename T>
Grid3D<T>::~Grid3D(){
	cout << "Grid3D destruction" << this << endl;
	cout << "Grid3D data destruction" << data << endl;
	if(data!=NULL)delete[] data;
	data=NULL;
	cout << "Grid3D destroyed" << endl;
};



#pragma endregion

#pragma region class Grid4D
const int LAYERS_COUNT_TO_MAINTAIN=5;

template<typename T>
class Grid4D{

private:

    string name;

	Grid3D<T>* layers[LAYERS_COUNT_TO_MAINTAIN];
	double itLargestAccessable;
	void switchToNextLayer();

public:

	Grid4D(const string& _name);
	~Grid4D();

	int nodesCountX;
	int nodesCountY;
	int nodesCountZ;

	double cellsCountX;
	double cellsCountY;
	double cellsCountZ;

	double dx;
	double dy;
	double dz;
    double dt;

	double minIndexT,maxIndexT;
	double minIndexX,maxIndexX;
	double minIndexY,maxIndexY;
	double minIndexZ,maxIndexZ;

	double minT,maxT;
    double minX,maxX;
	double minY,maxY;
	double minZ,maxZ;

	void setExtentT(double _minT,double _maxT);
	void setExtentX(double _minX,double _maxX);
	void setExtentY(double _minY,double _maxY);
	void setExtentZ(double _minZ,double _maxZ);

	void setExtentIndexT(double _minIndexT,double _maxIndexT);
	void setExtentIndexX(double _minIndexX,double _maxIndexX);
	void setExtentIndexY(double _minIndexY,double _maxIndexY);
	void setExtentIndexZ(double _minIndexZ,double _maxIndexZ);

	template<typename CL>
	void iterateWhole(const CL& closure);

	template<typename CL>
	void iterateInternal(int BORDER_NX,int BORDER_NY,int BORDER_NZ,const CL& closure);

    template<typename CL>
	void iterateInternal(const CL& closure);

	template<typename CL>
	void iterateBorder(const CL& closure);

    template<typename CL>
    void iterateBorder(int BORDER_WIDTH,const CL& closure);

    template<typename CL>
	void iterateBorder(int BORDER_NX,int BORDER_NY,int BORDER_NZ,const CL& closure);

    template<typename CL>
	void iterateBorder(int BORDER_OUTER_NX,int BORDER_OUTER_NY,int BORDER_OUTER_NZ,int BORDER_INNER_NX,int BORDER_INNER_NY,int BORDER_INNER_NZ,const CL& closure);

    template<typename CL>
	void iterateBorderMinX(int BORDER_WIDTH,const CL& closure);

    template<typename CL>
	void iterateBorderMaxX(int BORDER_WIDTH,const CL& closure);

    template<typename CL>
	void iterateBorderMinY(int BORDER_WIDTH,const CL& closure);

    template<typename CL>
	void iterateBorderMaxY(int BORDER_WIDTH,const CL& closure);

    template<typename CL>
	void iterateBorderMinZ(int BORDER_WIDTH,const CL& closure);

    template<typename CL>
	void iterateBorderMaxZ(int BORDER_WIDTH,const CL& closure);

    template<typename CL>
	void iterateBorderMinX(const CL& closure);

    template<typename CL>
	void iterateBorderMaxX(const CL& closure);

    template<typename CL>
	void iterateBorderMinY(const CL& closure);

    template<typename CL>
	void iterateBorderMaxY(const CL& closure);

    template<typename CL>
	void iterateBorderMinZ(const CL& closure);

    template<typename CL>
	void iterateBorderMaxZ(const CL& closure);


	Grid4D<T> clone();
	Grid4D<T> wrap(double n);
	void build();

	T& operator()(double it,double ix,double iy,double iz);
	Grid3D<T>& operator[](double it);

    string getDescription();
};

template<typename T>
Grid4D<T>::Grid4D(const string& _name){
    name=_name;
	for(int i=0;i<LAYERS_COUNT_TO_MAINTAIN;i++)layers[i]=NULL;
}

template<typename T>
void Grid4D<T>::setExtentT(double _minT,double _maxT){
	minT=_minT;
	maxT=_maxT;
};

template<typename T>
void Grid4D<T>::setExtentX(double _minX,double _maxX){
	minX=_minX;
	maxX=_maxX;
};

template<typename T>
void Grid4D<T>::setExtentY(double _minY,double _maxY){
	minY=_minY;
	maxY=_maxY;
};

template<typename T>
void Grid4D<T>::setExtentZ(double _minZ,double _maxZ){
	minZ=_minZ;
	maxZ=_maxZ;
};

template<typename T>
void Grid4D<T>::setExtentIndexT(double _minIndexT,double _maxIndexT){
	minIndexT=_minIndexT;
	maxIndexT=_maxIndexT;
};

template<typename T>
void Grid4D<T>::setExtentIndexX(double _minIndexX,double _maxIndexX){
	minIndexX=_minIndexX;
	maxIndexX=_maxIndexX;
};

template<typename T>
void Grid4D<T>::setExtentIndexY(double _minIndexY,double _maxIndexY){
	minIndexY=_minIndexY;
	maxIndexY=_maxIndexY;
};

template<typename T>
void Grid4D<T>::setExtentIndexZ(double _minIndexZ,double _maxIndexZ){
	minIndexZ=_minIndexZ;
	maxIndexZ=_maxIndexZ;
};

template<typename T>
void Grid4D<T>::switchToNextLayer(){
	//cout<< "switchToNextLayer th:" << omp_get_thread_num() << endl;
	itLargestAccessable++;
	Grid3D<T>* tmp=layers[LAYERS_COUNT_TO_MAINTAIN-1];
	for(int i=LAYERS_COUNT_TO_MAINTAIN-1;i>0;i--)layers[i]=layers[i-1];
	layers[0]=tmp;
}

template<typename T>
Grid3D<T>& Grid4D<T>::operator[](double it){

    if(!isInteger(it-minIndexT)){
        cout << "wrong time index: " << it << endl;
        cout << getDescription() << endl;
        cout << "Exiting program..." << endl;
        pressAnyKey();
        exit(1);
    }

	if(it>itLargestAccessable)switchToNextLayer();
		if((itLargestAccessable-it+1)>LAYERS_COUNT_TO_MAINTAIN || (itLargestAccessable-it<0)){
			cout << "ERROR: cannot access this layer " << it << endl;
            cout << getDescription();
			pressAnyKey();
			exit(1);
	}

	return *layers[static_cast<int>(itLargestAccessable-it)];
}

template<typename T>
T& Grid4D<T>::operator()(double it,double ix,double iy,double iz){
	return operator[](it)(ix,iy,iz);
}



template<typename T>
template<typename CL>
void Grid4D<T>::iterateWhole(const CL& closure){

    #pragma omp parallel for collapse(3)
    for(int nix=0;nix<nodesCountX;nix++)
		for(int niy=0;niy<nodesCountY;niy++)
			for(int niz=0;niz<nodesCountZ;niz++){
				double ix=minIndexX+nix;
				double iy=minIndexY+niy;
				double iz=minIndexZ+niz;

				closure(ix,iy,iz);
		}
}


template<typename T>
template<typename CL>
void Grid4D<T>::iterateInternal(const CL& closure){

	#pragma omp parallel for collapse(3)
	for(int nix=1;nix<nodesCountX-1;nix++)
		for(int niy=1;niy<nodesCountY-1;niy++)
			for(int niz=1;niz<nodesCountZ-1;niz++){
				double ix=minIndexX+nix;
				double iy=minIndexY+niy;
				double iz=minIndexZ+niz;

				closure(ix,iy,iz);
		}
}

template<typename T>
template<typename CL>
void Grid4D<T>::iterateInternal(int BORDER_NX,int BORDER_NY,int BORDER_NZ,const CL& closure){
    #pragma omp parallel for collapse(3)
	for(int nix=BORDER_NX;nix<nodesCountX-BORDER_NX;nix++)
		for(int niy=BORDER_NY;niy<nodesCountY-BORDER_NY;niy++)
			for(int niz=BORDER_NZ;niz<nodesCountZ-BORDER_NZ;niz++){
				double ix=minIndexX+nix;
				double iy=minIndexY+niy;
				double iz=minIndexZ+niz;

				closure(ix,iy,iz);
		}
}

template<typename T>
template<typename CL>
void Grid4D<T>::iterateBorder(const CL& closure){

	#pragma omp parallel for collapse(2)
	for(int nix=0;nix<nodesCountX;nix++)
		for(int niy=0;niy<nodesCountY;niy++){
			double ix=minIndexX+nix;
			double iy=minIndexY+niy;

			closure(ix,iy,minIndexZ);
			closure(ix,iy,maxIndexZ);
		}


	#pragma omp parallel for collapse(2)
    for(int nix=0;nix<nodesCountX;nix++)
		for(int niz=0;niz<nodesCountZ;niz++){
			double ix=minIndexX+nix;
			double iz=minIndexZ+niz;

			closure(ix,minIndexY,iz);
			closure(ix,maxIndexY,iz);
		}

    #pragma omp parallel for collapse(2)
    for(int niy=0;niy<nodesCountY;niy++)
		for(int niz=0;niz<nodesCountZ;niz++){
			double iy=minIndexY+niy;
			double iz=minIndexZ+niz;

			closure(minIndexX,iy,iz);
			closure(maxIndexX,iy,iz);
		}
}


template<typename T>
template<typename CL>
void Grid4D<T>::iterateBorder(int BORDER_WIDTH,const CL& closure){


	//region1
	#pragma omp parallel for collapse(3)
	for(int nix=0;nix<nodesCountX;nix++)
		for(int niy=0;niy<BORDER_WIDTH;niy++)
		   for(int niz=0;niz<nodesCountZ;niz++){

			double ix=minIndexX+nix;
			double iy=minIndexY+niy;
			double iz=minIndexZ+niz;

			closure(ix,iy,iz);
		}

	//region2
	#pragma omp parallel for collapse(3)
	for(int nix=0;nix<nodesCountX;nix++)
		for(int niy=nodesCountY-BORDER_WIDTH;niy<nodesCountY;niy++)
		   for(int niz=0;niz<nodesCountZ;niz++){

			double ix=minIndexX+nix;
			double iy=minIndexY+niy;
			double iz=minIndexZ+niz;

			closure(ix,iy,iz);
		}



    //region3
	#pragma omp parallel for collapse(3)
	for(int nix=0;nix<BORDER_WIDTH;nix++)
		for(int niy=BORDER_WIDTH/*+1*/;niy<nodesCountY-BORDER_WIDTH;niy++)
		   for(int niz=0;niz<nodesCountZ;niz++){

			double ix=minIndexX+nix;
			double iy=minIndexY+niy;
			double iz=minIndexZ+niz;

			closure(ix,iy,iz);
		}

    //region4
	#pragma omp parallel for collapse(3)
	for(int nix=nodesCountX-BORDER_WIDTH;nix<nodesCountX;nix++)
		for(int niy=BORDER_WIDTH/*+1*/;niy<nodesCountY-BORDER_WIDTH;niy++)
		   for(int niz=0;niz<nodesCountZ;niz++){

			double ix=minIndexX+nix;
			double iy=minIndexY+niy;
			double iz=minIndexZ+niz;

			closure(ix,iy,iz);
		}



	//region5
	#pragma omp parallel for collapse(3)
	for(int nix=BORDER_WIDTH/*+1*/;nix<nodesCountX-BORDER_WIDTH;nix++)
		for(int niy=BORDER_WIDTH/*+1*/;niy<nodesCountY-BORDER_WIDTH;niy++)
		   for(int niz=0;niz<BORDER_WIDTH;niz++){

			double ix=minIndexX+nix;
			double iy=minIndexY+niy;
			double iz=minIndexZ+niz;

			closure(ix,iy,iz);
		}


	//region6
	#pragma omp parallel for collapse(3)
	for(int nix=BORDER_WIDTH/*+1*/;nix<nodesCountX-BORDER_WIDTH;nix++)
		for(int niy=BORDER_WIDTH/*+1*/;niy<nodesCountY-BORDER_WIDTH;niy++)
		   for(int niz=nodesCountZ-BORDER_WIDTH;niz<nodesCountZ;niz++){

			double ix=minIndexX+nix;
			double iy=minIndexY+niy;
			double iz=minIndexZ+niz;

			closure(ix,iy,iz);
		}

}

template<typename T>
template<typename CL>
void Grid4D<T>::iterateBorder(int BORDER_NX,int BORDER_NY,int BORDER_NZ,const CL& closure){


	//region1
	#pragma omp parallel for collapse(3)
	for(int nix=0;nix<nodesCountX;nix++)
		for(int niy=0;niy<BORDER_NY;niy++)
		   for(int niz=0;niz<nodesCountZ;niz++){

			double ix=minIndexX+nix;
			double iy=minIndexY+niy;
			double iz=minIndexZ+niz;

			closure(ix,iy,iz);
		}

	//region2
	#pragma omp parallel for collapse(3)
	for(int nix=0;nix<nodesCountX;nix++)
		for(int niy=nodesCountY-BORDER_NY;niy<nodesCountY;niy++)
		   for(int niz=0;niz<nodesCountZ;niz++){

			double ix=minIndexX+nix;
			double iy=minIndexY+niy;
			double iz=minIndexZ+niz;

			closure(ix,iy,iz);
		}



    //region3
	#pragma omp parallel for collapse(3)
	for(int nix=0;nix<BORDER_NX;nix++)
		for(int niy=BORDER_NY;niy<nodesCountY-BORDER_NY;niy++)
		   for(int niz=0;niz<nodesCountZ;niz++){

			double ix=minIndexX+nix;
			double iy=minIndexY+niy;
			double iz=minIndexZ+niz;

			closure(ix,iy,iz);
		}

    //region4
	#pragma omp parallel for collapse(3)
	for(int nix=nodesCountX-BORDER_NX;nix<nodesCountX;nix++)
		for(int niy=BORDER_NY;niy<nodesCountY-BORDER_NY;niy++)
		   for(int niz=0;niz<nodesCountZ;niz++){

			double ix=minIndexX+nix;
			double iy=minIndexY+niy;
			double iz=minIndexZ+niz;

			closure(ix,iy,iz);
		}



	//region5
	#pragma omp parallel for collapse(3)
	for(int nix=BORDER_NX;nix<nodesCountX-BORDER_NX;nix++)
		for(int niy=BORDER_NY;niy<nodesCountY-BORDER_NY;niy++)
		   for(int niz=0;niz<BORDER_NZ;niz++){

			double ix=minIndexX+nix;
			double iy=minIndexY+niy;
			double iz=minIndexZ+niz;

			closure(ix,iy,iz);
		}


	//region6
	#pragma omp parallel for collapse(3)
	for(int nix=BORDER_NX;nix<nodesCountX-BORDER_NX;nix++)
		for(int niy=BORDER_NY;niy<nodesCountY-BORDER_NY;niy++)
		   for(int niz=nodesCountZ-BORDER_NZ;niz<nodesCountZ;niz++){

			double ix=minIndexX+nix;
			double iy=minIndexY+niy;
			double iz=minIndexZ+niz;

			closure(ix,iy,iz);
		}

}



template<typename T>
template<typename CL>
void Grid4D<T>::iterateBorder(int BORDER_OUTER_NX,int BORDER_OUTER_NY,int BORDER_OUTER_NZ,int BORDER_INNER_NX,int BORDER_INNER_NY,int BORDER_INNER_NZ,const CL& closure){


	//region1
	#pragma omp parallel for collapse(3)
	for(int nix=BORDER_OUTER_NX;nix<nodesCountX-BORDER_OUTER_NX;nix++)
		for(int niy=BORDER_OUTER_NY;niy<BORDER_INNER_NY;niy++)
		   for(int niz=BORDER_OUTER_NZ;niz<nodesCountZ-BORDER_OUTER_NZ;niz++){

			double ix=minIndexX+nix;
			double iy=minIndexY+niy;
			double iz=minIndexZ+niz;

			closure(ix,iy,iz);
		}

	//region2
	#pragma omp parallel for collapse(3)
	for(int nix=BORDER_OUTER_NX;nix<nodesCountX-BORDER_OUTER_NX;nix++)
		for(int niy=nodesCountY-BORDER_INNER_NY;niy<nodesCountY-BORDER_OUTER_NY;niy++)
		   for(int niz=BORDER_OUTER_NZ;niz<nodesCountZ-BORDER_OUTER_NZ;niz++){

			double ix=minIndexX+nix;
			double iy=minIndexY+niy;
			double iz=minIndexZ+niz;

			closure(ix,iy,iz);
		}



    //region3
	#pragma omp parallel for collapse(3)
	for(int nix=BORDER_OUTER_NX;nix<BORDER_INNER_NX;nix++)
		for(int niy=BORDER_INNER_NY;niy<nodesCountY-BORDER_INNER_NY;niy++)
		   for(int niz=BORDER_OUTER_NZ;niz<nodesCountZ-BORDER_OUTER_NZ;niz++){

			double ix=minIndexX+nix;
			double iy=minIndexY+niy;
			double iz=minIndexZ+niz;

			closure(ix,iy,iz);
		}

    //region4
	#pragma omp parallel for collapse(3)
	for(int nix=nodesCountX-BORDER_INNER_NX;nix<nodesCountX-BORDER_OUTER_NX;nix++)
		for(int niy=BORDER_INNER_NY;niy<nodesCountY-BORDER_INNER_NY;niy++)
		   for(int niz=BORDER_OUTER_NZ;niz<nodesCountZ-BORDER_OUTER_NZ;niz++){

			double ix=minIndexX+nix;
			double iy=minIndexY+niy;
			double iz=minIndexZ+niz;

			closure(ix,iy,iz);
		}



	//region5
	#pragma omp parallel for collapse(3)
	for(int nix=BORDER_INNER_NX;nix<nodesCountX-BORDER_INNER_NX;nix++)
		for(int niy=BORDER_INNER_NY;niy<nodesCountY-BORDER_INNER_NY;niy++)
		   for(int niz=BORDER_OUTER_NZ;niz<BORDER_INNER_NZ;niz++){

			double ix=minIndexX+nix;
			double iy=minIndexY+niy;
			double iz=minIndexZ+niz;

			closure(ix,iy,iz);
		}


	//region6
	#pragma omp parallel for collapse(3)
	for(int nix=BORDER_INNER_NX;nix<nodesCountX-BORDER_INNER_NX;nix++)
		for(int niy=BORDER_INNER_NY;niy<nodesCountY-BORDER_INNER_NY;niy++)
		   for(int niz=nodesCountZ-BORDER_INNER_NZ;niz<nodesCountZ-BORDER_OUTER_NZ;niz++){

			double ix=minIndexX+nix;
			double iy=minIndexY+niy;
			double iz=minIndexZ+niz;

			closure(ix,iy,iz);
		}

}



template<typename T>
template<typename CL>
void Grid4D<T>::iterateBorderMinX(int BORDER_WIDTH,const CL& closure){

    for(int nix=0;nix<BORDER_WIDTH;nix++)
        for(int niy=0;niy<nodesCountY;niy++)
		   for(int niz=0;niz<nodesCountZ;niz++){

			double ix=minIndexX+nix;
			double iy=minIndexY+niy;
			double iz=minIndexZ+niz;

			closure(ix,iy,iz);
		}
}


template<typename T>
template<typename CL>
void Grid4D<T>::iterateBorderMaxX(int BORDER_WIDTH,const CL& closure){

    for(int nix=nodesCountX-BORDER_WIDTH;nix<nodesCountX;nix++)
        for(int niy=0;niy<nodesCountY;niy++)
		   for(int niz=0;niz<nodesCountZ;niz++){

			double ix=minIndexX+nix;
			double iy=minIndexY+niy;
			double iz=minIndexZ+niz;

			closure(ix,iy,iz);
		}
}


template<typename T>
template<typename CL>
void Grid4D<T>::iterateBorderMinY(int BORDER_WIDTH,const CL& closure){

    for(int nix=0;nix<nodesCountX;nix++)
        for(int niy=0;niy<BORDER_WIDTH;niy++)
		   for(int niz=0;niz<nodesCountZ;niz++){

			double ix=minIndexX+nix;
			double iy=minIndexY+niy;
			double iz=minIndexZ+niz;

			closure(ix,iy,iz);
		}
}


template<typename T>
template<typename CL>
void Grid4D<T>::iterateBorderMaxY(int BORDER_WIDTH,const CL& closure){

    for(int nix=0;nix<nodesCountX;nix++)
        for(int niy=nodesCountY-BORDER_WIDTH;niy<nodesCountY;niy++)
		   for(int niz=0;niz<nodesCountZ;niz++){

			double ix=minIndexX+nix;
			double iy=minIndexY+niy;
			double iz=minIndexZ+niz;

			closure(ix,iy,iz);
		}
}



template<typename T>
template<typename CL>
void Grid4D<T>::iterateBorderMinZ(int BORDER_WIDTH,const CL& closure){

    for(int nix=0;nix<nodesCountX;nix++)
        for(int niy=0;niy<nodesCountY;niy++)
		   for(int niz=0;niz<BORDER_WIDTH;niz++){

			double ix=minIndexX+nix;
			double iy=minIndexY+niy;
			double iz=minIndexZ+niz;

			closure(ix,iy,iz);
		}
}


template<typename T>
template<typename CL>
void Grid4D<T>::iterateBorderMaxZ(int BORDER_WIDTH,const CL& closure){

    for(int nix=0;nix<nodesCountX;nix++)
        for(int niy=0;niy<nodesCountY;niy++)
		   for(int niz=nodesCountZ-BORDER_WIDTH;niz<nodesCountZ;niz++){

			double ix=minIndexX+nix;
			double iy=minIndexY+niy;
			double iz=minIndexZ+niz;

			closure(ix,iy,iz);
		}
}



template<typename T>
template<typename CL>
void Grid4D<T>::iterateBorderMinX(const CL& closure){
    iterateBorderMinX(1,closure);
}


template<typename T>
template<typename CL>
void Grid4D<T>::iterateBorderMaxX(const CL& closure){
    iterateBorderMaxX(1,closure);
}


template<typename T>
template<typename CL>
void Grid4D<T>::iterateBorderMinY(const CL& closure){
    iterateBorderMinY(1,closure);
}


template<typename T>
template<typename CL>
void Grid4D<T>::iterateBorderMaxY(const CL& closure){
    iterateBorderMaxY(1,closure);
}



template<typename T>
template<typename CL>
void Grid4D<T>::iterateBorderMinZ(const CL& closure){
    iterateBorderMinZ(1,closure);
}


template<typename T>
template<typename CL>
void Grid4D<T>::iterateBorderMaxZ(const CL& closure){
    iterateBorderMaxZ(1,closure);
}







template<typename T>
Grid4D<T> Grid4D<T>::clone(){

	Grid4D<T> res;

	res.setExtentT(minT,maxT);
	res.setExtentX(minX,maxX);
	res.setExtentY(minY,maxY);
	res.setExtentZ(minZ,maxZ);

	res.setExtentIndexT(minIndexT,maxIndexT);
	res.setExtentIndexX(minIndexX,maxIndexX);
	res.setExtentIndexY(minIndexY,maxIndexY);
	res.setExtentIndexZ(minIndexZ,maxIndexZ);

	return res;
}

template<typename T>
Grid4D<T> Grid4D<T>::wrap(double n){

	Grid4D<T> res("warped "+name);

	res.setExtentT(minT,maxT);
	res.setExtentX(minX-n*dx,maxX+n*dx);
	res.setExtentY(minY-n*dy,maxY+n*dy);
	res.setExtentZ(minZ-n*dz,maxZ+n*dz);

	res.setExtentIndexT(minIndexT,maxIndexT);
	res.setExtentIndexX(minIndexX-n,maxIndexX+n);
	res.setExtentIndexY(minIndexY-n,maxIndexY+n);
	res.setExtentIndexZ(minIndexZ-n,maxIndexZ+n);

	return res;
}

template<typename T>
void Grid4D<T>::build(){

    itLargestAccessable=minIndexT;

	nodesCountX=(int)(maxIndexX-minIndexX+1);
	nodesCountY=(int)(maxIndexY-minIndexY+1);
	nodesCountZ=(int)(maxIndexZ-minIndexZ+1);

	cellsCountX=(int)(maxIndexX-minIndexX);
	cellsCountY=(int)(maxIndexY-minIndexY);
	cellsCountZ=(int)(maxIndexZ-minIndexZ);

	dx=(maxX-minX)/cellsCountX;
	dy=(maxY-minY)/cellsCountY;
	dz=(maxZ-minZ)/cellsCountZ;
    dt=(maxT-minT)/(maxIndexT-minIndexT+1);

	//cout << "building Grid4D " << this << endl;
	for(int i=0;i<LAYERS_COUNT_TO_MAINTAIN;i++){
		layers[i]=new Grid3D<T>(name);

		layers[i]->setExtentX(minX,maxX);
		layers[i]->setExtentY(minY,maxY);
		layers[i]->setExtentZ(minZ,maxZ);

		layers[i]->setExtentIndexX(minIndexX,maxIndexX);
		layers[i]->setExtentIndexY(minIndexY,maxIndexY);
		layers[i]->setExtentIndexZ(minIndexZ,maxIndexZ);

		layers[i]->build();
	}
}

template<typename T>
string Grid4D<T>::getDescription(){
    stringstream ss;
    ss << "Grid4D " << name << "  <" << minX << "," << minY << "," << minZ << ">..<" << maxX << "," << maxY << "," << maxZ <<     ">    [" << minIndexX << "," << minIndexY << "," << minIndexZ << "]..[" << maxIndexX << "," << maxIndexY << "," << maxIndexZ << "]";
    return ss.str();
}

template<typename T>
Grid4D<T>::~Grid4D(){
	cout << "   Grid4D destruction this:" << this <<" thread" << omp_get_thread_num() << endl;
	for(int i=0;i<LAYERS_COUNT_TO_MAINTAIN;i++)if(layers[i]!=NULL)delete layers[i];
}

#pragma endregion

#pragma region iterators
#define GRID_ITERATOR [&](double ix,double iy,double iz)
#define GRID_CALCULATOR [&](double ix,double iy,double iz)->double
#pragma endregion


#pragma endregion


void mainFDTD();



#pragma region tests

void test1(){
    Grid4D<double> g("g");
	g.setExtentIndexT(0,100);
	g.setExtentIndexX(0,2);
	g.setExtentIndexY(0,2);
	g.setExtentIndexZ(0,2);

	g.setExtentT(0,100);
	g.setExtentX(0,16);
	g.setExtentY(0,16);
	g.setExtentZ(0,16);

	g.build();


	Grid4D<double> gg("gg"),ggg("ggg");


	gg=g.wrap(8);
	gg.build();

	ggg=g.wrap(8);
	ggg.build();


	for(double it=0;it<100;it++){
		gg.iterateWhole(GRID_ITERATOR{
			gg(0,ix,iy,iz)=5;
		});
	}

	pressAnyKey();
}

template <typename T>
class C{
private:
	T data;
public:
	template<typename CL>
	void function(const CL& closure);
};

template<typename T>
template<typename CL>
void C<T>::function(const CL& closure){
	closure(1,2,3);
}

void testClosure(){
	C<int> c;
	c.function([&](int x,int y,int z){cout << x << y << z ;});

	pressAnyKey();
}

void testIterateBorder(){

	Grid4D<double> grid("grid");
	grid.setExtentIndexT(0,10);
	grid.setExtentIndexX(0,20);
	grid.setExtentIndexY(0,30);
	grid.setExtentIndexZ(0,40);

	grid.setExtentT(0,10);
	grid.setExtentX(0,20);
	grid.setExtentY(0,20);
	grid.setExtentZ(0,20);

	grid.build();

	grid[0];

	grid.iterateWhole(GRID_ITERATOR{
		grid(0,ix,iy,iz)=1.0;
	});


	grid.iterateBorder(8,GRID_ITERATOR{
		grid(0,ix,iy,iz)=grid(0,ix,iy,iz)+1.0;
	});


	grid[0].saveVTK("c:\\results\\testIterateBorder.vtk");

	pressAnyKey();
}

void testIterateBorder2(){

	Grid4D<double> grid("grid");
	grid.setExtentIndexT(0,10);
	grid.setExtentIndexX(0,20);
	grid.setExtentIndexY(0,30);
	grid.setExtentIndexZ(0,40);

	grid.setExtentT(0,10);
	grid.setExtentX(0,20);
	grid.setExtentY(0,20);
	grid.setExtentZ(0,20);

	grid.build();

	grid[0];

	grid.iterateWhole(GRID_ITERATOR{
		grid(0,ix,iy,iz)=1.0;
	});


	grid.iterateBorder(0,1,1,GRID_ITERATOR{
		grid(0,ix,iy,iz)=grid(0,ix,iy,iz)+1.0;
	});


	grid[0].saveVTK("c:\\results\\testIterateBorder2.vtk");

	pressAnyKey();
}

void testIterateBorder3(){

	Grid4D<double> grid("grid");
	grid.setExtentIndexT(0,10);
	grid.setExtentIndexX(0,20);
	grid.setExtentIndexY(0,30);
	grid.setExtentIndexZ(0,40);

	grid.setExtentT(0,10);
	grid.setExtentX(0,20);
	grid.setExtentY(0,20);
	grid.setExtentZ(0,20);

	grid.build();

	grid[0];

	grid.iterateWhole(GRID_ITERATOR{
		grid(0,ix,iy,iz)=1.0;
	});


	grid.iterateBorder(1,1,1,
                       3,3,3,GRID_ITERATOR{
		grid(0,ix,iy,iz)=grid(0,ix,iy,iz)+1.0;
	});


	grid[0].saveVTK("c:\\results\\testIterateBorder3.vtk");

	pressAnyKey();
}

void testFibonacci(){
    /*
    Grid4D<double> fib;
	fib.setExtentIndexT(0,10);
	fib.setExtentIndexX(0,2);
	fib.setExtentIndexY(0,2);
	fib.setExtentIndexZ(0,2);

	fib.setExtentT(0,10);
	fib.setExtentX(0,20);
	fib.setExtentY(0,20);
	fib.setExtentZ(0,20);

	fib.build();

	fib[0].iterateInternal(GRID_ITERATOR{cout << ix << " " << iy << " " << iz << endl;});


	fib[0].iterateWhole(GRID_ITERATOR{
		cout << ix << " " << iy << " " << iz << endl;
		fib[0](ix,iy,iz)=0;
	});



	fib[1].iterateWhole(GRID_ITERATOR{
		cout << ix << " " << iy << " " << iz << endl;
		fib[1](ix,iy,iz)=1;
	});




	for(int it=2;it<10;it++){

		fib[it].iterateWhole(GRID_ITERATOR{
			cout << ix << " " << iy << " " << iz << endl;
			fib[it](ix,iy,iz)=fib[it-1](ix,iy,iz)+fib[it-2](ix,iy,iz);
		});


		cout << fib[it](1,1,1) <<endl;
	}
    */
}

void testConfiguration(){

    cout << "testConfiguration() start" << endl;

    Configuration cfg("C:\\diplom\\conf.ini");

    double C=cfg.getDouble("C");
    double DT=cfg.getDouble("DT");

    int Nx=cfg.getInt("Nx");
    int Ny=cfg.getInt("Ny");
    int Nz=cfg.getInt("Nz");

    string saveDir=cfg.getString("saveDir");


    cout << "C: " << C << endl;
    cout << "DT:" << DT << endl;

    cout << "Nx: " << Nx << endl;
    cout << "Ny: " << Ny << endl;
    cout << "Nz: " << Nz << endl;

    cout << "saveDir: "<< saveDir << endl;

    cout << "testConfiguration() end" << endl;

    pressAnyKey();
}

#pragma endregion

//global variables to be filled from conf.ini file
int SAVE_EVERY_N;
int THREAD_NUM;
string SAVE_DIR;
bool saveBMP;
bool saveVTK;

void mainFDTD(){

	double C=1.0;

	double Nx=64+1;
	double Ny=64+1;
	double Nz=64;
	double Nt=5000;


	double Sx=1.0;
	double Sy= (Sx/Nx) * Ny;
	double Sz= (Sx/Nx) * Nz;
	double T=1.0;

	double dx=Sx/Nx;
	double dy=Sy/Ny;
	double dz=Sz/Nz;
	double DT=T/Nt;

    double IMPULSE_TIME_WIDTH=0.5*0.25*Sx/C;
    double IMPULSE_PERIOD=IMPULSE_TIME_WIDTH*4.0;
    double IMPULSE_WAVE_LENGTH=IMPULSE_PERIOD*C;

    DT=0.5*dx/C;


    cout << "stability parameters:" << endl;
    cout << "WaveLength/dx = " << IMPULSE_WAVE_LENGTH/dx << "  (should be >32)" << endl;
    cout << "C*dt/dx = " << C*DT/dx << "  (should be <1/sqrt(3) ~ <0.57)" << endl;
    cout << "WavePeriod/dt = " << IMPULSE_PERIOD/DT << "  (should be >32)" << endl;

    pressAnyKey();


	int PML_N=8;
	double PML_WIDTH=PML_N*dx;


	Grid4D<double> Ex("Ex"),
                   Ey("Ey"),
                   Ez("Ez");
	Grid4D<double> Hx("Hx"),
                   Hy("Hy"),
                   Hz("Hz");

	/*
	Grid4D<double> Ex_y("Ex_y"),Ex_z("Ex_z");
	Grid4D<double> Ey_x("Ey_x"),Ey_z("Ey_x");
	Grid4D<double> Ez_x("Ez_x"),Ez_y("Ez_y");

	Grid4D<double> Hx_y("Hx_y"),Hx_z("Hx_z");
	Grid4D<double> Hy_x("Hy_x"),Hy_z("Hy_z");
	Grid4D<double> Hz_x("Hz_x"),Hz_y("Hz_y");
	*/

    Grid4D<double> epsilon("epsilon");
	Grid4D<double> mu("mu");



#pragma region buildGrids

#pragma region buildEx
	Ex.setExtentT(0,T);
	Ex.setExtentX(0.5*dx , Sx-0.5*dx);
	Ex.setExtentY(0,Sy);
	Ex.setExtentZ(0,Sz);

	Ex.setExtentIndexT(0,Nt);
	Ex.setExtentIndexX(0.5 , Nx-0.5);
	Ex.setExtentIndexY(0,Ny);
	Ex.setExtentIndexZ(0,Nz);

	Ex.build();
#pragma endregion

#pragma region buildEy
	Ey.setExtentT(0,T);
	Ey.setExtentX(0 , Sx);
	Ey.setExtentY(0.5*dy , Sy-0.5*dy);
	Ey.setExtentZ(0 , Sz);

	Ey.setExtentIndexT(0,Nt);
	Ey.setExtentIndexX(0,Nx);
	Ey.setExtentIndexY(0.5,Ny-0.5);
    Ey.setExtentIndexZ(0,Nz);

	Ey.build();
#pragma endregion

#pragma region buildEz
	Ez.setExtentT(0,T);
	Ez.setExtentX(0 , Sx);
	Ez.setExtentY(0 , Sy);
	Ez.setExtentZ(0.5*dz , Sz-0.5*dz);

	Ez.setExtentIndexT(0,Nt);
	Ez.setExtentIndexX(0 , Nx);
	Ez.setExtentIndexY(0 , Ny);
	Ez.setExtentIndexZ(0.5 , Nz-0.5);

	Ez.build();
#pragma endregion

#pragma region buildEpsilon
	epsilon.setExtentT(0,T);
	epsilon.setExtentX(0,Sx);
	epsilon.setExtentY(0,Sy);
	epsilon.setExtentZ(0,Sz);

	epsilon.setExtentIndexT(0,Nt);
	epsilon.setExtentIndexX(0,Nx);
	epsilon.setExtentIndexY(0,Ny);
	epsilon.setExtentIndexZ(0,Nz);

	epsilon.build();

    epsilon.iterateWhole(GRID_ITERATOR{
        epsilon(0,ix,iy,iz)=lens(ix*dx,iy*dy,iz*dz,
                                 Sx*0.25+0.3*Sx,0.5*Sy,0.5*Sz,
                                 Sx*0.75+0.3*Sx,0.5*Sy,0.5*Sz,
                                 0.3*Sx,
                                 1.5,
                                 dx
                                 );
        epsilon(0,ix,iy,iz) = 1.0;
    });
#pragma endregion


#pragma region buildMu
	mu.setExtentT(0,T);
	mu.setExtentX(-0.5*dx,Sx+0.5*dx);
	mu.setExtentY(-0.5*dy,Sy+0.5*dy);
	mu.setExtentZ(-0.5*dz,Sz+0.5*dz);

	mu.setExtentIndexT(0,Nt);
	mu.setExtentIndexX(-0.5,Nx+0.5);
	mu.setExtentIndexY(-0.5,Ny+0.5);
	mu.setExtentIndexZ(-0.5,Nz+0.5);

	mu.build();

    mu.iterateWhole(GRID_ITERATOR{
        mu(0,ix,iy,iz)=lens(ix*dx,iy*dy,iz*dz,
                                 Sx*0.25+0.3*Sx,0.5*Sy,0.5*Sz,
                                 Sx*0.75+0.3*Sx,0.5*Sy,0.5*Sz,
                                 0.32*Sx,
                                 1.5,
                                 dx
                                 );
        mu(0,ix,iy,iz)=1.0;
    });
#pragma endregion



#pragma region buildHx
	Hx.setExtentT(0,T);
	Hx.setExtentX(0 , Sx);
	Hx.setExtentY(0.5*dy , Sy-0.5*dy);
	Hx.setExtentZ(0.5*dz , Sz-0.5*dz);

	Hx.setExtentIndexT(0.5,Nt+0.5);
	Hx.setExtentIndexX(0   , Nx);
	Hx.setExtentIndexY(0.5 , Ny-0.5);
	Hx.setExtentIndexZ(0.5 , Nz-0.5);

	Hx.build();
#pragma endregion

#pragma region buildHy
	Hy.setExtentT(0,T);
	Hy.setExtentX(0.5*dx , Sx-0.5*dx);
	Hy.setExtentY(0      , Sy       );
	Hy.setExtentZ(0.5*dz , Sz-0.5*dz);

	Hy.setExtentIndexT(0.5,Nt+0.5);
	Hy.setExtentIndexX(0.5 , Nx-0.5);
	Hy.setExtentIndexY(0   , Ny    );
    Hy.setExtentIndexZ(0.5 , Nz-0.5);

	Hy.build();
#pragma endregion

#pragma region buildHz
	Hz.setExtentT(0,T);
	Hz.setExtentX(0.5*dx , Sx-0.5*dx);
	Hz.setExtentY(0.5*dy , Sy-0.5*dy);
	Hz.setExtentZ(0      , Sz);

	Hz.setExtentIndexT(0.5,Nt+0.5);
	Hz.setExtentIndexX(0.5 , Nx-0.5);
	Hz.setExtentIndexY(0.5 , Ny-0.5);
	Hz.setExtentIndexZ(0   , Nz    );

	Hz.build();
#pragma endregion

#pragma endregion


#pragma region buildSplitGrids
    /*
	Ex_y = Ex.wrap(PML_N);
	Ex_y.build();

	Ex_z = Ex.wrap(PML_N);
	Ex_z.build();



	Ey_x=Ey.wrap(PML_N);
	Ey_x.build();

	Ey_z=Ey.wrap(PML_N);
	Ey_z.build();



	Ez_x=Ez.wrap(PML_N);
	Ez_x.build();

	Ez_y=Ez.wrap(PML_N);
	Ez_y.build();





	Hx_y = Hx.wrap(PML_N);
	Hx_y.build();

	Hx_z = Hx.wrap(PML_N);
	Hx_z.build();



	Hy_x=Hy.wrap(PML_N);
	Hy_x.build();

	Hy_z=Hy.wrap(PML_N);
	Hy_z.build();



	Hz_x=Hz.wrap(PML_N);
	Hz_x.build();

	Hz_y=Hz.wrap(PML_N);
	Hz_y.build();
    */
#pragma endregion


	Timer mainTimer;
	for(double it=0;it<=Nt;it++){
		mainTimer.logTime("calculating layer "+ toString(it));


		Ex[it+1];
		Ey[it+1];
		Ez[it+1];
		Hx[it+1.5];
		Hy[it+1.5];
		Hz[it+1.5];

        /*
        Ex_y[it+1];
        Ex_z[it+1];

		Ey_x[it+1];
        Ey_z[it+1];

		Ez_x[it+1];
        Ez_y[it+1];

		Hx_y[it+1.5];
        Hx_z[it+1.5];

		Hy_x[it+1.5];
        Hy_z[it+1.5];

		Hz_x[it+1.5];
        Hz_y[it+1.5];
        */


        cout << "source\n" << endl;
		double t = it * DT;
        double sourceHz=gauss(t,IMPULSE_TIME_WIDTH*2,IMPULSE_TIME_WIDTH);
        Hz(it+0.5, Nx/2.0 , Ny/2.0 , Nz/2.0) = sourceHz;
        cout << "sourceHz=" << sourceHz << endl;


        cout << "Ex" << endl;
        Ex.iterateInternal(0,1,1,GRID_ITERATOR{
            double EPS = 0.5*(epsilon(0,ix-0.5,iy,iz)+epsilon(0,ix+0.5,iy,iz));

            Ex(it + 1, ix, iy, iz) =
                Ex(it,ix,iy,iz) + DT*C/EPS * (
                (Hz(it+0.5, ix, iy+0.5, iz) - Hz(it+0.5, ix, iy - 0.5, iz)) / dy
                -
                (Hy(it+0.5, ix, iy, iz+0.5) - Hy(it+0.5, ix, iy, iz - 0.5)) / dz
                );
        });

        cout << "Ey" << endl;
        Ey.iterateInternal(1,0,1,GRID_ITERATOR{
            double EPS = 0.5*(epsilon(0,ix,iy-0.5,iz)+epsilon(0,ix,iy+0.5,iz));

            Ey(it+1, ix, iy, iz)=
                Ey(it, ix, iy, iz) +
                DT*C/EPS * (
                (Hx(it+0.5, ix, iy, iz+0.5) - Hx(it+0.5, ix, iy, iz-0.5)) / dz
                -
                (Hz(it+0.5, ix+0.5, iy, iz) - Hz(it+0.5, ix-0.5, iy, iz)) / dx
                );
        });

        cout << "Ez" << endl;
        Ez.iterateInternal(1,1,0,GRID_ITERATOR{
            double EPS = 0.5*(epsilon(0,ix,iy,iz-0.5)+epsilon(0,ix,iy,iz+0.5));

            Ez(it + 1, ix, iy, iz)=
                Ez(it, ix, iy, iz) +
                DT*C/EPS * (
                (Hy(it+0.5, ix+0.5, iy, iz) - Hy(it+0.5, ix-0.5, iy, iz)) / dx
                -
                (Hx(it+0.5, ix, iy+0.5, iz) - Hx(it+0.5, ix, iy-0.5, iz)) / dy
                );
        });



        /*PEC for main region*/
        /*
        Ex.iterateBorder(0,1,1,GRID_ITERATOR{
            Ex(it+1,ix,iy,iz)=0;
        });

        Ey.iterateBorder(1,0,1,GRID_ITERATOR{
            Ey(it+1,ix,iy,iz)=0;
        });

        Ez.iterateBorder(1,1,0,GRID_ITERATOR{
            Ez(it+1,ix,iy,iz)=0;
        });*/

        /*Liao 3rd order absorbing boundary conditions*/
        double a1=0.995;
        double a2=0.99;
        //cout << "Ex1" << endl;
        Ex.iterateBorderMinY(GRID_ITERATOR{
            double u1=Ex(it,ix,iy+1,iz);
            double u2=Ex(it-1,ix,iy+2,iz);
            double u3=Ex(it-2,ix,iy+3,iz);
            double u4=Ex(it-3,ix,iy+4,iz);

            double d1=u1-u2;
            double d2=u2-u3;
            double d3=u3-u4;

            double dd1=d1-d2;
            double dd2=d2-d3;

            double ddd1=dd1-dd2;

            if(it>90 && iz==32 && ix==32.5){

                cout << "u1 u2 u3 u4    d1 d2 d3   dd1 dd2  ddd1: " << u1 << " "<< u2 << " "<< u3 << " "<< u4 << "    " <<endl
                                                                    << d1 << " "<< d2 << " "<< d3 << "     "<<endl
                                                                    << dd1 << " "<<dd2  << "   "<<endl
                                                                    << ddd1 << "        " << u1+d1+dd1+ddd1 <<endl;
            }

            Ex(it+1,ix,iy,iz)=u1+d1+dd1+ddd1;
        });
        //cout << "Ex2" << endl;
        Ex.iterateBorderMaxY(GRID_ITERATOR{
            double u1=Ex(it,ix,iy-1,iz);
            double u2=Ex(it-1,ix,iy-2,iz);
            double u3=Ex(it-2,ix,iy-3,iz);
            double u4=Ex(it-3,ix,iy-4,iz);

            double d1=u1-u2;
            double d2=u2-u3;
            double d3=u3-u4;

            double dd1=d1-d2;
            double dd2=d2-d3;

            double ddd1=dd1-dd2;

            Ex(it+1,ix,iy,iz)=u1+d1+dd1+ddd1;
        });
        //cout << "Ex3" << endl;
        Ex.iterateBorderMinZ(GRID_ITERATOR{
            double u1=Ex(it,ix,iy,iz+1);
            double u2=Ex(it-1,ix,iy,iz+2);
            double u3=Ex(it-2,ix,iy,iz+3);
            double u4=Ex(it-3,ix,iy,iz+4);

            double d1=u1-u2;
            double d2=u2-u3;
            double d3=u3-u4;

            double dd1=d1-d2;
            double dd2=d2-d3;

            double ddd1=dd1-dd2;

            Ex(it+1,ix,iy,iz)=u1+d1+dd1+ddd1;
        });
        //cout << "Ex4" << endl;
        Ex.iterateBorderMaxZ(GRID_ITERATOR{
            double u1=Ex(it,ix,iy,iz-1);
            double u2=Ex(it-1,ix,iy,iz-2);
            double u3=Ex(it-2,ix,iy,iz-3);
            double u4=Ex(it-3,ix,iy,iz-4);

            double d1=u1-u2;
            double d2=u2-u3;
            double d3=u3-u4;

            double dd1=d1-d2;
            double dd2=d2-d3;

            double ddd1=dd1-dd2;

            Ex(it+1,ix,iy,iz)=u1+d1+dd1+ddd1;
        });


        //cout << "Ex5" << endl;
        Ey.iterateBorderMinX(GRID_ITERATOR{
            double u1=Ey(it,ix+1,iy,iz);
            double u2=Ey(it-1,ix+2,iy,iz);
            double u3=Ey(it-2,ix+3,iy,iz);
            double u4=Ey(it-3,ix+4,iy,iz);

            double d1=u1-u2;
            double d2=u2-u3;
            double d3=u3-u4;

            double dd1=d1-d2;
            double dd2=d2-d3;

            double ddd1=dd1-dd2;

            Ey(it+1,ix,iy,iz)=u1+d1+dd1+ddd1;
        });

        //cout << "Ex6" << endl;
        Ey.iterateBorderMaxX(GRID_ITERATOR{
            double u1=Ey(it,ix-1,iy,iz);
            double u2=Ey(it-1,ix-2,iy,iz);
            double u3=Ey(it-2,ix-3,iy,iz);
            double u4=Ey(it-3,ix-4,iy,iz);

            double d1=u1-u2;
            double d2=u2-u3;
            double d3=u3-u4;

            double dd1=d1-d2;
            double dd2=d2-d3;

            double ddd1=dd1-dd2;

            Ey(it+1,ix,iy,iz)=u1+d1+dd1+ddd1;
        });

        //cout << "Ex7" << endl;
        Ey.iterateBorderMinZ(GRID_ITERATOR{
            double u1=Ey(it,ix,iy,iz+1);
            double u2=Ey(it-1,ix,iy,iz+2);
            double u3=Ey(it-2,ix,iy,iz+3);
            double u4=Ey(it-3,ix,iy,iz+4);

            double d1=u1-u2;
            double d2=u2-u3;
            double d3=u3-u4;

            double dd1=d1-d2;
            double dd2=d2-d3;

            double ddd1=dd1-dd2;

            Ey(it+1,ix,iy,iz)=u1+d1+dd1+ddd1;
        });

        //cout << "Ex8" << endl;
        Ey.iterateBorderMaxZ(GRID_ITERATOR{
            double u1=Ey(it,ix,iy,iz-1);
            double u2=Ey(it-1,ix,iy,iz-2);
            double u3=Ey(it-2,ix,iy,iz-3);
            double u4=Ey(it-3,ix,iy,iz-4);

            double d1=u1-u2;
            double d2=u2-u3;
            double d3=u3-u4;

            double dd1=d1-d2;
            double dd2=d2-d3;

            double ddd1=dd1-dd2;

            Ey(it+1,ix,iy,iz)=u1+d1+dd1+ddd1;
        });


        //cout << "Ex9" << endl;
        Ez.iterateBorderMinX(GRID_ITERATOR{
            double u1=Ez(it,ix+1,iy,iz);
            double u2=Ez(it-1,ix+2,iy,iz);
            double u3=Ez(it-2,ix+3,iy,iz);
            double u4=Ez(it-3,ix+4,iy,iz);

            double d1=u1-u2;
            double d2=u2-u3;
            double d3=u3-u4;

            double dd1=d1-d2;
            double dd2=d2-d3;

            double ddd1=dd1-dd2;

            Ez(it+1,ix,iy,iz)=u1+d1+dd1+ddd1;
        });
        //cout << "Ex10" << endl;
        Ez.iterateBorderMaxX(GRID_ITERATOR{
            double u1=Ez(it,ix-1,iy,iz);
            double u2=Ez(it-1,ix-2,iy,iz);
            double u3=Ez(it-2,ix-3,iy,iz);
            double u4=Ez(it-3,ix-4,iy,iz);

            double d1=u1-u2;
            double d2=u2-u3;
            double d3=u3-u4;

            double dd1=d1-d2;
            double dd2=d2-d3;

            double ddd1=dd1-dd2;

            Ez(it+1,ix,iy,iz)=u1+d1+dd1+ddd1;
        });


        //cout << "Ex11" << endl;
        Ez.iterateBorderMinY(GRID_ITERATOR{
            double u1=Ez(it,ix,iy+1,iz);
            double u2=Ez(it-1,ix,iy+2,iz);
            double u3=Ez(it-2,ix,iy+3,iz);
            double u4=Ez(it-3,ix,iy+4,iz);

            double d1=u1-u2;
            double d2=u2-u3;
            double d3=u3-u4;

            double dd1=d1-d2;
            double dd2=d2-d3;

            double ddd1=dd1-dd2;

            Ez(it+1,ix,iy,iz)=u1+d1+dd1+ddd1;
        });
        //cout << "Ex12" << endl;
        Ez.iterateBorderMaxY(GRID_ITERATOR{
            double u1=Ez(it,ix,iy-1,iz);
            double u2=Ez(it-1,ix,iy-2,iz);
            double u3=Ez(it-2,ix,iy-3,iz);
            double u4=Ez(it-3,ix,iy-4,iz);

            double d1=u1-u2;
            double d2=u2-u3;
            double d3=u3-u4;

            double dd1=d1-d2;
            double dd2=d2-d3;

            double ddd1=dd1-dd2;

            Ez(it+1,ix,iy,iz)=u1+d1+dd1+ddd1;
        });




        cout << "Hx" << endl;
        Hx.iterateWhole(GRID_ITERATOR{
            double MU = 0.5*(mu(0,ix-0.5,iy,iz)+mu(0,ix+0.5,iy,iz));

            Hx(it+1.5, ix, iy, iz)=
                Hx(it+0.5, ix, iy, iz) -
                DT*C/MU  * (
                (Ez(it+1, ix, iy+0.5, iz) - Ez(it+1, ix, iy-0.5, iz)) / dy
                -
                (Ey(it+1, ix, iy, iz+0.5) - Ey(it+1, ix, iy, iz-0.5)) / dz
                );
        });

        cout << "Hy" << endl;
        Hy.iterateWhole(GRID_ITERATOR{
            double MU = 0.5*(mu(0,ix,iy-0.5,iz)+mu(0,ix,iy+0.5,iz));

            Hy(it+1.5, ix, iy, iz)=
                Hy(it+0.5, ix, iy, iz) -
                DT*C/MU  * (
                (Ex(it+1, ix, iy, iz+0.5) - Ex(it+1, ix, iy, iz-0.5)) / dz
                -
                (Ez(it+1, ix+0.5, iy, iz) - Ez(it+1, ix-0.5, iy, iz)) / dx
                );
        });


        cout << "Hz" << endl;
        Hz.iterateWhole(GRID_ITERATOR{
            double MU = 0.5*(mu(0,ix,iy,iz-0.5)+mu(0,ix,iy,iz+0.5));

            Hz(it+1.5, ix, iy, iz)=
                Hz(it+0.5, ix, iy, iz) -
                DT*C/MU  * (
                (Ey(it+1, ix+0.5, iy, iz) - Ey(it+1, ix-0.5, iy, iz)) / dx
                -
                (Ex(it+1, ix, iy+0.5, iz) - Ex(it+1, ix, iy-0.5, iz)) / dy
                );
        });


















		if(((int)it % SAVE_EVERY_N) == 0){

            if(saveVTK){
                Hz[it+0.5].saveVTK(SAVE_DIR+"Hz"+toZPadString((int)it,5)+".vtk");
            }

            if(saveBMP){
                Hz[it+0.5].saveZSlice(Nz/2,SAVE_DIR+"Hz"+toZPadString((int)it,5)+".bmp");
            }
		}

	}
	mainTimer.logTime("calculation finished");


	pressAnyKey();
}

void mainHD();





int main(int argc, char *argv[]){

    Configuration cfg("conf.ini");

    THREAD_NUM=cfg.getInt("threadNum");
    SAVE_DIR=cfg.getString("saveDir");
    SAVE_EVERY_N=cfg.getInt("saveEvery");

    saveBMP=cfg.getBoolean("saveBMP");
    saveVTK=cfg.getBoolean("saveVTK");

    cout << "saveDir: " << SAVE_DIR << endl;
    cout << "saveEvery: " << SAVE_EVERY_N << endl;
    cout << "threadNum: " << THREAD_NUM << endl;

    omp_set_num_threads(THREAD_NUM);

	mainFDTD();
    //testConfiguration();
	//testClosure();
	//testIterateBorder();
    //testIterateBorder2();
    //testIterateBorder3();
    //mainHD();

    cout<<"calculation completed"<<endl;
    pressAnyKey();
	return 0;

}


const double GAMMA=5.0/3.0;

struct HydroDynVector3D{
    double density;
    double momentumX;
    double momentumY;
    double momentumZ;
    double fullEnergy;

    HydroDynVector3D(){
        density=0.0;
        momentumX=0.0;
        momentumY=0.0;
        momentumZ=0.0;
        fullEnergy=0.0;
    }

    /*
    HydroDynVector3D(double _density,double _e,double _Vx,double _Vy,double _Vz){
        density=_density;
        Vx=_Vx;
        Vy=_Vy;
        Vz=_Vz;
        e=_e;
    } */
};


 struct Vec5D{
    double c0;
    double c1;
    double c2;
    double c3;
    double c4;

    Vec5D(){
        c0=0.0;
        c1=0.0;
        c2=0.0;
        c3=0.0;
        c4=0.0;
    }

    Vec5D(double _c0,double _c1,double _c2,double _c3,double _c4){
        c0=_c0;
        c1=_c1;
        c2=_c2;
        c3=_c3;
        c4=_c4;
    }
 };

 Vec5D U_from_Rho_P_V(double rho,double P,double Vx,double Vy,double Vz){
        double eps=P/rho/(GAMMA-1.0);
        double e=rho*eps + rho*sqr(Vx,Vy,Vz)/2.0;
        return Vec5D(rho,rho*Vx,rho*Vy,rho*Vz,e);
 }



 Vec5D operator+(Vec5D a,Vec5D b){
     Vec5D res;
     res.c0 = a.c0 + b.c0;
     res.c1 = a.c1 + b.c1;
     res.c2 = a.c2 + b.c2;
     res.c3 = a.c3 + b.c3;
     res.c4 = a.c4 + b.c4;
     return res;
 }

 Vec5D operator-(Vec5D a,Vec5D b){
     Vec5D res;
     res.c0 = a.c0 - b.c0;
     res.c1 = a.c1 - b.c1;
     res.c2 = a.c2 - b.c2;
     res.c3 = a.c3 - b.c3;
     res.c4 = a.c4 - b.c4;
     return res;
 }

 Vec5D operator*(Vec5D v,double a){
     Vec5D res;
     res.c0 = v.c0*a;
     res.c1 = v.c1*a;
     res.c2 = v.c2*a;
     res.c3 = v.c3*a;
     res.c4 = v.c4*a;
     return res;
 }

 Vec5D operator*(double a,Vec5D v){
     return v*a;
 }

 Vec5D operator/(Vec5D v,double a){
     Vec5D res;
     res.c0 = v.c0/a;
     res.c1 = v.c1/a;
     res.c2 = v.c2/a;
     res.c3 = v.c3/a;
     res.c4 = v.c4/a;
     return res;
 }

 double rho_from_U(Vec5D v){
     return v.c0;
 }

 double Vx_from_U(Vec5D v){
     return v.c1/v.c0;
 }

 double Vy_from_U(Vec5D v){
     return v.c2/v.c0;
 }

 double Vz_from_U(Vec5D v){
     return v.c3/v.c0;
 }

 double e_from_U(Vec5D v){
     return v.c4;
 }

 double eps_from_U(Vec5D v){
     double e=e_from_U(v);
     double rho=rho_from_U(v);
     double sqrV=sqr(v.c1,v.c2,v.c3);
     return e/rho - sqrV/2.0;
 }

 double P_from_U(Vec5D v){
     double e=e_from_U(v);
     double rho=rho_from_U(v);
     double Vx=Vx_from_U(v);
     double Vy=Vy_from_U(v);
     double Vz=Vz_from_U(v);

     return (GAMMA-1.0)*(e-rho*sqr(Vx,Vy,Vz)/2.0);//(GAMMA-1.0)*rho*eps;
 }

 Vec5D U_to_Fx(const Vec5D& u){
     double rho=rho_from_U(u);
     double Vx=Vx_from_U(u);
     double Vy=Vy_from_U(u);
     double Vz=Vz_from_U(u);
     double e=e_from_U(u);
     double eps=eps_from_U(u);
     double P=(GAMMA-1.0)*(e-rho*sqr(Vx,Vy,Vz)/2.0);//(GAMMA-1.0)*rho*eps;

     Vec5D res(rho*Vx,
               rho*Vx*Vx+P,
               rho*Vx*Vy,
               rho*Vx*Vz,
               (e+P)*Vx);

     return res;
 }

  Vec5D U_to_Fy(Vec5D u){
     double rho=rho_from_U(u);
     double Vx=Vx_from_U(u);
     double Vy=Vy_from_U(u);
     double Vz=Vz_from_U(u);
     double e=e_from_U(u);
     double eps=eps_from_U(u);
     double P=(GAMMA-1.0)*rho*eps;

     Vec5D res(rho*Vy,
               rho*Vy*Vx,
               rho*Vy*Vy+P,
               rho*Vy*Vz,
               (e+P)*Vy);

     return res;
 }

  Vec5D U_to_Fz(Vec5D u){
     double rho=rho_from_U(u);
     double Vx=Vx_from_U(u);
     double Vy=Vy_from_U(u);
     double Vz=Vz_from_U(u);
     double e=e_from_U(u);
     double eps=eps_from_U(u);
     double P=(GAMMA-1.0)*rho*eps;

     Vec5D res(rho*Vz,
               rho*Vz*Vx,
               rho*Vy*Vz,
               rho*Vz*Vz+P,
               (e+P)*Vz);

     return res;
 }


 double soundVelocity(double gamma,double P,double rho){
     return sqrt(gamma*P/rho);
 }

