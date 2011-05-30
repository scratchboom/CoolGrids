#pragma once

bool MEGADEBUG = false;

#define DBGLN {std::cout << "DBGLN: " << __FILE__ <<"(" << __LINE__ <<")"<< std::endl;};
#define DBGVAL(val) {if(MEGADEBUG){std::cout << #val << ": " << (val) <<std::endl;}};
#define DBGLNVAL(val) {{if(MEGADEBUG){std::cout << "DBGLN: " << __FILE__ <<": " << __LINE__ << "    " << #val << ": " << (val) << std::endl;}};

#define CHECK(what,op,val) {std::cout << "checking " << #what << #op << #val; std::cout << ":    " #what << " = " << what; if(what op val) std::cout << " - OK (" << #what << #op << val << ")" << std::endl;else {std::cout << " - FAILED (should be " << #what << #op << val << ")" << std::endl; exit(0);}};

template<typename T>
std::string toString(T v){
	std::stringstream s;
	s << v;
	return s.str();
}

template<typename T>
std::string toZPadString(T v,int w){
	std::stringstream s;
	s << std::setw(w) << std::setfill('0') << v;
	return s.str();
}

template<typename T>
bool isInteger(T x){
	return x==floor(x);
}

template<typename T>
bool inRangeInclusive(T x,T xmin,T xmax){
	return (x>=xmin)&&(x<=xmax);
}

double rnd(){
    return rand()/(double(RAND_MAX)+1);
}

