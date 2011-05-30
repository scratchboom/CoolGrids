#include "GridsCommon.hpp"

int main(){

	auto f=[&](double t)->double{
		return sin(t);
	};

	DBGVAL(solve(f,1.0,0.01));

	return 0;
}
