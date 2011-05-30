#include <iostream>



using namespace std;

void printFunctionTable(int (*func)(int)){
	for(int i=1;i<=10;i++) cout << func(i) << "  ";
	cout << endl;
}

int square(int x){
	return x*x;
}

template <typename CL>
int closureToFunction(int x){
	CL cl;
	return cl(x);
}

#define CLOSURE_TO_FUNCTION(cl) closureToFunction<decltype(cl)>

int main(void)
{
	printFunctionTable(square);

	auto cube=[](int x){return x*x*x;};
	auto quad=[](int x){return x*x*x*x;};

	int a=2;
	auto closure=[=](int x){return a*x;};


	struct structure{
		static int calc(int x){return x*x;};
	} rhs;


	struct St{
		double a;

		St(int a){
			this->a = a;
		};

		int calc(int x){
			return a*x*x;
		};

		static int staticCalc(St *zis,int x){
			return this->calc(x)
		};

	} st(5);

	printFunctionTable(CLOSURE_TO_FUNCTION(cube));
	printFunctionTable(CLOSURE_TO_FUNCTION(quad));
	printFunctionTable(CLOSURE_TO_FUNCTION(closure));

	printFunctionTable(structure::calc);
	printFunctionTable(st.calc);


	return 0;
}
