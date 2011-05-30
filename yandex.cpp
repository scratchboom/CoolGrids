#include <iostream>
#include <vector>

using namespace std;


class Base{
public:

	Base(){
			cout << "Base constructor" << endl;
	}

	void method(){
		cout << "Base" << endl;
	}

	virtual ~Base(){
		cout << "Base destructor" << endl;
	}
};

class Derived : Base{
public:

	Derived(){
		cout << "Derived constructor" << endl;
	}

	void method(){
		cout << "Derived" << endl;
	}

	virtual ~Derived(){
			cout << "Derived destructor" << endl;
	}
};


int main(){

//	Base base;
//	Derived derived;
//
//	base.method();
//	derived.method();
//
//	vector<Base> values;
	//values.push_back(base);
	//values.push_back(derived);


	//Base *b=new Base();
	//delete b;

//	Derived *d=new Derived();
//	delete d;

	Derived *dArr=new Derived[5];
	delete[] dArr;


	return 0;
}
