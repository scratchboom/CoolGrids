#include <iostream>
#include <vector>

using namespace std;


class Base{
public:

	Base(){
		cout << "Base constructor" << endl;
	}

	virtual void method(){
		cout << "Base method" << endl;
	}

	virtual ~Base(){
		cout << "Base destructor" << endl;
	}
};

class Derived : public Base{
public:

	Derived(){
		cout << "Derived constructor" << endl;
	}

	void method(){
		cout << "Derived method" << endl;
	}

	virtual ~Derived(){
			cout << "Derived destructor" << endl;
	}
};


int main(){

	Base base;
	Derived derived;

	Base* ref = &derived;

	ref->method();

	ref = &base;

	ref->method();


	return 0;
}
