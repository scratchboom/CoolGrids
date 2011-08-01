#include <iostream>
#include <vector>
#include <boost/filesystem.hpp>

using namespace std;
using namespace boost::filesystem;

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

	path p = "/home";

	if (exists(p)) // does p actually exist?
			{
		if (is_regular_file(p)) // is p a regular file?
			cout << p << " size is " << file_size(p) << '\n';
		else if (is_directory(p)) // is p a directory?
			cout << p << " is a directory\n";
		else
			cout
					<< p
					<< " exists, but is neither a regular file nor a directory\n";
	} else
		cout << p << " does not exist\n";



	return 0;
}
