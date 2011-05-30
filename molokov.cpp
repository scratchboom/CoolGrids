#include <iostream>

using namespace std;


class TestClass {

private:
	double tmp;

public:
	TestClass() {
		//tmp = new double[3];//зачем тут память динамически выделять? 3 элемента можно и статически выделить
	}

public:
	void remove() {
		//delete[] tmp;//если выделять статически - эта строчка будет не нужна, если всё же хочешь динамически, перенеси в деструктор, или вызови в деструкторе этот метод, если он тебе нужен отдельно
		//delete this;//зачем пилить сук на котором сидишь?
	}

	~TestClass(){
		remove();
	}

};



int main(void) {
	int N = 5000000;

	TestClass** array = new TestClass*[N];

	for (int i = 0; i < N; i++) {
		array[i] = new TestClass();
	}


	sleep(10);


	//for (int i = 0; i < N; i++)array[i]->remove();
	for (int i = 0; i < N; i++) delete array[i];


	delete[] array;

	sleep(100);

	cout << "Done!" << endl;

}
;
