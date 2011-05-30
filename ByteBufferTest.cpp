
#include "GridsCommon.hpp"

int main(){

	ByteBuffer buf;

	buf.write(1);
	std::cout << buf.getSize() << std::endl;
	buf.write("a");
	std::cout << buf.getSize() << std::endl;
	buf.write(1.123);
	std::cout << buf.getSize() << std::endl;
	buf.write(1.123f);
    std::cout << buf.getSize() << std::endl;
    buf.clear();
    std::cout << buf.getSize() << std::endl;


	return 0;
}
