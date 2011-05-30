#pragma once

class ByteBuffer{

private:

	size_t capacity;
	size_t size;
	char* buffer;

public:

	ByteBuffer(){
		buffer=NULL;
		clear();
	}

	void clear(){
		capacity=1024;
        size=0;
        delete buffer;
        buffer=new char[capacity];
	}

	template<typename T>
	void write(T value){
		while(size+sizeof(value)>capacity){
			capacity*=2;
			if(size+sizeof(value)<=capacity){
				char* newbuffer=new char[capacity];
				delete buffer;
				buffer=newbuffer;
			}
		}
		memcpy(buffer+size,&value,sizeof(value));
		size+=sizeof(value);
	}

	char* getBuffer(){
		return buffer;
	}

	size_t getSize(){
		return size;
	}


};
