#pragma once

long getTickCount(){
	struct timeval tv;
    gettimeofday(&tv,NULL);
    return (tv.tv_sec*1000+tv.tv_usec/1000);
}

class Timer{

private:

    unsigned long startTick;

public:

    Timer(){
    	reset();
    }

    double getTime(){
    	return (getTickCount()-startTick)/1000.0;
    }
    void reset(){
        startTick=getTickCount();
    };

	void logTime(const std::string& message){
		std::cout<<getTime() <<": "<<message << std::endl;
	}
};



