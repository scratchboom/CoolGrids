#pragma once

template<typename T>
class Accumulator{

private:
	T minValue;
	T maxValue;

public:

	void add(const T& value){
		maxValue=max(maxValue,value);
		minValue=min(minValue,value);
	}

	T getMinValue(){
		return minValue;
	}

    T getMaxValue(){
		return maxValue;
	}

	T getMaxAbs(){
		return maxAbs(minValue,maxValue);
	}

};
