#pragma once

template<typename T>
class readonly{
private:
	T value;

public:
	void setValue(const T& value){
		this->value=value;
	}

	operator T const &() const{
		return value;
	}
};

