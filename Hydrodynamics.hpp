#pragma once

template<typename T,int N>
struct Vector{
public:
	T array[N];

public:
	Vector(){

	}

	Vector(T arr[]){
		for(int i=0;i<N;i++)array[i]=arr[i];
	}

	Vector(T v0,T v1,T v2){
		array[0]=v0;
		array[1]=v1;
		array[2]=v2;
	}

	T& operator[](int i){
		return array[i];
	}

	void clear(){
		for(int i=0;i<N;i++)array[i]=T();
	}
};

template<typename T,int N>
std::string toString(const Vector<T,N>& vec){
	std::string res("{");
	for(int i=0;i<N-1;i++)res.append(toString(vec.array[i])).append(",");
	res.append(toString(vec.array[N-1]));
	res.append("}");
	return res;
}


template<typename T,int N>
Vector<T,N> operator+(const Vector<T,N>& a,const Vector<T,N>& b){
	Vector<T,N> res;
	for(int i=0;i<N;i++)res[i]=a.array[i]+b.array[i];
	return res;
}

template<typename T,int N>
Vector<T,N> operator-(const Vector<T,N>& a,const Vector<T,N>& b){
	Vector<T,N> res;
	for(int i=0;i<N;i++)res[i]=a.array[i]-b.array[i];
	return res;
}

template<typename T,int N>
Vector<T,N> operator-(const Vector<T,N>& a){
	Vector<T,N> res;
	res.clear();
	for(int i=0;i<N;i++)res[i]=-a.array[i];
	return res;
}

template<typename T,int N>
Vector<T,N> operator*(double a,const Vector<T,N>& v){
	Vector<T,N> res;
	for(int i=0;i<N;i++)res[i]=a*v.array[i];
	return res;
}

template<typename T,int N>
Vector<T,N> operator*(const Vector<T,N>& v,double a){
	Vector<T,N> res;
	for(int i=0;i<N;i++)res[i]=a*v.array[i];
	return res;
}

template<typename T,int N>
Vector<T,N> operator/(const Vector<T,N>& v,double a){
	Vector<T,N> res;
	for(int i=0;i<N;i++)res[i]=v.array[i]/a;
	return res;
}

typedef Vector<double,1> Vec1D;
typedef Vector<double,2> Vec2D;
typedef Vector<double,3> Vec3D;

const double GAMMA=5.0/3.0;

struct HydroDynVec{
	static Vec3D from_Density_Pressure_Velocity(double rho,double p,double v){
		Vec3D res(1,1,1);
		res[0]=rho;
		res[1]=rho*v;
		res[2]=rho*(1.0/(GAMMA-1.0)*p/rho);
		return res;
	}

	//rho
	static double density(const Vec3D& U){
		return U.array[0];
	}

	//p
	static double pressure(const Vec3D& U){
		double rho=density(U);
		double v=velocityX(U);
		double E=fullEnergyPerVolumeUnit(U);
		return (E-rho*sqr(v)/2.0)*(GAMMA-1.0);
	}

	//Vx
	static double velocityX(const Vec3D& U){
		return U.array[1]/density(U);
	}

	//E
	static double fullEnergyPerVolumeUnit(const Vec3D& U){
		return U.array[2];
	}

	//rho*V^2/2
	static double kineticEnergyPerVolumeUnit(const Vec3D& U){
		return density(U)*sqr(velocityX(U))/2.0;
	}

	//e
	static double internalEnergyPerVolumeUnit(const Vec3D& U){
		return fullEnergyPerVolumeUnit(U)-kineticEnergyPerVolumeUnit(U);
	}

	//eps
	static double internalEnergyPerMassUnit(const Vec3D& U){
		return internalEnergyPerVolumeUnit(U)/density(U);
	}


	static Vec3D toFlow(const Vec3D& U){
		double rho=density(U);
		double v=velocityX(U);
		double E=fullEnergyPerVolumeUnit(U);
		double p=pressure(U);

		Vec3D res;
		res[0]=rho*v;
		res[1]=rho*sqr(v)+p;
		res[2]=(E+p)*v;
		return res;
	}
};

struct FlowVec{
	static Vec3D from_Density_Pressure_Velocity(double rho,double p,double v){
		Vec3D res;
		res[0]=rho*v;
		res[1]=rho*sqr(v)+p;
		res[2]=(rho*(1.0/(GAMMA-1.0)*p/rho +v*v/2.0)+p)*v;
		return res;
	}
};


double soundVelocity(double gamma,double P,double rho){
    return sqrt(gamma*P/rho);
}



double interpLiao3(double u1,double u2,double u3,double u4){
	double d1=u1-u2;
	double d2=u2-u3;
	double d3=u3-u4;

	double dd1=d1-d2;
	double dd2=d2-d3;

	double ddd1=dd1-dd2;

	return u1+d1+dd1+ddd1;
}

