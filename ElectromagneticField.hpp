
#pragma once


struct EH{
	double Ex,Ey,Ez;
	double Hx,Hy,Hz;
};


class ElectromagneticField{
public:
	virtual EH getField(double t,double DT,double x,double y,double z){return EH();};//TODO some HACK
};

class Superposition:public ElectromagneticField{
private:
	ElectromagneticField* f1;
	ElectromagneticField* f2;

public:

	Superposition(ElectromagneticField* f1,ElectromagneticField* f2){
		this->f1= f1;
		this->f2= f2;
	}

	EH getField(double t,double DT,double x,double y,double z){
		EH res;

		EH f1Val=f1->getField(t,DT,x,y,z);
		EH f2Val=f2->getField(t,DT,x,y,z);

		res.Ex=f1Val.Ex + f2Val.Ex;
		res.Ey=f1Val.Ey + f2Val.Ey;
		res.Ez=f1Val.Ez + f2Val.Ez;

		res.Hx=f1Val.Hx + f2Val.Hx;
		res.Hy=f1Val.Hy + f2Val.Hy;
		res.Hz=f1Val.Hz + f2Val.Hz;

		return res;
	}
};

ElectromagneticField* superposition(ElectromagneticField* f1,ElectromagneticField* f2){
	return new Superposition(f1,f2);
}

template <typename FX,typename FY,typename FZ>
class PointwiseChargeField:public ElectromagneticField{
public:
	double e;
	FX sourceX;
	FY sourceY;
	FZ sourceZ;

	PointwiseChargeField(double e,FX sourceX,FY sourceY,FZ sourceZ){
		this->e=e;
		this->sourceX=sourceX;
		this->sourceY=sourceY;
		this->sourceZ=sourceZ;
	}

	EH getField(double t,double DT,double x,double y,double z){

		double C=1;//TODO change this - VERY IMPORTANT!!!

		auto RRR=[&](double tt){
			return C*(tt-t) + length(x-sourceX(tt),y-sourceY(tt),z-sourceZ(tt));
		};

		double t_=solve(RRR,length(x-sourceX(t),y-sourceY(t),z-sourceZ(t)),DT);

		EH res;

		double X=x-sourceX(t_);
		double Y=y-sourceY(t_);
		double Z=z-sourceZ(t_);

		double R=length(X,Y,Z);

		double Vx= (sourceX(t_+DT)-sourceX(t_-DT))/2.0/DT;
		double Vy= (sourceY(t_+DT)-sourceY(t_-DT))/2.0/DT;
		double Vz= (sourceZ(t_+DT)-sourceZ(t_-DT))/2.0/DT;

		double ax= ((sourceX(t_+DT)-2.0*sourceX(t_)+sourceX(t_-DT))/sqr(DT)) * (1.0 - Vx*X/R/C);
		double ay= ((sourceY(t_+DT)-2.0*sourceY(t_)+sourceY(t_-DT))/sqr(DT)) * (1.0 - Vy*Y/R/C);
		double az= ((sourceZ(t_+DT)-2.0*sourceZ(t_)+sourceZ(t_-DT))/sqr(DT)) * (1.0 - Vz*Z/R/C);



		double RV=X*Vx+Y*Vy+Z*Vz;

		double sc=(X-Vx*R/C)*ax + (Y-Vy*R/C)*ay + (Z-Vz*R/C)*az;

		res.Ex=e*(1.0-sqr(Vx,Vy,Vz)/sqr(C))/cube(R-RV/C)*(X-Vx*R/C)
			+e/sqr(C)/cube(R-RV/C)*X*sc;

		res.Ey=e*(1.0-sqr(Vx,Vy,Vz)/sqr(C))/cube(R-RV/C)*(Y-Vy*R/C)
				+e/sqr(C)/cube(R-RV/C)*Y*sc;

		res.Ez=e*(1.0-sqr(Vx,Vy,Vz)/sqr(C))/cube(R-RV/C)*(Z-Vz*R/C)
				+e/sqr(C)/cube(R-RV/C)*Z*sc;

		/*
		i  j  k
		X  Y  Z
		Ex Ey Ez
		*/

		res.Hx=(Y*res.Ez-Z*res.Ey)/R;
		res.Hy=(Z*res.Ex-X*res.Ez)/R;
		res.Hz=(X*res.Ey-Y*res.Ex)/R;

		return res;
	}


};

template<typename FX,typename FY,typename FZ>
ElectromagneticField* createPointwiseChargeField(double e,FX x,FY y,FZ z){
	return new PointwiseChargeField<FX,FY,FZ>(e,x,y,z);
}
