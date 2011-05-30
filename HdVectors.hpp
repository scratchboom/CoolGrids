#pragma once

typedef boost::numeric::ublas::vector<double, boost::numeric::ublas::bounded_array<double, 2> > Vector2D;
typedef boost::numeric::ublas::vector<double, boost::numeric::ublas::bounded_array<double, 3> > Vector3D;
typedef boost::numeric::ublas::vector<double, boost::numeric::ublas::bounded_array<double, 4> > Vector4D;
typedef boost::numeric::ublas::vector<double, boost::numeric::ublas::bounded_array<double, 5> > Vector5D;

namespace HdVec1D {

Vector3D fromDensityPressureVelocity(double rho, double p, double v) {
	Vector3D res(3);
	res[0] = rho;
	res[1] = rho * v;
	res[2] = rho * (1.0 / (GAMMA - 1.0) * p / rho) + rho * sqr(v) / 2.0;
	return res;
}

//rho
static double density(const Vector3D& U) {
	return U[0];
}

//Vx
double velocityX(const Vector3D& U) {
	return U[1] / density(U);
}

//E
double fullEnergyPerVolumeUnit(const Vector3D& U) {
	return U[2];
}

//p
double pressure(const Vector3D& U) {
	double rho = density(U);
	double v = velocityX(U);
	double E = fullEnergyPerVolumeUnit(U);
	return (E - rho * sqr(v) / 2.0) * (GAMMA - 1.0);
}

//rho*V^2/2
double kineticEnergyPerVolumeUnit(const Vector3D& U) {
	return density(U) * sqr(velocityX(U)) / 2.0;
}

//e
double internalEnergyPerVolumeUnit(const Vector3D& U) {
	return fullEnergyPerVolumeUnit(U) - kineticEnergyPerVolumeUnit(U);
}

//eps
double internalEnergyPerMassUnit(const Vector3D& U) {
	return internalEnergyPerVolumeUnit(U) / density(U);
}

Vector3D toFlow(const Vector3D& U) {
	double rho = density(U);
	double v = velocityX(U);
	double E = fullEnergyPerVolumeUnit(U);
	double p = pressure(U);

	Vector3D res(3);
	res[0] = rho * v;
	res[1] = rho * sqr(v) + p;
	res[2] = (E + p) * v;
	return res;
}

}

namespace HdFlowVec1D {

namespace X{

Vector3D fromDensityPressureVelocity(double rho, double p, double v) {
	Vector3D res(3);
	res[0] = rho * v;
	res[1] = rho * sqr(v) + p;
	res[2] = (rho * (1.0 / (GAMMA - 1.0) * p / rho + v * v / 2.0) + p) * v;
	return res;
}


//returns X Flux vector between two cells
Vector3D riemannFlux(const Vector3D& u1,const Vector3D& u2){

	double pressure,density,energy,velocity;

	double rho1=HdVec1D::density(u1);
	double p1=HdVec1D::pressure(u1);
	double s1=soundVelocity(GAMMA,p1,rho1);

	double rho2=HdVec1D::density(u2);
	double p2=HdVec1D::pressure(u2);
	double s2=soundVelocity(GAMMA,p2,rho2);

	Real4 result=SolveRiemannProblem(HdVec1D::density(u1),
			HdVec1D::internalEnergyPerMassUnit(u1),
			HdVec1D::velocityX(u1),
					GAMMA,s1,
					HdVec1D::density(u2),
					HdVec1D::internalEnergyPerMassUnit(u2),
					HdVec1D::velocityX(u2),
					GAMMA,s2);
    pressure=Pressure(result);
    density=Density(result);
    velocity=Velocity(result);
    energy=Energy(result);

	return HdFlowVec1D::X::fromDensityPressureVelocity(density,pressure,velocity);
}

}

}

namespace HdVec2D {

Vector4D fromDensityPressureVelocity(double rho, double p, double vx, double vy) {
	Vector4D res(4);
	res[0] = rho;
	res[1] = rho * vx;
	res[2] = rho * vy;
	res[3] = rho * (1.0 / (GAMMA - 1.0) * p / rho) + rho * sqr(vx, vy) / 2.0;
	return res;
}

//rho
static double density(const Vector4D& U) {
	return U[0];
}

//Vx
double velocityX(const Vector4D& U) {
	return U[1] / density(U);
}

//Vy
double velocityY(const Vector4D& U) {
	return U[2] / density(U);
}

//E
double fullEnergyPerVolumeUnit(const Vector4D& U) {
	return U[3];
}

//p
double pressure(const Vector4D& U) {
	double rho = density(U);
	double vx = velocityX(U);
	double vy = velocityY(U);
	double E = fullEnergyPerVolumeUnit(U);
	return (E - rho * sqr(vx, vy) / 2.0) * (GAMMA - 1.0);
}

//rho*V^2/2
double kineticEnergyPerVolumeUnit(const Vector4D& U) {
	return density(U) * sqr(velocityX(U), velocityY(U)) / 2.0;
}

//e
double internalEnergyPerVolumeUnit(const Vector4D& U) {
	return fullEnergyPerVolumeUnit(U) - kineticEnergyPerVolumeUnit(U);
}

//eps
double internalEnergyPerMassUnit(const Vector4D& U) {
	return internalEnergyPerVolumeUnit(U) / density(U);
}

Vector4D toFlowX(const Vector4D& U) {
	double rho = density(U);
	double vx = velocityX(U);
	double vy = velocityY(U);
	double E = fullEnergyPerVolumeUnit(U);
	double p = pressure(U);

	Vector4D res(4);
	res[0] = rho * vx;
	res[1] = rho * sqr(vx) + p;
	res[2] = rho * vx * vy;
	res[3] = (E + p) * vx;
	return res;
}

Vector4D toFlowY(const Vector4D& U) {
	double rho = density(U);
	double vx = velocityX(U);
	double vy = velocityY(U);
	double E = fullEnergyPerVolumeUnit(U);
	double p = pressure(U);

	Vector4D res(4);
	res[0] = rho * vy;
	res[1] = rho * vy * vx;
	res[2] = rho * sqr(vy) + p;
	res[3] = (E + p) * vy;
	return res;
}

}

namespace HdFlowVec2D {

namespace X {

Vector4D fromDensityPressureVelocity(double rho, double p, double vx, double vy) {
	Vector4D res(4);
	res[0] = rho * vx;
	res[1] = rho * sqr(vx) + p;
	res[2] = rho * vx * vy;
	res[3] = (rho * (1.0 / (GAMMA - 1.0) * p / rho + sqr(vx, vy) / 2.0) + p) * vx;
	return res;
}

//returns X Flux vector between two cells
Vector4D riemannFlux(const Vector4D& u1,const Vector4D& u2){

	double pressure,density,energy,velocity;

	double rho1=HdVec2D::density(u1);
	double p1=HdVec2D::pressure(u1);
	double s1=soundVelocity(GAMMA,p1,rho1);

	double rho2=HdVec2D::density(u2);
	double p2=HdVec2D::pressure(u2);
	double s2=soundVelocity(GAMMA,p2,rho2);

	Real4 result=SolveRiemannProblem(HdVec2D::density(u1),
			HdVec2D::internalEnergyPerMassUnit(u1),
			HdVec2D::velocityX(u1),
					GAMMA,s1,
					HdVec2D::density(u2),
					HdVec2D::internalEnergyPerMassUnit(u2),
					HdVec2D::velocityX(u2),
					GAMMA,s2);
    pressure=Pressure(result);
    density=Density(result);
    velocity=Velocity(result);
    energy=Energy(result);

	return HdFlowVec2D::X::fromDensityPressureVelocity(density,pressure,velocity,0);
}

}

namespace Y {

Vector4D fromDensityPressureVelocity(double rho, double p, double vx, double vy) {
	Vector4D res(4);
	res[0] = rho * vy;
	res[1] = rho * vy * vx;
	res[2] = rho * sqr(vy) + p;
	res[3] = (rho * (1.0 / (GAMMA - 1.0) * p / rho + sqr(vx, vy) / 2.0) + p) * vy;
	return res;
}

//returns Y Flux vector between two cells
Vector4D riemannFlux(const Vector4D& u1,const Vector4D& u2){

	double pressure,density,energy,velocity;

	double rho1=HdVec2D::density(u1);
	double p1=HdVec2D::pressure(u1);
	double s1=soundVelocity(GAMMA,p1,rho1);

	double rho2=HdVec2D::density(u2);
	double p2=HdVec2D::pressure(u2);
	double s2=soundVelocity(GAMMA,p2,rho2);

	Real4 result=SolveRiemannProblem(HdVec2D::density(u1),
			HdVec2D::internalEnergyPerMassUnit(u1),
			HdVec2D::velocityY(u1),
					GAMMA,s1,
					HdVec2D::density(u2),
					HdVec2D::internalEnergyPerMassUnit(u2),
					HdVec2D::velocityY(u2),
					GAMMA,s2);
    pressure=Pressure(result);
    density=Density(result);
    velocity=Velocity(result);
    energy=Energy(result);

	return HdFlowVec2D::Y::fromDensityPressureVelocity(density,pressure,0,velocity);
}



}

}

namespace HdVec3D {

Vector5D fromDensityPressureVelocity(double rho, double p, double vx, double vy, double vz) {
	Vector5D res(5);
	res[0] = rho;
	res[1] = rho * vx;
	res[2] = rho * vy;
	res[3] = rho * vz;
	res[4] = rho * (1.0 / (GAMMA - 1.0) * p / rho) + rho * sqr(vx, vy, vz) / 2.0;
	return res;
}

Vector5D fromDesityVelocityInternalEnergyPerVolumeUnit(double rho, double vx, double vy, double vz,double internalEnergyPerVolumeUnit){
	Vector5D res(5);
	res[0] = rho;
	res[1] = rho * vx;
	res[2] = rho * vy;
	res[3] = rho * vz;
	res[4] = rho*sqr(vx,vy,vz)/2.0 + internalEnergyPerVolumeUnit;
	return res;
}

//rho
static double density(const Vector5D& U) {
	return U[0];
}

//Vx
double velocityX(const Vector5D& U) {
	return U[1] / density(U);
}

//Vy
double velocityY(const Vector5D& U) {
	return U[2] / density(U);
}

//Vz
double velocityZ(const Vector5D& U) {
	return U[3] / density(U);
}

//E
double fullEnergyPerVolumeUnit(const Vector5D& U) {
	return U[4];
}

//p
double pressure(const Vector5D& U) {
	double rho = density(U);
	double vx = velocityX(U);
	double vy = velocityY(U);
	double vz = velocityZ(U);
	double E = fullEnergyPerVolumeUnit(U);
	return (E - rho * sqr(vx, vy, vz) / 2.0) * (GAMMA - 1.0);
}

//rho*V^2/2
double kineticEnergyPerVolumeUnit(const Vector5D& U) {
	return density(U) * sqr(velocityX(U), velocityY(U), velocityZ(U)) / 2.0;
}

//e
double internalEnergyPerVolumeUnit(const Vector5D& U) {
	return fullEnergyPerVolumeUnit(U) - kineticEnergyPerVolumeUnit(U);
}

//eps
double internalEnergyPerMassUnit(const Vector5D& U) {
	return internalEnergyPerVolumeUnit(U) / density(U);
}

Vector5D toFlowX(const Vector5D& U) {
	double rho = density(U);
	double vx = velocityX(U);
	double vy = velocityY(U);
	double vz = velocityZ(U);
	double E = fullEnergyPerVolumeUnit(U);
	double p = pressure(U);

	Vector5D res(5);
	res[0] = rho * vx;
	res[1] = rho * sqr(vx) + p;
	res[2] = rho * vx * vy;
	res[3] = rho * vx * vz;
	res[4] = (E + p) * vx;
	return res;
}

Vector5D toFlowY(const Vector5D& U) {
	double rho = density(U);
	double vx = velocityX(U);
	double vy = velocityY(U);
	double vz = velocityZ(U);
	double E = fullEnergyPerVolumeUnit(U);
	double p = pressure(U);

	Vector5D res(5);
	res[0] = rho * vy;
	res[1] = rho * vy * vx;
	res[2] = rho * sqr(vy) + p;
	res[3] = rho * vy * vz;
	res[4] = (E + p) * vy;
	return res;
}

Vector5D toFlowZ(const Vector5D& U) {
	double rho = density(U);
	double vx = velocityX(U);
	double vy = velocityY(U);
	double vz = velocityZ(U);
	double E = fullEnergyPerVolumeUnit(U);
	double p = pressure(U);

	Vector5D res(5);
	res[0] = rho * vz;
	res[1] = rho * vz * vx;
	res[2] = rho * vz * vy;
	res[3] = rho * sqr(vz) + p;
	res[4] = (E + p) * vz;
	return res;
}

}

namespace HdFlowVec3D {

namespace X {

Vector5D fromDensityPressureVelocity(double rho, double p, double vx, double vy, double vz) {
	Vector5D res(5);
	res[0] = rho * vx;
	res[1] = rho * sqr(vx) + p;
	res[2] = rho * vx * vy;
	res[3] = rho * vx * vz;
	res[4] = (rho * (1.0 / (GAMMA - 1.0) * p / rho + sqr(vx, vy, vz) / 2.0) + p) * vx;
	return res;
}

Vector5D riemannFlux(const Vector5D& u1,const Vector5D& u2){

	double pressure,density,energy,velocity;

	double rho1=HdVec3D::density(u1);
	double p1=HdVec3D::pressure(u1);
	double s1=soundVelocity(GAMMA,p1,rho1);

	double rho2=HdVec3D::density(u2);
	double p2=HdVec3D::pressure(u2);
	double s2=soundVelocity(GAMMA,p2,rho2);


	Real4 result=SolveRiemannProblem(HdVec3D::density(u1),
			HdVec3D::internalEnergyPerMassUnit(u1),
			HdVec3D::velocityX(u1),
					GAMMA,s1,
					HdVec3D::density(u2),
					HdVec3D::internalEnergyPerMassUnit(u2),
					HdVec3D::velocityX(u2),
					GAMMA,s2);

    pressure=Pressure(result);
    density=Density(result);
    velocity=Velocity(result);
    energy=Energy(result);


	return HdFlowVec3D::X::fromDensityPressureVelocity(density,pressure,velocity,0,0);
}

}

namespace Y {

Vector5D fromDensityPressureVelocity(double rho, double p, double vx, double vy, double vz) {
	Vector5D res(5);
	res[0] = rho * vy;
	res[1] = rho * vy * vx;
	res[2] = rho * sqr(vy) + p;
	res[3] = rho * vy * vz;
	res[4] = (rho * (1.0 / (GAMMA - 1.0) * p / rho + sqr(vx, vy, vz) / 2.0) + p) * vy;
	return res;
}

Vector5D riemannFlux(const Vector5D& u1,const Vector5D& u2){

	double pressure,density,energy,velocity;

	double rho1=HdVec3D::density(u1);
	double p1=HdVec3D::pressure(u1);
	double s1=soundVelocity(GAMMA,p1,rho1);

	double rho2=HdVec3D::density(u2);
	double p2=HdVec3D::pressure(u2);
	double s2=soundVelocity(GAMMA,p2,rho2);

	Real4 result=SolveRiemannProblem(HdVec3D::density(u1),
			HdVec3D::internalEnergyPerMassUnit(u1),
			HdVec3D::velocityX(u1),
					GAMMA,s1,
					HdVec3D::density(u2),
					HdVec3D::internalEnergyPerMassUnit(u2),
					HdVec3D::velocityX(u2),
					GAMMA,s2);
    pressure=Pressure(result);
    density=Density(result);
    velocity=Velocity(result);
    energy=Energy(result);

	return HdFlowVec3D::Y::fromDensityPressureVelocity(density,pressure,0,velocity,0);
}

}

namespace Z {

Vector5D fromDensityPressureVelocity(double rho, double p, double vx, double vy, double vz) {
	Vector5D res(5);
	res[0] = rho * vz;
	res[1] = rho * vz * vx;
	res[2] = rho * vz * vy;
	res[3] = rho * sqr(vz) + p;
	res[4] = (rho * (1.0 / (GAMMA - 1.0) * p / rho + sqr(vx, vy, vz) / 2.0) + p) * vz;
	return res;
}

Vector5D riemannFlux(const Vector5D& u1,const Vector5D& u2){

	double pressure,density,energy,velocity;

	double rho1=HdVec3D::density(u1);
	double p1=HdVec3D::pressure(u1);
	double s1=soundVelocity(GAMMA,p1,rho1);

	double rho2=HdVec3D::density(u2);
	double p2=HdVec3D::pressure(u2);
	double s2=soundVelocity(GAMMA,p2,rho2);

	Real4 result=SolveRiemannProblem(HdVec3D::density(u1),
			HdVec3D::internalEnergyPerMassUnit(u1),
			HdVec3D::velocityX(u1),
					GAMMA,s1,
					HdVec3D::density(u2),
					HdVec3D::internalEnergyPerMassUnit(u2),
					HdVec3D::velocityX(u2),
					GAMMA,s2);
    pressure=Pressure(result);
    density=Density(result);
    velocity=Velocity(result);
    energy=Energy(result);

	return HdFlowVec3D::Z::fromDensityPressureVelocity(density,pressure,0,0,velocity);
}

}

}

