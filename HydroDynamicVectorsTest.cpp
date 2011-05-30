#include "GridsCommon.hpp"

using namespace std;

int main() {

	Vector3D vec = HdVec1D::fromDensityPressureVelocity(1, 2, 3);
	Vector3D flow = HdFlowVec1D::fromDensityPressureVelocity(1, 2, 3);

	Vector4D hdVec2Dx = HdVec2D::fromDensityPressureVelocity(1, 2, 3, 0);
	Vector4D hdVec2Dy = HdVec2D::fromDensityPressureVelocity(1, 2, 0, 3);

	Vector4D flowVec2Dx = HdFlowVec2D::X::fromDensityPressureVelocity(1, 2, 3, 0);
	Vector4D flowVec2Dy = HdFlowVec2D::Y::fromDensityPressureVelocity(1, 2, 0, 3);

	Vector5D hdVec3Dx = HdVec3D::fromDensityPressureVelocity(1, 2, 3, 0, 0);
	Vector5D hdVec3Dy = HdVec3D::fromDensityPressureVelocity(1, 2, 0, 3, 0);
	Vector5D hdVec3Dz = HdVec3D::fromDensityPressureVelocity(1, 2, 0, 0, 3);

	Vector5D flowVec3Dx = HdFlowVec3D::X::fromDensityPressureVelocity(1, 2, 3, 0, 0);
	Vector5D flowVec3Dy = HdFlowVec3D::Y::fromDensityPressureVelocity(1, 2, 0, 3, 0);
	Vector5D flowVec3Dz = HdFlowVec3D::Z::fromDensityPressureVelocity(1, 2, 0, 0, 3);

	cout << vec << endl;
	cout << HdVec1D::toFlow(vec) << endl;
	cout << flow << endl;

	cout << hdVec2Dx << endl;
	cout << hdVec2Dy << endl;

	cout << HdVec2D::toFlowX(hdVec2Dx) << endl;
	cout << HdVec2D::toFlowY(hdVec2Dy) << endl;

	cout << flowVec2Dx << endl;
	cout << flowVec2Dy << endl;

	cout << hdVec3Dx << endl;
	cout << hdVec3Dy << endl;
	cout << hdVec3Dz << endl;

	cout << HdVec3D::toFlowX(hdVec3Dx) << endl;
	cout << HdVec3D::toFlowY(hdVec3Dy) << endl;
	cout << HdVec3D::toFlowZ(hdVec3Dz) << endl;

	cout << flowVec3Dx << endl;
	cout << flowVec3Dy << endl;
	cout << flowVec3Dz << endl;

	return 0;
}
