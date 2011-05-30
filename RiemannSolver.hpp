#pragma once

typedef double Real;

#define Real(a) ((Real)a)
const unsigned int GNumLevelCount = 2;
struct Real4 {
	Real x, y, z, w;
};

struct Real2 {
	Real x, y;
};

inline Real4 make_Real4(Real x, Real y, Real z, Real w) {
	Real4 Result = { x, y, z, w };
	return Result;
};

#define Sqr(a) ((a)*(a))
#define Sign(a) (a > 0 ? 1 : (a < 0 ? -1 : 0))
#define Max(a, b) (a > b ? a : b)
#define Min(a, b) (a < b ? a : b)
#define Len2(a) (a[0] * a[0] + a[1] * a[1] + a[2] * a[2])
#define Len2Splited(a, b, c) (Sqr(a) + Sqr(b) + Sqr(c))
#define Len(a) (sqrt(Len2(a)))
#define LenSplited(a, b, c) (sqrt(Len2Splited(a, b, c)))

#define Epsilon                   Real(1e-20)    //точность вычислений
#define Bolcman                   Real(1.38e-23)//постоянная больцмана
#define UniversalGasConstant      Real(8.31)    //Универсальная газовая постоянная
#define DefGamma                  Real(1.4)     //показатель адиабаты
#define DefMolarMass              Real(0.029)   //молярная масса воздуха
#define DefImpirityHeatCapacity   Real(4200)    //Теплоемкость примеси
#define AirViscosityRatio         Real(2.43e-7) //коэффициент в вязкости воздуха M=AVR*T*0.75

#ifndef M_PI
#define M_PI Real(3.142715)
#endif

#define Pressure(AReal4) (AReal4).x
#define Velocity(AReal4) (AReal4).y
#define Density(AReal4) (AReal4).z
#define Energy(AReal4) (AReal4).w

#define CUDA "CUDA";//todo remove this

//New Riemann Solver
Real AdiabatVelocity(Real Pressure0, Real Velocity0, Real Density0,
		Real Sound0, Real Gamma, Real Pressure, bool Left) {
	bool Puasson = Pressure <= Pressure0;
	Real Sign_ = Left ? Real(-1) : Real(1);
	if (Puasson)
		return Velocity0 + Sign_ * 2 * Sound0 / (Gamma - 1) * (pow(Pressure	/ Pressure0, (Gamma - 1) / (2 * Gamma)) - 1);
	else
		return Velocity0 + Sign_ * (Pressure - Pressure0) * sqrt(2 / (Density0 * (Pressure * (Gamma + 1) + Pressure0 * (Gamma - 1))));
};

Real AdiabatVelocityDerivative(Real Pressure0, Real Velocity0, Real Density0,
		Real Sound0, Real Gamma, Real Pressure, bool Left) {
	bool Puasson = Pressure <= Pressure0;
	Real Sign_ = Left ? Real(-1) : Real(1);
	if (Puasson)
		return Sign_ * Sound0 / (Gamma * Pressure0) * pow(Pressure / Pressure0,	-(Gamma + 1) / (2 * Gamma));
	else {
		Real PowExp = Pressure * (Gamma + 1) + Pressure0 * (Gamma - 1);
		return Sign_ * (Pressure * (Gamma + 1) + Pressure0 * (3 * Gamma - 1))/ sqrt(2 * Density0 * PowExp * PowExp * PowExp);
	};
};

//return Pressure and Velocity
Real2 NewtonVelocity(Real PressureLeft, Real VelocityLeft, Real DensityLeft,
		Real GammaLeft, Real SoundLeft, Real PressureRight, Real VelocityRight,
		Real DensityRight, Real GammaRight, Real SoundRight) {
	int CountIterations = 0;
	Real2 Result;
	Pressure(Result) = (-VelocityRight + VelocityLeft + SoundLeft / GammaLeft
			+ SoundRight / GammaRight) / (SoundLeft
			/ (GammaLeft * PressureLeft) + SoundRight / (GammaRight
			* PressureRight));
	Real IterateVelocityLeft = AdiabatVelocity(PressureLeft, VelocityLeft,
			DensityLeft, SoundLeft, GammaLeft, Pressure(Result), true),
			IterateVelocityRight = AdiabatVelocity(PressureRight,
					VelocityRight, DensityRight, SoundRight, GammaRight,
					Pressure(Result), false);
	while ((fabs(IterateVelocityLeft - IterateVelocityRight) > Epsilon * (fabs(
			IterateVelocityLeft) - fabs(IterateVelocityRight)))
			&& (CountIterations < 1000)) {
		Real DerivatLeft = AdiabatVelocityDerivative(PressureLeft,
				VelocityLeft, DensityLeft, SoundLeft, GammaLeft,
				Pressure(Result), true);
		Real DerivatRight = AdiabatVelocityDerivative(PressureRight,
				VelocityRight, DensityRight, SoundRight, GammaRight,
				Pressure(Result), false);
		Real PressureStep = -(IterateVelocityLeft - IterateVelocityRight)
				/ (DerivatLeft - DerivatRight);
		if (fabs(PressureStep / Pressure(Result)) < Epsilon)
			break;
		Pressure(Result) += PressureStep;
		IterateVelocityRight = AdiabatVelocity(PressureRight, VelocityRight,
				DensityRight, SoundRight, GammaRight, Pressure(Result), false);
		IterateVelocityLeft = AdiabatVelocity(PressureLeft, VelocityLeft,
				DensityLeft, SoundLeft, GammaLeft, Pressure(Result), true);
		CountIterations++;
	};


    #ifndef CUDA
	if (CountIterations >= 1000)
		Velocity(Result) = Real(CountIterations / 0);
    else
	#endif
		Velocity(Result) = Real(0.5) * (IterateVelocityLeft + IterateVelocityRight);
		return Result;
	};

//return Pressure, Velocity, Density
Real4 HalfProblemDepressionWave(Real PressureOut, Real VelocityOut,
		Real DensityOut, Real GammaOut, Real SoundOut, Real PressureContact,
		Real VelocityContact) {
	Real4 Result;
	Real WaveLeftBoundarySpeed = VelocityOut - SoundOut,
			WaveRightBoundarySpeed = VelocityContact * Real(0.5) * (1
					+ GammaOut) - VelocityOut * Real(0.5) * (GammaOut - 1)
					- SoundOut;
	if (WaveLeftBoundarySpeed > 0) {
		Pressure(Result) = PressureOut;
		Velocity(Result) = VelocityOut;
		Density(Result) = DensityOut;
	} else {
		Real SoundContact = SoundOut + Real(0.5) * (GammaOut - 1)
				* (VelocityOut - VelocityContact);
		Real DensityContact = GammaOut * PressureContact / Sqr(SoundContact);
		if (WaveRightBoundarySpeed < 0) {
			Velocity(Result) = VelocityContact;
			Pressure(Result) = PressureContact;
			Density(Result) = DensityContact;
		} else {
			Real Ratio = WaveLeftBoundarySpeed / (WaveLeftBoundarySpeed
					- WaveRightBoundarySpeed);
			Pressure(Result) = PressureOut + (PressureContact - PressureOut)
					* Ratio;
			Density(Result) = DensityOut + (DensityContact - DensityOut)
					* Ratio;
			Velocity(Result) = VelocityOut + (VelocityContact - VelocityOut)
					* Ratio;
		};
	};
	return Result;
}
;

//return Pressure, Velocity, Density
Real4 HalfProblemShockWave(Real PressureOut, Real VelocityOut, Real DensityOut,
		Real GammaOut, Real SoundOut, Real PressureContact,
		Real VelocityContact) {
	Real WaveSpeed;
	Real4 Result;
	Real DPressure = PressureContact / PressureOut;
	if (fabs(1 - DPressure) < Epsilon)
		WaveSpeed = VelocityOut - SoundOut;
	else {
		Real DGamma = (GammaOut - 1) / (2 * GammaOut);
		WaveSpeed = VelocityOut - DGamma * SoundOut * (1 - DPressure) / (1
				- pow(DPressure, DGamma));
	};
	if (WaveSpeed >= 0) {
		Pressure(Result) = PressureOut;
		Velocity(Result) = VelocityOut;
		Density(Result) = DensityOut;
	} else {
		Pressure(Result) = PressureContact;
		Velocity(Result) = VelocityContact;
		Density(Result) = DensityOut * ((GammaOut + 1) * PressureContact
				+ (GammaOut - 1) * PressureOut) / ((GammaOut - 1)
				* PressureContact + (GammaOut + 1) * PressureOut);
	};
	return Result;
}
;

//return Pressure, Velocity, Density
Real4 HalfProblemContactDiscontinuty(Real PressureOut, Real VelocityOut,
		Real DensityOut, Real GammaOut, Real SoundOut, Real PressureContact,
		Real VelocityContact) {
	Real4 Result;
	Pressure(Result) = PressureContact;
	Density(Result) = DensityOut;
	Velocity(Result) = VelocityContact;
	return Result;
}
;

//return Pressure, Velocity, Density, Energy
Real4 SolveHalfProblem(Real PressureOut, Real VelocityOut, Real DensityOut,
		Real GammaOut, Real SoundOut, Real PressureContact,
		Real VelocityContact, bool Left) {
	Real Sign = Left ? Real(1) : Real(-1);
	Real4 Result;
	if (PressureOut > PressureContact)
		Result = HalfProblemDepressionWave(PressureOut, Sign * VelocityOut,
				DensityOut, GammaOut, SoundOut, PressureContact, Sign
						* VelocityContact);
	else if (PressureOut < PressureContact)
		Result = HalfProblemShockWave(PressureOut, Sign * VelocityOut,
				DensityOut, GammaOut, SoundOut, PressureContact, Sign
						* VelocityContact);
	else
		Result = HalfProblemContactDiscontinuty(PressureOut,
				Sign * VelocityOut, DensityOut, GammaOut, SoundOut,
				PressureContact, Sign * VelocityContact);
	Energy(Result) = Pressure(Result) / ((GammaOut - 1) * Density(Result));
	Velocity(Result) = Sign * Velocity(Result);
	return Result;
}
;

//return Pressure, Velocity, Density, Energy
Real4 CaseNotVacuum(Real PressureLeft, Real VelocityLeft, Real DensityLeft,
		Real GammaLeft, Real SoundLeft, Real PressureRight, Real VelocityRight,
		Real DensityRight, Real GammaRight, Real SoundRight) {
	Real2 Contact = NewtonVelocity(PressureLeft, VelocityLeft, DensityLeft,
			GammaLeft, SoundLeft, PressureRight, VelocityRight, DensityRight,
			GammaRight, SoundRight);
	if (Velocity(Contact) > 0)
		return SolveHalfProblem(PressureLeft, VelocityLeft, DensityLeft,
				GammaLeft, SoundLeft, Pressure(Contact), Velocity(Contact),
				true);
	else
		return SolveHalfProblem(PressureRight, VelocityRight, DensityRight,
				GammaRight, SoundRight, Pressure(Contact), Velocity(Contact),
				false);
}
;

//return Pressure, Velocity, Density, Energy
Real4 CaseVacuum(Real PressureLeft, Real VelocityLeft, Real DensityLeft,
		Real GammaLeft, Real SoundLeft, Real PressureRight, Real VelocityRight,
		Real DensityRight, Real GammaRight, Real SoundRight) {
	Real4 Result;
	if (VelocityLeft - SoundLeft >= 0) {
		Density(Result) = DensityLeft;
		Pressure(Result) = PressureLeft;
		Velocity(Result) = VelocityLeft;
		Energy(Result) = Pressure(Result)
				/ ((GammaLeft - 1) * Density(Result));
	} else if (VelocityLeft + 2 * SoundLeft / (GammaLeft - 1) >= 0) {
		Real Ratio = 2 / (GammaLeft + 1);
		Density(Result) = DensityLeft * Ratio;
		Pressure(Result) = PressureLeft * Ratio;
		Velocity(Result) = VelocityLeft + SoundLeft * Ratio;
		Energy(Result) = Pressure(Result)
				/ ((GammaLeft - 1) * Density(Result));
	}
	if (VelocityRight + SoundRight <= 0) {
		Density(Result) = DensityRight;
		Pressure(Result) = PressureRight;
		Velocity(Result) = VelocityRight;
		Energy(Result) = Pressure(Result) / ((GammaRight - 1)
				* Density(Result));
	} else if (VelocityLeft - 2 * SoundRight / (GammaRight - 1) <= 0) {
		Real Ratio = 2 / (GammaRight + 1);
		Density(Result) = DensityRight * Ratio;
		Pressure(Result) = PressureRight * Ratio;
		Velocity(Result) = VelocityRight - SoundLeft * Ratio;
		Energy(Result) = Pressure(Result) / ((GammaRight - 1)
				* Density(Result));
	} else {
		Density(Result) = 0;
		Pressure(Result) = 0;
		Velocity(Result) = 0;
		Energy(Result) = 0;
	};
	return Result;
}
;

Real4 SolveRiemannProblem(Real DensityLeft, Real EnergyLeft, Real VelocityLeft,
		Real GammaLeft, Real SoundLeft, Real DensityRight, Real EnergyRight,
		Real VelocityRight, Real GammaRight, Real SoundRight) {
	Real PressureLeft = DensityLeft * (GammaLeft - 1) * EnergyLeft,
			PressureRight = DensityRight * (GammaRight - 1) * EnergyRight;
	if (PressureLeft < Epsilon && PressureRight < Epsilon)
		return make_Real4(0, 0, 0, 0);
	else if (VelocityLeft + 2 * SoundLeft / (GammaLeft - 1) > VelocityRight - 2
			* SoundRight / (GammaRight - 1))
		return CaseNotVacuum(PressureLeft, VelocityLeft, DensityLeft,
				GammaLeft, SoundLeft, PressureRight, VelocityRight,
				DensityRight, GammaRight, SoundRight);
	else
		return CaseVacuum(PressureLeft, VelocityLeft, DensityLeft, GammaLeft,
				SoundLeft, PressureRight, VelocityRight, DensityRight,
				GammaRight, SoundRight);

}
;
