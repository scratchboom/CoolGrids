//---------------------------------------------------------------------------

#ifndef URiemannH
#define URiemannH
//#include "UUtilits.h"
//#include "UGasDynamicTypes.h"
//---------------------------------------------------------------------------------------------
//extern void SolveRiemannProblem(RGasData *GasData, RGasDynamicVars1D *LeftVars, RGasDynamicVars1D *RightVars, RGasDynamicVars1D *Result);
//extern void SolveRiemannProblem(RGasData *GasData, char *Axis, RGasDynamicVars3D *LeftVars, RGasDynamicVars3D *RightVars, RGasDynamicVars3D *Result);
extern void SolveRiemannProblem( double Density1,
                                double Energy1,        // эпсилон
                                double Velocity1,
                                double Gamma1,
                                double Density2,
                                double Energy2,        // эпсилон
                                double Velocity2,
                                double Gamma2,
                                double &Pressure,
                                double &Density,
                                double &Energy,        // эпсилон
                                double &Velocity);
#endif
