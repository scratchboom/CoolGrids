#pragma once

double COULUMB_TO_CGS(double q){
	return q/2997924580.0;
}

double KG_TO_G(double m){
	return m*1000.0;
}

double EV_TO_J(double energyInEV){
	return energyInEV * 1.60217646E-19;
}

double J_TO_EV(double energyInJ){
	return energyInJ / 1.60217646E-19;
}

double CELSIUS_TO_KELVIN(double T){
	return T + 273.15;
}



//E = 3/2 kT
double J_TO_KELVIN(double E){
	const double k_b=1.38E-23;// Дж/К TODO переместить в физические константы
	return 2.0/3.0*E / k_b;
}
