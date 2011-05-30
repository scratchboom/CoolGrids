
#include "GridsCommon.hpp"

const double PI_A0_SQR = 0.88E-16     *1E-4;// [м^2]
const double Ry=13.6;// [эВ]

const double I = 14.86;// Константа I в ф-ле для сечения ионизации, эВ

const double MU = 1.0;


//Глобальные физические константы
const double EPS0=8.854E-12;// Ф/м
const double MU0=1.256E-6;// Гн/м
const double SPEED_OF_LIGHT=299792458.0;// [м/с]
const double VACUUM_WAVE_RESISTIVITY=sqrt(MU0/EPS0);// волновое сопротивление вакуума [Ом]

const double C=SPEED_OF_LIGHT;// м/с    C=sqrt(1/EPS0/MU0)

const double ATOMIC_MASS_UNIT = 1.660538782E-27;// [кг]
const double AVERAGE_AIR_ION_MASS = 29.0 * ATOMIC_MASS_UNIT;// [кг]  для высот < 100-120 км молекулярная масса почти постоянна
const double ELECTRON_MASS = 9.10938215E-31;// Масса электрона, [кг]
const double ELECTRON_CHARGE = 1.602E-19;// Заряд электрона, [Кл]


//параметры атмосферы на высоте 40км
const double ATMOSPHERE_O2_CONCENTRATION_40KM = 1.7E16         * 1E6;// [1/м^3]
const double ATMOSPHERE_ELECTRON_CONCENTRATION_40KM = 3.16E4   * 1E6;// [1/м^3]
const double ATMOSPHERE_TEMPERATURE_40KM = 250.0;// [К]
const double ATMOSPHERE_DENSITY_40KM = 0.003851;// [кг/м^3]


const double Wevlt = 1.602176487E-19;// Джоулей в одном э­лектрон-Вольте

const double cs = SPEED_OF_LIGHT;//Скороcть cвета, m/s [cs]

const double cze = 1.602E-19;//Заряд электрона, Кл. [cze]
const double cme = 9.10938215E-31;//Масса электрона, кг ? [cme]
const double ckb = 0.861E-4;//Константа Больцмана, э­В/Кград [ckb]
const double ci = 14.86;//Константа I в ф-ле для сечения ионизации, эВ [ci]
const double cih = 13.6;//Консттанта Ih в ф-ле для сечения ионизации, эВ [cih]
const double ca2 = 3.52E-16;//Константа 4*Pi*alfa2 в ф-ле для сечения ионизации, sm2 [ca2]
const double ce0 = 0.025;//Константа e0 в ф-ле для Qei, эВ [ce0]
const double cnw = 1.14e23;//Концентрация воздуха у поверхности Земли, 1/м3 ? [cnw]
const double cbv = -0.348E-4;//Показатель компоненты в распр. концентрации воздуха по высоте, ? [cbv]

const double cno2 = 0.21;//Доля кислорода в атмосфере, ? [cno2]

const double czi = 900000.0;//Высота начала ионосферы, м ? [czi]   60-90км

const double czb = 1.0;//Степень ионизации [czb]

const double cne0 = 3E4;//28E6;//начальная концентрация электронов, 1/см3 ? [cne0]
const double cni0 = cni0;//28E6;//начальная концентрация ионов, 1/см3 ? [cni0]

//!!! НЕ ИСПОЛЬЗУЕТСЯ! double ct0=200000.0;// ЖЕСТКОВАТО!!! 1000. 115070.   Температура воздуха у поверхности Земли, гр.К [ct0]
//!!!double cti=200000.0;//  1000. 115070.   наачальная температура электронов, гр.К [cti] TODO какова на самом деле температура электронов в атмосфере
const double chi = 4E9;//  8 5.5e8 !!!!!0.55e9 !!!3.e9          Опорная частота (Гц) [chi]
const double csi = 4E5;//  1.e8 1.5e8 !!!1.e6    Частота следований импульсов (скважность) (Гц) [csi]
const double cpi = 1.0;//  25 50 5  !!!!!3 1 !!!50.                Число импульсов в "пачке" [cpi]
const double cep = 1E7;//  !!!5.e5 !!!1.e7 1.e4       амплитуда напряженности электрического поля (В/м) [cep]
const double Z0 = 42000.0;//  40000. !!!50000. Nachalnaya vysota   176121372031.662 [Z0]

const double WPI=M_PI;
const double WZ = 1.0;

const double WM = AVERAGE_AIR_ION_MASS;//= 4.84d-26 ! average air ion mass = ??????? in SI
const double WEt0 = 0.025;
const double Wmu0 = 4*WPI*1E-7;


double k_b=1.38E-23;// Дж/К

double NE0=ATMOSPHERE_ELECTRON_CONCENTRATION_40KM;
double ATM_T=200.0; //температура атмоcферы K


struct UserData{
	double E_x;
	double E_y;
	double E_z;

	double H_x;
	double H_y;
	double H_z;

	double V_SQR;
};

void rhs(UserData* userData , int *ptrN,double *t,double *u,double *f, int i);

void rhs2(UserData* userData , int *ptrN,double *t,double *u,double *f, int i){
	f[0] = (*t);
	f[1] = (*t)*(*t);
	f[2] = (*t)*(*t)*(*t);
	f[3] = sin((*t));
	f[4] = cos((*t));
}



Grid1D<double> gridM_N_V_ei_Vx("-(v_ei     )*V_x");
Grid1D<double> gridM_N_V_e0_Vx("-(     v_e0)*V_x");
Grid1D<double> grid_e_Ex_m("-e*E_x/m");
Grid1D<double> grid_e_MU_MU0_m_Vz_Hy("e*MU*MU0/m*V_z*H_y");
Grid1D<double> grid_e_MU_MU0_m_Vy_Hz("-e*MU*MU0/m*V_y*H_z");



//Grid1D<double> gridM_N_V_ei_Vx;
//Grid1D<double> gridMNV_e0_Vx;

int NT = 10000;
double dt = 1E-16;

//int NT = 100;
//double dt = 1E-2;


int main(){



	double E_AMPLITUDE = 1E6;
	double H_AMPLITUDE = E_AMPLITUDE*VACUUM_WAVE_RESISTIVITY;

	Grid1D<double> gridN;
	gridN.setIndexRangeX(0,NT);
	gridN.setRangeX(0,1);
	gridN.build();
	gridN.fill(0);
	DBGLN;

	DBGLN;


	gridM_N_V_ei_Vx.setIndexRangeX(0,NT);
	gridM_N_V_ei_Vx.setRangeX(0,1);
	gridM_N_V_ei_Vx.build();
	gridM_N_V_ei_Vx.fill(0);

	gridM_N_V_e0_Vx.setIndexRangeX(0,NT);
	gridM_N_V_e0_Vx.setRangeX(0,1);
	gridM_N_V_e0_Vx.build();
	gridM_N_V_e0_Vx.fill(0);

	grid_e_Ex_m.setIndexRangeX(0,NT);
	grid_e_Ex_m.setRangeX(0,1);
	grid_e_Ex_m.build();
	grid_e_Ex_m.fill(0);

	grid_e_MU_MU0_m_Vz_Hy.setIndexRangeX(0,NT);
	grid_e_MU_MU0_m_Vz_Hy.setRangeX(0,1);
	grid_e_MU_MU0_m_Vz_Hy.build();
	grid_e_MU_MU0_m_Vz_Hy.fill(0);

	grid_e_MU_MU0_m_Vy_Hz.setIndexRangeX(0,NT);
	grid_e_MU_MU0_m_Vy_Hz.setRangeX(0,1);
	grid_e_MU_MU0_m_Vy_Hz.build();
	grid_e_MU_MU0_m_Vy_Hz.fill(0);





	Grid1D<double> gridVx;
	gridVx.setIndexRangeX(0,NT);
	gridVx.setRangeX(0,1);
	DBGLN;
	gridVx.build();
	gridVx.fill(0);
	DBGLN;



	Grid1D<double> gridVy;
	gridVy.setIndexRangeX(0, NT);
	gridVy.setRangeX(0, 1);
	DBGLN;
	gridVy.build();
	gridVy.fill(0);

	DBGLN;

	Grid1D<double> gridVz;
	gridVz.setIndexRangeX(0, NT);
	gridVz.setRangeX(0, 1);
	gridVz.build();
	gridVz.fill(0);

	Grid1D<double> gridT;
	gridT.setIndexRangeX(0, NT);
	gridT.setRangeX(0, 1);
	gridT.build();
	gridT.fill(0);



	GnuPlotSaver1D gnuPlotSaverN;
	gnuPlotSaverN.setLineColor("#FF0000");
	//gnuPlotSaverN.setValueRange(-2.0*H_AMPLITUDE , 2.0*H_AMPLITUDE)
	gnuPlotSaverN.setAutoScale();
	gnuPlotSaverN.setSize(600,250);

	GnuPlotSaver1D gnuPlotSaverVx;
	gnuPlotSaverVx.setLineColor("#FF0000");
		//gnuPlotSaverN.setValueRange(-2.0*H_AMPLITUDE , 2.0*H_AMPLITUDE)
	gnuPlotSaverVx.setAutoScale();
	gnuPlotSaverVx.setSize(600,250);

	DBGLN;
	GnuPlotSaver1D gnuPlotSaverVy;
	gnuPlotSaverVy.setLineColor("#FF0000");
		//gnuPlotSaverN.setValueRange(-2.0*H_AMPLITUDE , 2.0*H_AMPLITUDE)
	gnuPlotSaverVy.setAutoScale();
	gnuPlotSaverVy.setSize(600,250);


	GnuPlotSaver1D gnuPlotSaverVz;
	gnuPlotSaverVz.setLineColor("#FF0000");
	//gnuPlotSaverN.setValueRange(-2.0*H_AMPLITUDE , 2.0*H_AMPLITUDE)
	gnuPlotSaverVz.setAutoScale();
	gnuPlotSaverVz.setSize(600, 250);

	GnuPlotSaver1D gnuPlotSaverT;
	gnuPlotSaverT.setLineColor("#FF0000");
	//gnuPlotSaverN.setValueRange(-2.0*H_AMPLITUDE , 2.0*H_AMPLITUDE)
	gnuPlotSaverT.setAutoScale();
	gnuPlotSaverT.setSize(600, 250);
	DBGLN;


	Grid1D<double> gridF0;
	gridF0.setIndexRangeX(0, NT);
	gridF0.setRangeX(0, 1);
	gridF0.build();
	gridF0.fill(0);

	Grid1D<double> gridF1;
	gridF1.setIndexRangeX(0, NT);
	gridF1.setRangeX(0, 1);
	gridF1.build();
	gridF1.fill(0);

	Grid1D<double> gridF2;
	gridF2.setIndexRangeX(0, NT);
	gridF2.setRangeX(0, 1);
	gridF2.build();
	gridF2.fill(0);

	Grid1D<double> gridF3;
	gridF3.setIndexRangeX(0, NT);
	gridF3.setRangeX(0, 1);
	gridF3.build();
	gridF3.fill(0);

	Grid1D<double> gridF4;
	gridF4.setIndexRangeX(0, NT);
	gridF4.setRangeX(0, 1);
	gridF4.build();
	gridF4.fill(0);


	gridM_N_V_ei_Vx.setIndexRangeX(0, NT);
	gridM_N_V_ei_Vx.setRangeX(0, 1);
	gridM_N_V_ei_Vx.build();
	gridM_N_V_ei_Vx.fill(0);

	gridM_N_V_e0_Vx.setIndexRangeX(0, NT);
	gridM_N_V_e0_Vx.setRangeX(0, 1);
	gridM_N_V_e0_Vx.build();
	gridM_N_V_e0_Vx.fill(0);


	GnuPlotSaver1D gnuPlotSaverAuto;
	gnuPlotSaverAuto.setLineColor("#FF0000");
	//gnuPlotSaverN.setValueRange(-2.0*H_AMPLITUDE , 2.0*H_AMPLITUDE)
	gnuPlotSaverAuto.setAutoScale();
	gnuPlotSaverAuto.setSize(1000, 250);
	DBGLN;






	int N = 5;//размерность системы уравнений
	double u[N];
	double f[N];
	DBGLN;

	double kT = ckb * ATM_T;


	u[0] = ATMOSPHERE_ELECTRON_CONCENTRATION_40KM;
	u[1] = 0;
	u[2] = 0;
	u[3] = 0;
	u[4] = EV_TO_J(3.0/2.0*kT);

    /*
	u[0] = 0;
	u[1] = 0;
	u[2] = 0;
	u[3] = 0;
	u[4] = 0;
*/

	UserData ud;
	ud.E_x = 0.0;
	ud.E_y = E_AMPLITUDE;
	ud.E_z = 0.0;

	ud.H_x = 0.0;
	ud.H_y = 0.0;
	ud.H_z = H_AMPLITUDE;

	ud.V_SQR = 0.0;

	double t=0;

	DBGLN;

	double tmp[5];

	double y0[5];
	double y1[5];
	double y2[5];
	double y3[5];

	double k1[5];
	double k2[5];
	double k3[5];
	double k4[5];

	for(int i=0;i<NT-1;i++){


		t = i *dt;

		//euler
		//rhs2(&ud,&N,&t,u,f,i);
		//for(int j=0;j<5;j++) u[j]+=f[j]*dt;



		//RK4
		for (int j = 0; j < 5; j++)	y0[j] = u[j];

		double tt = t;
		rhs(&ud, &N, &tt, y0, k1, i);
		for (int j = 0; j < 5; j++)	y1[j] = y0[j] + k1[j] * dt / 2.0;

		tt = t + dt / 2.0;
		rhs(&ud, &N, &tt, y1, k2, i);
		for (int j = 0; j < 5; j++)	y2[j] = y0[j] + k2[j] * dt / 2.0;

		tt = t + dt / 2.0;
		rhs(&ud, &N, &tt, y2, k3, i);
		for (int j = 0; j < 5; j++) y3[j] = y0[j] + k3[j] * dt;

		tt = t + dt;
		rhs(&ud, &N, &tt, y3, k4, i);
		for (int j = 0; j < 5; j++)	u[j] = y0[j] + (k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j]) * dt / 6.0;




		gridF0(i) = f[0];
		gridF1(i) = f[1];
		gridF2(i) = f[2];
		gridF3(i) = f[3];
		gridF4(i) = f[4];

		gridN(i) = u[0];
		gridVx(i) = u[1];
		gridVy(i) = u[2];
		gridVz(i) = u[3];
		gridT(i) = u[4];
	}

	gnuPlotSaverT.save(gridT,"T.plot");
	gnuPlotSaverN.save(gridN,"N.plot");
	gnuPlotSaverVx.save(gridVx,"Vx.plot");
	gnuPlotSaverVy.save(gridVy,"Vy.plot");
	gnuPlotSaverVz.save(gridVz,"Vz.plot");

	gnuPlotSaverAuto.save(gridF0,"F0.plot");
	gnuPlotSaverAuto.save(gridF1,"F1.plot");
	gnuPlotSaverAuto.save(gridF2,"F2.plot");
	gnuPlotSaverAuto.save(gridF3,"F3.plot");
    gnuPlotSaverAuto.save(gridF4,"F4.plot");


    gnuPlotSaverAuto.save(gridM_N_V_ei_Vx,"gridM_N_V_ei_Vx.plot");
	gnuPlotSaverAuto.save(gridM_N_V_e0_Vx,"gridM_N_V_e0_Vx.plot");
	gnuPlotSaverAuto.save(grid_e_Ex_m,"grid_e_Ex_m.plot");
	gnuPlotSaverAuto.save(grid_e_MU_MU0_m_Vz_Hy,"grid_e_MU_MU0_m_Vz_Hy.plot");
	gnuPlotSaverAuto.save(grid_e_MU_MU0_m_Vy_Hz,"grid_e_MU_MU0_m_Vy_Hz.plot");


    DBGVAL("_________");

    DBGVAL(u[0]);
    DBGVAL(u[1]);
    DBGVAL(u[2]);
    DBGVAL(u[3]);
    DBGVAL(u[4]);
    DBGVAL(J_TO_EV(u[4]));

	return 0;
}






void rhs(UserData* userdata , int *ptrN,double *t,double *u,double *f,int i){



	int N = *ptrN;

	double n_e = u[0];

	double V_x = u[1];
	double V_y = u[2];
	double V_z = u[3];

	double Et = u[4];//тепловая энергия одного электрона [Дж]
	double WEt = J_TO_EV(Et);//тепловая энергия одного электрона [эВ]
	double WEfull = WEt + ELECTRON_MASS*sqr(V_x,V_y,V_z)/2.0;//полная энергия одного электрона [эВ]

//	DBGVAL(WEt);

	const double n_i = n_e; //концентрация положительных ионов [м^-3] TODO why n_i ~= n_e ?

    const double n_all = ATMOSPHERE_DENSITY_40KM / AVERAGE_AIR_ION_MASS;
	const double n_O2 = ATMOSPHERE_O2_CONCENTRATION_40KM;
	const double n = n_all - n_O2; //концентрация нейтралов + положительных ионов TODO[] отсылка к Ступицкому
	const double n_0 = n - n_i;//т.к n = n_i + n_0 (такие обозначения)


	//const UserData* userdata = (UserData*)(&u[N]);
	//const UserData userdata = (UserData*)(&u[N]);

    const double E_x = userdata->E_x;
    const double E_y = userdata->E_y;
    const double E_z = userdata->E_z;

    const double H_x = userdata->H_x;
    const double H_y = userdata->H_y;
    const double H_z = userdata->H_z;

    //const double V_SQR = userdata->V_SQR;
    const double V_SQR = sqr(u[1],u[2],u[3]);

    const double m_e = ELECTRON_MASS;
    const double m = m_e;
    const double M = AVERAGE_AIR_ION_MASS;
    const double e = ELECTRON_CHARGE;


    const double f_ = 11.67*(1.0+0.64*(WEt-11.35)/(88.65+WEt))*(1.0-exp(-0.0083*(WEt-11.35)));// [безразмерно]
    const double sgm =  WEt>ci  ?  f_*4.0*PI_A0_SQR* sqr(Ry/WEt)*(WEt-ci)/ci :  0;// [m^2] TODO: Порог (WEt должно быть > ci) Se=0 или sgm=0
    const double V = 5.3E7*sqrt(WEt)                          * 1E-2;//[m/s]
//    DBGVAL(f_);
//    DBGVAL(sgm);
//    DBGVAL(V);

    const double j_0e = sgm*V;//[m^3/s]
    const double j_g = 1.16E-8 / WEt                          *1E-6;// [m^3/s]
    const double j_v = 2.7E-13/pow(2.0/3.0*WEt,3.0/4.0)       *1E-6;// [m^3/s]

    const double j_ei = 8.75E-27/ pow((2.0/3.0*WEt),9.0/2.0)  *1E-12;// [m^6/s]
    const double j_p = 3.8E-31/WEt * exp(-0.103/WEt)          *1E-12;// [m^6/s]

    const double kT = 2.0/3.0*WEt;//температура электронов в эВ [эВ]


    const double L = 25.2 + log(2.0*WEt/3.0/ckb)- 0.5*log(n_e        * 1E-6);//TODO check for range 15..18; decimal or natural logorythm
//    DBGVAL(L);
    const double Z =1.0;

    const double sigma_e0 = 12.47 * PI_A0_SQR * (0.4 + 0.84*WEt/(0.5 + WEt));// [м^2]

    const double v_ei = 16.0*sqrt(M_PI)/3.0 * quad( COULUMB_TO_CGS(cze) ) * L * n_i * sqr(Z)/sqrt(m_e)/pow(2.0*kT,3.0/2.0);
    const double v_e0 = 8.0*sigma_e0/3.0/sqrt(M_PI)*sqrt(2.0* EV_TO_J(kT)/m_e)*n_0;


    const double S_e = WEt>ci  ?  (j_0e*n_e*n_0 - j_ei*sqr(n_e)*n_i) - j_v*n_e*n_i - j_g*n_e*n_i - j_p*n_e*n_O2*n   :  0.0;//TODO [?] пороги!!!

    const double fi = 0.64 + 0.11*log(I/kT);//TODO natural or decimal
    const double S_ee = -(I+3.0/2.0*kT)*Wevlt*(n_e*n_0*j_0e - sqr(n_e)*n_i*j_ei)+
    		            + (3.0/2.0-fi)*kT*Wevlt*n_e*n_i*j_v - 3.0/2.0*kT*Wevlt*j_g*n_e*n_i - 3.0/2.0*kT*Wevlt*j_p*n_e*n_O2*n;//TODO заменить Wevlt на EV_TO_J()

    const double WEt_ = 0.025;//TODO что это? сколько это градусов?
    const double Q_ei = -2.0 * v_ei*n_e*(WEt-WEt_)*Wevlt*m_e/M + m*V_SQR/2.0;//TODO заменить Wevlt на EV_TO_J() ? порог??? WEt должно быть > WEt_

    const double nu = 6E-8 * n_0*sqrt(WEt)*(0.4 + 0.84*WEt/(0.5+WEt))      *1E-6;// [s^-1]   1E6-коэффициент перевода для концентрации из СИ в СГС (так как оригинальная формула написана в СГС)
    const double delta = 1.7E-3 * (1.0 + 0.2*pow(WEt/0.9,5)) / (1.0 + 3.7E-2*(1.0 + 0.2*(pow(WEt/0.9,5))) );// безразмерно

    const double Q_e0 = -n_e*EV_TO_J(WEt-WEt_)*nu*delta;//TODO см диплом  порог

    const double Q_e = Q_ei + Q_e0;
    //const double Q_w = e*n_e*(abs(E_x*V_x) + abs(E_y*V_y) + abs(E_z*V_z));//TODO не симметрично относительно поворота СК, что имел ввиду Ступитский?
    //const double Q_w = e*n_e*(E_x*V_x + E_y*V_y + E_z*V_z);//TODO не симметрично относительно поворота СК, размерность нормальная? [Дж/с/м^3]
    //const double Q_w = e*n_e*sqrt( sqr(E_x*V_x , E_y*V_y , E_z*V_z) );
    const double Q_w = e*n_e*( E_x*V_x + E_y*V_y + E_z*V_z );
    DBGVAL(Q_w);
    DBGVAL(e*n_e*E_x*V_x  /n_e);
    DBGVAL(e*n_e*E_y*V_y  /n_e);
    DBGVAL(e*n_e*E_z*V_z  /n_e);





	f[0] = S_e;
    f[1] = -(v_ei+v_e0)*V_x  -  e*E_x/m  +  e*MU*MU0/m*V_z*H_y - e*MU*MU0/m*V_y*H_z;
    f[2] = -(v_ei+v_e0)*V_y  -  e*E_y/m  +  e*MU*MU0/m*V_x*H_z - e*MU*MU0/m*V_z*H_x;
    f[3] = -(v_ei+v_e0)*V_z  -  e*E_z/m  +  e*MU*MU0/m*V_y*H_x - e*MU*MU0/m*V_x*H_y;
    f[4] = (S_ee + Q_e + Q_w)/n_e;



    gridM_N_V_ei_Vx(i) = -(v_ei     ) ;//*V_x;
	gridM_N_V_e0_Vx(i) = -(     v_e0);//*V_x;
	grid_e_Ex_m(i) = -e*E_y/m;
	grid_e_MU_MU0_m_Vz_Hy(i) = e*MU*MU0/m*V_z*H_y;
	grid_e_MU_MU0_m_Vy_Hz(i) = -e*MU*MU0/m*V_y*H_z;



    DBGVAL(-(v_ei )*V_y);
    DBGVAL(-(     v_e0)*V_y);
    DBGVAL(-  e*E_y/m);
    DBGVAL(-E_y);
    DBGVAL(e);
    DBGVAL(m);



    #if(0)
    if((isnan(u[0])==0) && (isinf(u[0])==0)) acc_U0(u[0]);
    if((isnan(u[1])==0) && (isinf(u[1])==0)) acc_U1(u[1]);
    if((isnan(u[2])==0) && (isinf(u[2])==0)) acc_U2(u[2]);
    if((isnan(u[3])==0) && (isinf(u[3])==0)) acc_U3(u[3]);
    if((isnan(u[4])==0) && (isinf(u[4])==0)) acc_U4(u[4]);

    if((isnan(f[0])==0) && (isinf(f[0])==0)) acc_f0(f[0]);
    if((isnan(f[1])==0) && (isinf(f[1])==0)) acc_f1(f[1]);
    if((isnan(f[2])==0) && (isinf(f[2])==0)) acc_f2(f[2]);
    if((isnan(f[3])==0) && (isinf(f[3])==0)) acc_f3(f[3]);
    if((isnan(f[4])==0) && (isinf(f[4])==0)) acc_f4(f[4]);
    #endif

    #if(1)
    DBGVAL(dt / *t);

    DBGVAL(j_0e );
    DBGVAL(n_e);
    DBGVAL(n_0);

    DBGVAL(- j_ei*sqr(n_e)*n_i);
    DBGVAL(j_0e*n_e*n_0);
    DBGVAL(- j_v*n_e*n_i);
    DBGVAL(- j_g*n_e*n_i);
    DBGVAL(- j_p*n_e*n_O2*n);


    DBGVAL(v_ei);
    DBGVAL(v_e0);
    DBGVAL(-(v_ei+v_e0)*V_x);
    DBGVAL(-  e*E_x/m);
    DBGVAL(e*MU*MU0/m*V_z*H_y);
    DBGVAL(- e*MU*MU0/m*V_z*H_x);











    DBGVAL(-(I+3.0/2.0*kT)*Wevlt*(n_e*n_0*j_0e  ));
    DBGVAL(-(I+3.0/2.0*kT)*Wevlt*(             - sqr(n_e)*n_i*j_ei));
    DBGVAL((3.0/2.0-fi)*kT*Wevlt*n_e*n_i*j_v);
    DBGVAL(-3.0/2.0*kT*Wevlt*j_g*n_e*n_i);
    DBGVAL(-3.0/2.0*kT*Wevlt*j_p*n_e*n_O2*n);//TODO заменить Wevlt на EV_TO_J()


    DBGVAL(S_ee/n_e);
    DBGVAL(Q_e/n_e);
    DBGVAL(Q_w/n_e);

    DBGVAL(Q_ei/n_e);
    DBGVAL(Q_e0/n_e);

    DBGVAL(u[0]);
    DBGVAL(u[1]);
    DBGVAL(u[2]);
    DBGVAL(u[3]);
    DBGVAL(u[4]);
    DBGVAL(J_TO_EV(u[4]));

    DBGVAL(f[0]);
    DBGVAL(f[1]);
    DBGVAL(f[2]);
    DBGVAL(f[3]);
    DBGVAL(f[4]);

    DBGVAL(f[0]*dt);
    DBGVAL(f[1]*dt);
    DBGVAL(f[2]*dt);
    DBGVAL(f[3]*dt);
    DBGVAL(f[4]*dt);

    DBGVAL(dt);

//    pressAnyKey();





//    DBGVAL(-(v_ei+v_e0)*V_x);
//    DBGVAL(-e*E_x/m);
//    DBGVAL(e*MU*MU0/m*V_z*H_y);
//    DBGVAL(- e*MU*MU0/m*V_y*H_z);
//
//
//    DBGVAL("  -----  ");
//
//    DBGVAL(-(v_ei+v_e0)*V_x);
//    DBGVAL(-e*E_x/m);
//    DBGVAL(e*MU*MU0/m*V_z*H_y);
//    DBGVAL(- e*MU*MU0/m*V_y*H_z);
//
//    DBGVAL("  -----  ");
//
//    DBGVAL(-(v_ei+v_e0)*V_y);
//    DBGVAL(-  e*E_y/m);
//    DBGVAL(e*MU*MU0/m*V_x*H_z);
//    DBGVAL(- e*MU*MU0/m*V_z*H_x);
//
//    DBGVAL("  -----  ");
//
//    DBGVAL(-(v_ei+v_e0)*V_z);
//    DBGVAL(-  e*E_z/m);
//    DBGVAL(e*MU*MU0/m*V_y*H_x);
//    DBGVAL(- e*MU*MU0/m*V_x*H_y);


    #endif


}
