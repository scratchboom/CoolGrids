#pragma once

struct RGB{
	double R;
	double G;
	double B;

	RGB(double r,double g,double b){
		R=r;
		G=g;
		B=b;
	}
};

struct HSV{
	double H;
	double S;
	double V;
};

RGB HSVtoRGB(double H,double S,double V){

	double R,G,B;

	double p,q,t,f;

	if(H < 1.0/6.0){
		f=(H - 0.0/6.0)*6.0;
		p=V*(1.0-S);
		q=V*(1.0-f*S);
		t=V*(1.0-(1.0-f)*S);
		R = V;G = t;B = p;

	}else if(H < 2.0/6.0){
		f=(H - 1.0/6.0)*6.0;
		p=V*(1.0-S);
		q=V*(1.0-f*S);
		t=V*(1.0-(1.0-f)*S);
		R = q;G = V;B = p;

	}else if(H < 3.0/6.0){
		f=(H - 2.0/6.0)*6.0;
		p=V*(1.0-S);
		q=V*(1.0-f*S);
		t=V*(1.0-(1.0-f)*S);
        R = p;G = V;B = t;

	}else if(H < 4.0/6.0){
		f=(H - 3.0/6.0)*6.0;
		p=V*(1.0-S);
		q=V*(1.0-f*S);
		t=V*(1.0-(1.0-f)*S);
		R = p;G = q;B = V;

	}else if(H < 5.0/6.0){
		f=(H - 4.0/6.0)*6.0;
		p=V*(1.0-S);
		q=V*(1.0-f*S);
		t=V*(1.0-(1.0-f)*S);
		R = t;G = p;B = V;
	}else if(H <= 6.0/6.0){
		f=(H - 5.0/6.0)*6.0;
		p=V*(1.0-S);
		q=V*(1.0-f*S);
		t=V*(1.0-(1.0-f)*S);
		R = V;G = p;B = q;
	}

	return RGB(R,G,B);
}

class ColorMap{
	RGB getColor(double value);
};




class RainbowColorMap:ColorMap{

public:

	RGB getColor(double value){
		return HSVtoRGB(value,0.9,0.9);
	}

} rainbowColorMap;



class GrayscaleColorMap:ColorMap{

public:

	RGB getColor(double value){
		return RGB(value,value,value);
	}

} grayscaleColorMap;



class RedBlueColorMap:ColorMap{

public:

	RGB getColor(double value){
		if(value>0.5)return RGB((value-0.5)*2,0,0);
		else return RGB(0,0,(0.5-value)*2);
	}

} redBlueColorMap;



/*
<Point x="0.81779" o="1" r="1" g="1" b="1"/>
  <Point x="1.09565" o="1" r="0.392157" g="0" b="0.00392157"/>
  <Point x="1.38896" o="1" r="1" g="0" b="1"/>
  <Point x="1.6591" o="1" r="0.392157" g="0" b="0.392157"/>
  <Point x="1.92925" o="1" r="0" g="0" b="1"/>
  <Point x="2.21483" o="1" r="0" g="0" b="0.392157"/>
  <Point x="2.50041" o="1" r="0" g="1" b="1"/>
  <Point x="2.79371" o="1" r="0" g="0.392157" b="0.392157"/>
  <Point x="3.05614" o="1" r="0" g="1" b="0"/>
  <Point x="3.36488" o="1" r="0" g="0.392157" b="0"/>
  <Point x="3.61959" o="1" r="1" g="1" b="0"/>
  <Point x="3.79711" o="1" r="0.392157" g="0.392157" b="0"/>
  <Point x="3.92832" o="1" r="1" g="0" b="0"/>
  <Point x="4.05182" o="1" r="0.392157" g="0" b="0"/>
 */
