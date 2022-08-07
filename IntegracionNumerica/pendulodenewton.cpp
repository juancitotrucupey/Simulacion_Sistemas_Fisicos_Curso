// Péndulo de newton
#include <iostream>
#include <cmath>
#include "Vector.h"
using namespace std;

const double g=980; //gravedad en cgs
const double CHI   =0.1786178958448091;
const double LAMBDA=-0.2123418310626054;
const double XI    =-0.06626458266981849;
const double uno_m2LAMBDA_S2=(1-2*LAMBDA)/2;
const double uno_m2_XIplusCHI=(1-2*(XI+CHI));

const int N=20;
const double k=1e9;

class Cuerpo;
class Colisionador;

//-------------------Clase cuerpo------------------------------
class Cuerpo{
private:
  double theta, omega, tau;
  double I,R,L,x0;
public:
  void Inicie(double theta0,double omega0,double m0,double L0,double R0, double x00);
  void CalculeTorque(void);
  void AdicioneTorque(double F0);
  void Mueva_r(double dt,double CTE);
  void Mueva_V(double dt,double CTE);
  void Dibujese(void);
  double Getx(void){return x0+L*sin(theta);}; //Funcion inline (macro)
  double Gety(void){return -L*cos(theta);}; //Funcion inline (macro)
  double Getz(void){return tau;}; //Funcion inline (macro)

friend class Colisionador;
};

void Cuerpo::Inicie(double theta0,double omega0,double m0,double L0,double R0, double x00){
  theta=theta0;  omega=omega0; L=L0; I=m0*L0*L0;  R=R0;   x0=x00;
}
void Cuerpo::CalculeTorque(void){
	tau=-I/L*g*sin(theta);
	}
void Cuerpo::AdicioneTorque(double F0){
	tau+=L*F0;
	}
void Cuerpo::Mueva_r(double dt,double CTE){
  theta+=omega*(CTE*dt);
}
void Cuerpo::Mueva_V(double dt,double CTE){
  omega+=tau*(CTE*dt/I);
}
void Cuerpo::Dibujese(void){
  cout<<", "<<Getx()<<"+"<<R<<"*cos(t),"<<Gety()<<"+"<<R<<"*sin(t)";
}


//-------------Clase colisionador--------------
class Colisionador{
private:
public:  
    void CalculeTodasLasFuerzas(Cuerpo*Pendulo);
    void CalculeFuerzaEntre(Cuerpo & Pendulo1, Cuerpo & Pendulo2);  
};

void Colisionador::CalculeTodasLasFuerzas(Cuerpo*Pendulo){
	int i,j;
	//Calcular el torque hecho por la gravedad sobre cada pendulo
	for(i=0;i<N;i++)
		Pendulo[i].CalculeTorque();
	//Cada péndulo (menos el último) calcula el choque con el siguiente
	for(i=0;i<N-1;i++)
			CalculeFuerzaEntre(Pendulo[i],Pendulo[i+1]);
}
	
	
void Colisionador::CalculeFuerzaEntre(Cuerpo & Pendulo1, Cuerpo & Pendulo2){ //calcula fuerza de hertz
	double s,x1,x2,F2;
	x1=Pendulo1.Getx();   x2=Pendulo2.Getx();
	s=(Pendulo2.R+Pendulo1.R)-fabs(x2-x1);
	
	if(s>0){  //hay choque
		double F=k*pow(s,1.5);
		if(x2>x1)   F2=F;   else F2=-F;
		Pendulo2.AdicioneTorque(F2);    Pendulo1.AdicioneTorque(-F2);
	}

}



//-------------------------------------
void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'MiBalon.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange [-4:42]"<<endl;
  cout<<"set yrange [-12:0]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
}
void TermineCuadro(void){
    cout<<endl;
}


int main(void){
  Cuerpo Pendulo[N];   int i;
  Colisionador GravedadDeNewton;
  double t,Deltat=1e-6;
  double L=10,m0=20,R=1; 
  double omega,T, tmax,tdibujo;
  int ndibujos=1000;
  
  omega=sqrt(g/L);   T=2*M_PI/omega;
  
  tmax=5*T;



  InicieAnimacion();

  //------------(theta0, omega0, m0, L0, R0, x00);
  for(i=0;i<N-1;i++)
     Pendulo[i].Inicie(0,0,m0,L,R,i*2*R);
     
  Pendulo[N-1].Inicie(15*M_PI/180, 0 , m0,  L, R, 2*(N-1)*R);
  for(t=0,tdibujo=0;t<tmax;t+=Deltat,tdibujo+=Deltat){
    //Animacion
    
    if(tdibujo>tmax/ndibujos){
      InicieCuadro();
      for(i=0;i<N;i++)  Pendulo[i].Dibujese();
      TermineCuadro();
      tdibujo=0;
    }
    
    for(i=0;i<N;i++)  Pendulo[i].Mueva_r(Deltat,CHI);
    GravedadDeNewton.CalculeTodasLasFuerzas(Pendulo); for(i=0;i<N;i++) Pendulo[i].Mueva_V(Deltat,uno_m2LAMBDA_S2);
    for(i=0;i<N;i++) Pendulo[i].Mueva_r(Deltat,XI);
    GravedadDeNewton.CalculeTodasLasFuerzas(Pendulo); for(i=0;i<N;i++) Pendulo[i].Mueva_V(Deltat,LAMBDA); 
    for(i=0;i<N;i++) Pendulo[i].Mueva_r(Deltat,uno_m2_XIplusCHI);
    GravedadDeNewton.CalculeTodasLasFuerzas(Pendulo); for(i=0;i<N;i++) Pendulo[i].Mueva_V(Deltat,LAMBDA); 
    for(i=0;i<N;i++) Pendulo[i].Mueva_r(Deltat,XI);
    GravedadDeNewton.CalculeTodasLasFuerzas(Pendulo); for(i=0;i<N;i++) Pendulo[i].Mueva_V(Deltat,uno_m2LAMBDA_S2);
    for(i=0;i<N;i++) Pendulo[i].Mueva_r(Deltat,CHI);

  }

  return 0;
}
