// UnPlaneta[i] por LeapFrog
#include <iostream>
#include <cmath>
#include "Vector.h"
using namespace std;

const double G=1.0;
const double CHI   =0.1786178958448091;
const double LAMBDA=-0.2123418310626054;
const double XI    =-0.06626458266981849;
const double uno_m2LAMBDA_S2=(1-2*LAMBDA)/2;
const double uno_m2_XIplusCHI=(1-2*(XI+CHI));

const int N=2; //número de Planetas
//............................
class Cuerpo;
class Colisionador;
//............................
class Cuerpo{
private:
  vector3D r, V, F;
  double m,R;
public:
  void Inicie(double x0,double y0, double z0, double Vx0,double Vy0, double Vz0, double m0,double R0);
  void CalculeFuerza(void);
  void BorreFuerza(void);
  void AdicioneFuerza(vector3D F0);
  void Mueva_r(double dt,double CTE);
  void Mueva_V(double dt,double CTE);
  void Dibujese(void);
  double Getx(void){return r.x();}; //Funcion inline (macro)
  double Gety(void){return r.y();}; //Funcion inline (macro)
  
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0, double m0,double R0){
  r.cargue(x0,y0,z0); V.cargue(Vx0,Vy0,Vz0); m=m0;  R=R0;
}
void Cuerpo::CalculeFuerza(void){
  double aux=-G*m*pow(norma2(r),-1.5);
  F=r*aux;
}
void Cuerpo::BorreFuerza(void){
	F.cargue(0,0,0);
}

void Cuerpo::AdicioneFuerza(vector3D F0){
	F+=F0;
}
void Cuerpo::Mueva_r(double dt,double CTE){
  r+=V*CTE*dt;
}
void Cuerpo::Mueva_V(double dt,double CTE){
  V+=F*(CTE*dt/m);
}
void Cuerpo::Dibujese(void){
  cout<<", "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}
//-------------------------------------
void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'MiBalon.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange [-110:110]"<<endl;
  cout<<"set yrange [-110:110]"<<endl;
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
//.-.--.-.-.-.-.-.-.-.-.CLASE COLISIONADOR.-.-.-.-.-.-.-.-.-.-.--.-.-.-.
class Colisionador{
	private:
	public:
		void CalculeFuerzaEntre(Cuerpo & Planeta1, Cuerpo & Planeta2);
		void CalculeTodasLasFuerzas (Cuerpo * Planeta);
};

void Colisionador::CalculeTodasLasFuerzas (Cuerpo * Planeta){
	int i, j;
	//Borrar las fuerzas de todos los planetas
	for(i=0; i<N; i++)
		Planeta[i].BorreFuerza();
		//Hacer todas las parejas entre planetas y calcular la fuerza en cada pareja
		for(i=0; i<N; i++)
			for(j=0;j<i;j++)
				CalculeFuerzaEntre(Planeta[i],Planeta[j]);		
}

void Colisionador::CalculeFuerzaEntre(Cuerpo & Planeta1, Cuerpo & Planeta2){
	vector3D dr=Planeta2.r-Planeta1.r;
	double aux=-G*Planeta1.m*Planeta2.m*pow(norma2(dr), -1.5);
	vector3D F2=dr*aux;
	Planeta2.AdicioneFuerza(F2); Planeta1.AdicioneFuerza(F2*(-1));
}

int main(void){
  Cuerpo Planeta[N];
  Colisionador GravedadDeNewton;
  double t,dt=0.1;
  double r=100, m0=1047, m1=1;
  double Omega, M, T, tmax;
  double tdibujo, ndibujos=10000;
  double x0, x1, V0, V1;
  int i;

  M=m0+m1; Omega=sqrt(G*M/(r*r*r)); T=2*M_PI/Omega;
  x0=-r*m1/M; x1=r*m0/M; V0=Omega*x0; V1=Omega*x1;
    
  tmax=20.5*T;
  
  //   InicieAnimacion();
  
  Planeta[0].Inicie(x0,0,0,0,V0,0,m0,20);//Sol
  Planeta[1].Inicie(x1,0,0,0,V1,0,m1,10);//Jupiter

  
  //------------(x0,y0,Vx0,Vy0,m0   ,R0);
	for(t=0,tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
	   cout<<Planeta[0].Getx()<<" "<<Planeta[0].Gety()<<endl;
		//Animacion

	  /*	  	       
		if(tdibujo>tmax/ndibujos){
			InicieCuadro();
			for(i=0; i<N; i++) Planeta[i].Dibujese();
			TermineCuadro();
			tdibujo=0;
			}
	  */
		for(i=0; i<N; i++) Planeta[i].Mueva_r(dt,CHI);
		GravedadDeNewton.CalculeTodasLasFuerzas(Planeta); for(i=0; i<N; i++) Planeta[i].Mueva_V(dt,uno_m2LAMBDA_S2);
		for(i=0; i<N; i++) Planeta[i].Mueva_r(dt,XI);
		GravedadDeNewton.CalculeTodasLasFuerzas(Planeta); for(i=0; i<N; i++) Planeta[i].Mueva_V(dt,LAMBDA);
		for(i=0; i<N; i++) Planeta[i].Mueva_r(dt,uno_m2_XIplusCHI);
		GravedadDeNewton.CalculeTodasLasFuerzas(Planeta); for(i=0; i<N; i++) Planeta[i].Mueva_V(dt,LAMBDA);
		for(i=0; i<N; i++) Planeta[i].Mueva_r(dt,XI);
		GravedadDeNewton.CalculeTodasLasFuerzas(Planeta); for(i=0; i<N; i++) Planeta[i].Mueva_V(dt,uno_m2LAMBDA_S2);
		for(i=0; i<N; i++) Planeta[i].Mueva_r(dt,CHI);
	}
  return 0;
}
