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



const double g=980;
const double k=1e9;
const double Lx=100;
const double Ly=100;
const int N=1; //número de Granoa
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
  cout<<"set xrange [-10:110]"<<endl;
  cout<<"set yrange [-10:110]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
      cout<<" , "<<100/7<<"*t,0";
     cout<<" , "<<100/7<<"*t,100";
      cout<<" ,0,"<<100/7<<"*t";
       cout<<" ,100,"<<100/7<<"*t";
}
void TermineCuadro(void){
    cout<<endl;
  
}
//.-.--.-.-.-.-.-.-.-.-.CLASE COLISIONADOR.-.-.-.-.-.-.-.-.-.-.--.-.-.-.
class Colisionador{
	private:
	public:
		void CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2);
		void CalculeTodasLasFuerzas (Cuerpo * Grano);
};

void Colisionador::CalculeTodasLasFuerzas (Cuerpo * Grano){
	int i, j;
	//Borrar las fuerzas de todos los planetas
	for(i=0; i<N+4; i++)
		Grano[i].BorreFuerza();
		//Hacer todas las parejas entre planetas y calcular la fuerza en cada pareja
		for(i=0; i<N; i++)
			for(j=i+1;j<N+4;j++)
				CalculeFuerzaEntre(Grano[i],Grano[j]);		
}

void Colisionador::CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2){
  vector3D r21,n,F2;
  double d21,h,Fn;
  r21=Grano2.r-Grano1.r; d21=norma(r21);h=(Grano2.R+Grano1.R)-d21;
  if(h>0)
    Fn=k*pow(h,1.5);
  n=r21/d21;
  F2=n*Fn;

	Grano2.AdicioneFuerza(F2); Grano1.AdicioneFuerza(F2*(-1));
}

int main(void){
  Cuerpo Grano[N+4];int i;
  Colisionador GravedadDeNewton;
  double t,dt=0.1;
  double tdibujo,tmax, ndibujos=1000;
  double R0=1, m0=1;
  double T,VEL0=10,Alpha=M_PI/4,dx,dy;

  T =Lx/VEL0;
  tmax=5*T;
 

 
  
  InicieAnimacion();
  //--------------x0,y0,z0,Vx0,Vy0,Vz0,m0,R0

  //PAREDES

  //Iz
  Grano[N].Inicie(-10000,Ly/2,0,0,0,0,1,10000);//Sol

  //der
   Grano[N+1].Inicie(Lx+10000,Ly/2,0,0,0,0,1,10000);//Sol

   // abajo
   Grano[N+2].Inicie(Lx,-10000,0,0,0,0,1,10000);//Sol
   // arriba
   Grano[N+3].Inicie(Lx,Ly/2+10000,0,0,0,0,1,10000);//Sol
   //__________________________________________________
   Grano[0].Inicie(Lx/3,Ly/2,0,VEL0*cos(Alpha),VEL0*sin(Alpha),0,m0,R0);//Sol
 
// Planeta[1].Inicie(x1,0,0,0,V1,0,m1,1);//Tierra

  
  //------------(x0,y0,Vx0,Vy0,m0   ,R0);
	for(t=0,tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
		//cout<<Planeta[i].Getx()<<" "<<Planeta[i].Gety()<<endl;
		//Animacion

	  
		if(tdibujo>tmax/ndibujos){
			InicieCuadro();
			for(i=0; i<N+4; i++) Grano[i].Dibujese();
			TermineCuadro();
			tdibujo=0;
			}

		for(i=0; i<N; i++) Grano[i].Mueva_r(dt,CHI);
		GravedadDeNewton.CalculeTodasLasFuerzas(Grano); for(i=0; i<N; i++) Grano[i].Mueva_V(dt,uno_m2LAMBDA_S2);
		for(i=0; i<N; i++) Grano[i].Mueva_r(dt,XI);
		GravedadDeNewton.CalculeTodasLasFuerzas(Grano); for(i=0; i<N; i++) Grano[i].Mueva_V(dt,LAMBDA);
		for(i=0; i<N; i++) Grano[i].Mueva_r(dt,uno_m2_XIplusCHI);
		GravedadDeNewton.CalculeTodasLasFuerzas(Grano); for(i=0; i<N; i++) Grano[i].Mueva_V(dt,LAMBDA);
		for(i=0; i<N; i++) Grano[i].Mueva_r(dt,XI);
		GravedadDeNewton.CalculeTodasLasFuerzas(Grano); for(i=0; i<N; i++) Grano[i].Mueva_V(dt,uno_m2LAMBDA_S2);
		for(i=0; i<N; i++) Grano[i].Mueva_r(dt,CHI);
	}
  return 0;
}
