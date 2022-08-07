// UnPlaneta por LeapFrog
#include <iostream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"
using namespace std;

const double CHI   =0.1786178958448091;
const double LAMBDA=-0.2123418310626054;
const double XI    =-0.06626458266981849;
const double uno_m2LAMBDA_S2=(1-2*LAMBDA)/2;
const double uno_m2_XIplusCHI=(1-2*(XI+CHI));

const double g=980;
const double K=1e4;
const double Lx=100,Ly=100;
const int Nx=5,Ny=5; const int N=Nx*Ny;

class Cuerpo;
class Colisionador;

//------------------ Clase Cuerpo -----------------------
class Cuerpo{
private:
  vector3D r,V,F;
  double m,R;
public:
  void Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,double m0,double R0);
  void BorreFuerza(void);
  void AdicioneFuerza(vector3D F0);
  void Mueva_r(double dt,double CTE);
  void Mueva_V(double dt,double CTE);
  void Dibujese(void);
  double Getx(void){return r.x();}; //Funcion inline (macro)
  double Gety(void){return r.y();}; //Funcion inline (macro)
  double Getz(void){return r.z();}; //Funcion inline (macro)
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,double m0,double R0){
  r.cargue(x0,y0,z0);   V.cargue(Vx0,Vy0,Vz0);    m=m0;  R=R0;
}
void Cuerpo::BorreFuerza(void){
  F.cargue(0,0,0);
}
void Cuerpo::AdicioneFuerza(vector3D F0){
  F+=F0;
}
void Cuerpo::Mueva_r(double dt,double CTE){
  r+=V*(CTE*dt);
}
void Cuerpo::Mueva_V(double dt,double CTE){
  V+=F*(CTE*dt/m);
}
void Cuerpo::Dibujese(void){
  cout<<", "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}

//------------------ Clase Colisionador -----------------------
class Colisionador{
private:
public:
  void CalculeTodasLasFuerzas(Cuerpo * Grano);
  void CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2);
};
void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Grano){
  int i,j;
  //Borrar las fuerzas de todos los planetas
  for(i=0;i<N+4;i++)    
    Grano[i].BorreFuerza();
  //Hacer todas las parejas entre planetas y calcular la fuerza en cada pareja
  for(i=0;i<N;i++)    
    for(j=i+1;j<N+4;j++)
       CalculeFuerzaEntre(Grano[i],Grano[j]);
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2){
  vector3D r21,n,F2; 
  double d21,h,Fn;
  r21=Grano2.r-Grano1.r; d21=norma(r21); h=(Grano2.R+Grano1.R)-d21;
  if(h>0){
    Fn=K*pow(h,1.5); n=r21/d21; F2=n*Fn;
    Grano2.AdicioneFuerza(F2);  Grano1.AdicioneFuerza(F2*(-1));
  }
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
    cout<<" , 0,"<<100/7<<"*t";
    cout<<" , 100,"<<100/7<<"*t";
}
void TermineCuadro(void){
    cout<<endl;
}


int main(void){
  Cuerpo Grano[N+4]; int i,ix,iy;
  Colisionador Hertz;
  Crandom ran2(0);
  double t,Deltat=1e-3;
  double tmax, tdibujo; int ndibujos=1000;
  double m0=1,R0=3; 
  double T,VEL0=10,Alpha,dx,dy;

  T=Lx/VEL0;  tmax=5*T;
  
  InicieAnimacion();

  //-------------  (x0      ,y0     ,z0,Vx0,Vy0,Vz0,m0,R0   );

  //PAREDES
  //Pared izquierda
    Grano[N].Inicie(-10000  ,Ly/2   , 0,  0,  0,  0, 1,10000); 
  //Pared derecha
  Grano[N+1].Inicie(Lx+10000,Ly/2   , 0,  0,  0,  0, 1,10000); 
  //Pared abajo
  Grano[N+2].Inicie(Lx/2    ,-10000 , 0,  0,  0,  0, 1,10000); 
  //Pared arriba
  Grano[N+3].Inicie(Lx/2    ,Ly+10000,0,  0,  0,  0, 1,10000); 

  //-------------(x0  ,y0  ,z0,Vx0            ,Vy0            ,Vz0,m0,R0);
  //GRANOS
  dx=Lx/(Nx+1);  dy=Ly/(Ny+1);
  for(ix=0;ix<Nx;ix++)
    for(iy=0;iy<Ny;iy++){
      Alpha=2*M_PI*ran2.r();
      Grano[iy*Nx+ix].Inicie(dx*(ix+1),dy*(iy+1), 0,VEL0*cos(Alpha),VEL0*sin(Alpha),  0,m0,R0);
    }
  for(t=0,tdibujo=0;t<tmax;t+=Deltat,tdibujo+=Deltat){
    //cout<<Grano[1].Getx()<<" "<<Grano[1].Gety()<<endl;
    //Animacion
    
    if(tdibujo>tmax/ndibujos){
      InicieCuadro();
      for(i=0;i<N;i++) Grano[i].Dibujese();
      TermineCuadro();
      tdibujo=0;
    }
    
    for(i=0;i<N;i++) Grano[i].Mueva_r(Deltat,CHI);
    Hertz.CalculeTodasLasFuerzas(Grano); for(i=0;i<N;i++) Grano[i].Mueva_V(Deltat,uno_m2LAMBDA_S2);
    for(i=0;i<N;i++) Grano[i].Mueva_r(Deltat,XI);
    Hertz.CalculeTodasLasFuerzas(Grano); for(i=0;i<N;i++) Grano[i].Mueva_V(Deltat,LAMBDA); 
    for(i=0;i<N;i++) Grano[i].Mueva_r(Deltat,uno_m2_XIplusCHI);
    Hertz.CalculeTodasLasFuerzas(Grano); for(i=0;i<N;i++) Grano[i].Mueva_V(Deltat,LAMBDA); 
    for(i=0;i<N;i++) Grano[i].Mueva_r(Deltat,XI);
    Hertz.CalculeTodasLasFuerzas(Grano); for(i=0;i<N;i++) Grano[i].Mueva_V(Deltat,uno_m2LAMBDA_S2);
    for(i=0;i<N;i++) Grano[i].Mueva_r(Deltat,CHI);

  }

  return 0;
}
