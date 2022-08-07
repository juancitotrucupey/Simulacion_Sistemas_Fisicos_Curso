// UnPlaneta por LeapFrog
#include <iostream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"
using namespace std;
// Paquete de software libre para dinamica molecular Gromacs
const double CHI   =0.1786178958448091;
const double LAMBDA=-0.2123418310626054;
const double XI    =-0.06626458266981849;
const double uno_m2LAMBDA_S2=(1-2*LAMBDA)/2;
const double uno_m2_XIplusCHI=(1-2*(XI+CHI));

const double g=1.0;
const double K=500, Gamma=50, Kcundall=5, MU=0.4;
const double Lx=100,Ly=100;
const int Nx=1,Ny=1; const int N=Nx*Ny;

class Cuerpo;
class Colisionador;

//------------------ Clase Cuerpo -----------------------
class Cuerpo{
private:
  vector3D r,V,F,omega,tau;
  double m,R,I,theta;
public:
  void Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,double omegax0,double omegay0
	      ,double omegaz0,double theta0,double m0,double R0);
  void BorreFuerzaYTorque(void);
  void AdicioneFuerza(vector3D F0);
  void AdicioneTorque(vector3D tau0);
  void AgregueFuerzaGravedad(void);
  void Mueva_r(double dt,double CTE);
  void Mueva_V(double dt,double CTE);
  void Dibujese(void);
  double Getx(void){return r.x();}; //Funcion inline (macro)
  double Gety(void){return r.y();}; //Funcion inline (macro)
  double Getz(void){return r.z();}; //Funcion inline (macro)
  double GetVx(void){return V.x();}; //Funcion inline (macro)
  double GetVy(void){return V.y();}; //Funcion inline (macro)
  double GetVz(void){return V.z();}; //Funcion inline (macro)
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,double omegax0,double omegay0
		    ,double omegaz0,double theta0,double m0,double R0){
  r.cargue(x0,y0,z0);   V.cargue(Vx0,Vy0,Vz0);   omega.cargue(omegax0,omegay0,omegaz0); theta=theta0;
  m=m0;  R=R0;  I=2.0/5*m*R*R; 
}
void Cuerpo::BorreFuerzaYTorque(void){
  F.cargue(0,0,0);  tau.cargue(0,0,0);
}
void Cuerpo::AgregueFuerzaGravedad(void){
  vector3D peso;
  peso.cargue(0,-m*g,0); F+=peso;
}
void Cuerpo::AdicioneFuerza(vector3D F0){
  F+=F0;
}
void Cuerpo::AdicioneTorque(vector3D tau0){
  tau+=tau0;
}
void Cuerpo::Mueva_r(double dt,double CTE){
  r+=V*(CTE*dt);  theta+=omega.z()*(CTE*dt);
}
void Cuerpo::Mueva_V(double dt,double CTE){
  V+=F*(CTE*dt/m);  omega+=tau*(CTE*dt/I);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t) , "
      <<r.x()<<"+"<<R*cos(theta)/7.0<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7.0<<"*t"; 
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
    Grano[i].BorreFuerzaYTorque();
  //Agregar Fuerza de Gravedad
  for(i=0;i<N;i++)    
    Grano[i].AgregueFuerzaGravedad();
  //Hacer todas las parejas entre planetas y calcular la fuerza en cada pareja
  for(i=0;i<N;i++)    
    for(j=i+1;j<N+4;j++)
       CalculeFuerzaEntre(Grano[i],Grano[j]);
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2){
  vector3D r21,F2,Tau2,Fn,Ft,Vc,Vcn,Vct,n,t; 
  double d21,h,normaFn,Ftmax,m1,m2,m12,R1,R2,normaVcn,normaVct,normaFt,ERFF=1e-9;
  r21=Grano2.r-Grano1.r; d21=norma(r21); h=(Grano2.R+Grano1.R)-d21;
  if(h>0){
    //Geometría y dinámica del contacto
    m1=Grano1.m;   m2=Grano2.m;   m12=(m1*m2)/(m1+m2);
    R1=Grano1.R;   R2=Grano2.R;
    n=r21/d21;
    //Calcular velocidad de contacto y el vector tangente
    Vc=(Grano2.V-Grano1.V)-(((Grano2.omega*R2)+(Grano1.omega*R1))^n);
    normaVcn=Vc*n; Vcn=n*normaVcn; Vct=Vc-Vcn;  normaVct=norma(Vct);
    if(normaVct<ERFF) t.cargue(0,0,0); else t=Vct/normaVct;

    //FUERZAS NORMALES
    //Fuerza de Hertz
    normaFn=K*pow(h,1.5); 
    //Disipacion
    normaFn-=m12*sqrt(h)*Gamma*normaVcn; if(normaFn<0) normaFn=0;
    Fn=n*normaFn;
    Ft.cargue(0,0,0);
    
    //Construir la fuerza total
    F2=Fn+Ft;
    Grano2.AdicioneFuerza(F2);      Grano2.AdicioneTorque((n*(-R2))^Ft);      
    Grano1.AdicioneFuerza(F2*(-1)); Grano1.AdicioneTorque((n*R1)^(Ft*(-1)));      
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
  Colisionador Cundall;
  Crandom ran2(0);
  double t,Deltat=1e-5;
  double tmax, tdibujo; int ndibujos=1000;
  double m0=1,R0=10;
  double T,VEL0=0,Alpha,dx,dy;

  T=sqrt(2*Lx/g);  tmax=3*T;
  
    InicieAnimacion();

  //-------------  (x0      ,y0     ,z0,Vx0,Vy0,Vz0,omegax0,omegay0,omegaz0,theta0,m0,R0   );
  //PAREDES
  //Pared izquierda
    Grano[N].Inicie(-10000  ,Ly/2   , 0,  0,  0,  0,      0,       0,      0,    0, 1,10000); 
  //Pared derecha
  Grano[N+1].Inicie(Lx+10000,Ly/2   , 0,  0,  0,  0,      0,       0,      0,    0, 1,10000); 
  //Pared abajo
  Grano[N+2].Inicie(Lx/2,  -10000   , 0,  0,  0,  0,      0,       0,      0,    0, 1,10000); 
  //Pared arriba
  Grano[N+3].Inicie(Lx/2,Ly+10000   , 0,  0,  0,  0,      0,       0,      0,    0, 1,10000); 

  //-------------  (x0      ,y0     ,z0,Vx0,Vy0,Vz0,omegax0,omegay0,omegaz0,theta0,m0,R0   );
 //GRANOS
  dx=Lx/(Nx+1); dy=Ly/(Ny+1); if(dx/3 > dy/3) R0=dy/3; else R0=dx/3;
  for(i=0;i<N;i++){
    Alpha=2*M_PI*ran2.r();
    Grano[i].Inicie(((i%Nx)+1)*dx,((i/Ny)+1)*dy,0,VEL0*cos(Alpha),VEL0*sin(Alpha),0,0,0,1,0,m0,R0);
  }
   

  for(t=0,tdibujo=0;t<tmax;t+=Deltat,tdibujo+=Deltat){
    //    cout<<t<<" "<<Grano[0].GetVx()*m0+Grano[1].GetVx()*m1+Grano[2].GetVx()*m2<<endl;
    //Animacion
    
    if(tdibujo>tmax/ndibujos){
      //    cout<<t<<" "<<Grano[0].Gety()<<endl;
            InicieCuadro();
            for(i=0;i<N;i++) Grano[i].Dibujese();
            TermineCuadro();
      tdibujo=0;
    }
    
    for(i=0;i<N;i++) Grano[i].Mueva_r(Deltat,CHI);
    Cundall.CalculeTodasLasFuerzas(Grano); for(i=0;i<N;i++) Grano[i].Mueva_V(Deltat,uno_m2LAMBDA_S2);
    for(i=0;i<N;i++) Grano[i].Mueva_r(Deltat,XI);
    Cundall.CalculeTodasLasFuerzas(Grano); for(i=0;i<N;i++) Grano[i].Mueva_V(Deltat,LAMBDA); 
    for(i=0;i<N;i++) Grano[i].Mueva_r(Deltat,uno_m2_XIplusCHI);
    Cundall.CalculeTodasLasFuerzas(Grano); for(i=0;i<N;i++) Grano[i].Mueva_V(Deltat,LAMBDA); 
    for(i=0;i<N;i++) Grano[i].Mueva_r(Deltat,XI);
    Cundall.CalculeTodasLasFuerzas(Grano); for(i=0;i<N;i++) Grano[i].Mueva_V(Deltat,uno_m2LAMBDA_S2);
    for(i=0;i<N;i++) Grano[i].Mueva_r(Deltat,CHI);

  }

  
  return 0;
}
