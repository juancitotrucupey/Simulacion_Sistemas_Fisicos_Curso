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

const double Deltat=1e-4;
const double g=1.0;
const double K=5000, Gamma=50, Kcundall=5, MU=0.4;
const double Lx=100,Ly=100;
const int Nx=5,Ny=5; const int N=Nx*Ny;

class Cuerpo;
class Colisionador;

//------------------ Clase Cuerpo -----------------------
class Cuerpo{
private:
  vector3D r,V,F,omega,tau;
  double m,R,I,theta;
public:
  void Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0
	      ,double omegax0,double omegay0,double omegaz0,double theta0,double m0,double R0);
  void BorreFuerzaYTorque(void);
  void AgregueFuerzaGravedad(void);
  void AdicioneFuerza(vector3D F0);
  void AdicioneTorque(vector3D Tau0);
  void Mueva_r(double dt,double CTE);
  void Mueva_V(double dt,double CTE);
  void Dibujese(void);
  double Getx(void){return r.x();}; //Funcion inline (macro)
  double Gety(void){return r.y();}; //Funcion inline (macro)
  double Getz(void){return r.z();}; //Funcion inline (macro)
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0
		    ,double omegax0,double omegay0,double omegaz0,double theta0,double m0,double R0){
  r.cargue(x0,y0,z0);   V.cargue(Vx0,Vy0,Vz0);    omega.cargue(omegax0,omegay0,omegaz0);
  m=m0;  R=R0;   I=2.0/5*m*R*R; theta=theta0;
}
void Cuerpo::BorreFuerzaYTorque(void){
  F.cargue(0,0,0);   tau.cargue(0,0,0); 
}
void Cuerpo::AgregueFuerzaGravedad(void){
  vector3D Fg; Fg.cargue(0,-m*g,0); F+=Fg;
}
void Cuerpo::AdicioneFuerza(vector3D F0){
  F+=F0;
}
void Cuerpo::AdicioneTorque(vector3D Tau0){
  tau+=Tau0;
}

void Cuerpo::Mueva_r(double dt,double CTE){
  r+=V*(CTE*dt); theta+=omega.z()*(CTE*dt);
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
  vector3D dcontacto[N+4][N+4]; double hold[N+4][N+4];
public:
  void Inicie(void);
  void CalculeTodasLasFuerzas(Cuerpo * Grano,double dt);
  void CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2,vector3D &dcontacto,double &hold,double dt);
  double GetEp(void);
};
void Colisionador::Inicie(void){
  for(int i=0;i<N+4;i++)  
    for(int j=0;j<N+4;j++){
      dcontacto[i][j].cargue(0,0,0);
      hold[i][j]=0;
    }
}
double Colisionador::GetEp(void){
  double Ep=0,h;
  for(int i=0;i<N;i++)  
    for(int j=i+1;j<N+4;j++){
      Ep+=Kcundall*0.5*norma2(dcontacto[i][j]);
      if(hold[i][j]>0) Ep+=K*pow(h,3.5)/3.5;
    }
  return Ep;
}
void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Grano,double dt){
  int i,j;
  //Borrar las fuerzas de todos los planetas
  for(i=0;i<N+4;i++)    
    Grano[i].BorreFuerzaYTorque();
  //Agregar Fuerza de Gravedad
  for(i=0;i<N;i++)    
    Grano[i].AgregueFuerzaGravedad();
  Grano[N].AgregueFuerzaGravedad(); //Pared superior m??vil
  //Hacer todas las parejas entre planetas y calcular la fuerza en cada pareja
  for(i=0;i<N;i++)    
    for(j=i+1;j<N+4;j++)
      CalculeFuerzaEntre(Grano[i],Grano[j],dcontacto[i][j],hold[i][j],dt);
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2,vector3D &dcontacto,double &hold,double dt){
  vector3D r21,F2,Tau2,Fn,Ft,Vc,Vcn,Vct,n,t; 
  double d21,h,normaFn,Ftmax,m1,m2,m12,R1,R2,normaVcn,normaVct,normaFt,ERFF=1e-9;
  r21=Grano2.r-Grano1.r; d21=norma(r21); h=(Grano2.R+Grano1.R)-d21;
  if(h>0){
    //Geometr??a y din??mica del contacto
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
    //FUERZAS TANGENCIALES
    //fuerza est??tica
    dcontacto+=(Vct*dt);
    Ft=dcontacto*(-Kcundall);
    //fuerza cin??tica
    Ftmax=MU*normaFn; normaFt=norma(Ft);
    if(normaFt>Ftmax) Ft=t*(-Ftmax);
    
    //Construir la fuerza total
    F2=Fn+Ft;
    Grano2.AdicioneFuerza(F2);      Grano2.AdicioneTorque((n*(-R2))^Ft);      
    Grano1.AdicioneFuerza(F2*(-1)); Grano1.AdicioneTorque((n*R1)^(Ft*(-1)));      
  }
  else if(hold>0){
    //Reiniciar dcontacto en ceros
    dcontacto.cargue(0,0,0);
  }
  hold=h;
}

//-------------------------------------
void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'Gas2D.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange [-10:110]"<<endl;
  cout<<"set yrange [-10:110]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(double hParedSup){
    cout<<"plot 0,0 ";
    cout<<" , "<<100/7<<"*t,0";
    cout<<" , "<<100/7<<"*t,"<<hParedSup;
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
  double t;
  double tmax, tdibujo; int ndibujos=1000;
  double m0,R00,R0,omega0;
  double T,VEL0=1,Alpha,dx,dy,hParedSup;
  double Presion=1000,rho0=1,InertialNumber=1e-1,OMEGOTA;
  bool test=false;

  T=sqrt(2*Lx/g);  tmax=3*T;
  dx=Lx/(Nx+1);  dy=Ly/(Ny+1); if(dx<dy) R00=dx/3; else R00=dy/3;
  OMEGOTA=InertialNumber*sqrt(Presion/rho0)/(10000*R00);

  InicieAnimacion();

  //-------------  (x0      ,y0     ,z0,Vx0,Vy0,Vz0,omegax0,omegay0,omegaz0,theta0,m0,R0   );

  //PAREDES
   //Pared arriba
    Grano[N].Inicie(Lx/2    ,Ly+10000,0,  0,  0,  0,      0,      0,OMEGOTA,    0,Presion/Lx,10000); 
 //Pared izquierda
  Grano[N+1].Inicie(-10000  ,Ly/2   , 0,  0,  0,  0,      0,      0,       0,    0,10,10000); 
  //Pared abajo
  Grano[N+2].Inicie(Lx/2    ,-10000 , 0,  0,  0,  0,      0,      0,       0,    0, 1,10000); 
  //Pared derecha
  Grano[N+3].Inicie(Lx+10000,Ly/2   , 0,  0,  0,  0,      0,      0,       0,    0, 1,10000); 

  //-------------  (x0      ,y0     ,z0,Vx0,Vy0,Vz0,omegax0,omegay0,omegaz0,theta0,m0,R0   );
  //GRANOS
  for(ix=0;ix<Nx;ix++)
    for(iy=0;iy<Ny;iy++){
      Alpha=2*M_PI*ran2.r();
      omega0=2*ran2.r()-1;
      R0=(0.8+0.4*ran2.r())*R00;
      m0=4.0/3*M_PI*rho0*pow(R0,3);
      Grano[iy*Nx+ix].Inicie(dx*(ix+1),dy*(iy+1), 0,VEL0*cos(Alpha),VEL0*sin(Alpha),0,0,0,omega0,0,m0,R0);
    }
  for(t=0,tdibujo=0;t<tmax;t+=Deltat,tdibujo+=Deltat){
    //Animacion
    if(tdibujo>tmax/ndibujos){
      hParedSup=Grano[N].Gety()-10000;
      InicieCuadro(hParedSup);
      for(i=0;i<N;i++) Grano[i].Dibujese();
      TermineCuadro();
      tdibujo=0;
    }
    
    for(i=0;i<=N;i++) Grano[i].Mueva_r(Deltat,CHI);
    Hertz.CalculeTodasLasFuerzas(Grano,Deltat); for(i=0;i<=N;i++) Grano[i].Mueva_V(Deltat,uno_m2LAMBDA_S2);
    for(i=0;i<=N;i++) Grano[i].Mueva_r(Deltat,XI);
    Hertz.CalculeTodasLasFuerzas(Grano,Deltat); for(i=0;i<=N;i++) Grano[i].Mueva_V(Deltat,LAMBDA); 
    for(i=0;i<=N;i++) Grano[i].Mueva_r(Deltat,uno_m2_XIplusCHI);
    Hertz.CalculeTodasLasFuerzas(Grano,Deltat); for(i=0;i<=N;i++) Grano[i].Mueva_V(Deltat,LAMBDA); 
    for(i=0;i<=N;i++) Grano[i].Mueva_r(Deltat,XI);
    Hertz.CalculeTodasLasFuerzas(Grano,Deltat); for(i=0;i<=N;i++) Grano[i].Mueva_V(Deltat,uno_m2LAMBDA_S2);
    for(i=0;i<=N;i++) Grano[i].Mueva_r(Deltat,CHI);

    
  }

  return 0;
}
