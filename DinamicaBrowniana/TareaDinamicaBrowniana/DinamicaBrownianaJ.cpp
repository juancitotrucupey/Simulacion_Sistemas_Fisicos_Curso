// Metodo de la lanzadera
#include <iostream>
#include <cmath>
#include "Random64.h"
using namespace std;

// const double Deltat=0.01 //pico segundos
const double L=100; //en Amstrongs
const double T=298; //en Kelvin
const double Masa= 22.8916; //en u.m.a
const double e=1; //carga del electrón
const double D=0.132; //en A²/ps
const double kb=0.826; //en esas unidades
const double Gamma = kb*T/(Masa*D); 
const int N=10000; //Cantidad de particulas


//_____________________________________Clases________________________________________________________________
class Cuerpo{
private:
  double x,Vx,Fx,Fxpunto,m,q,R;
public:
  void Inicie(double x0,double Vx0,double Fx0,double m0,double q0,double R0);
  void CalculeFuerza(double E, double dt);
  void Muevase(double dt, Crandom & ran);
  void Dibujese(void);
  double Getx(void){return x;}; //Funcion inline (macro)
  double Getj(void){return q*Vx;};
  };
void Cuerpo::Dibujese(void){
  cout<<", "<<x<<"+"<<R<<"*cos(t),0+"<<R<<"*sin(t)";
}
void Cuerpo::Inicie(double x0,double Vx0,double Fx0,double m0,double q0,double R0){
  x=x0;   Vx=Vx0; Fx=Fx0; m=m0; q=q0;  R=R0;
}
void Cuerpo::CalculeFuerza(double E, double dt){
 double Fxnew=q*E; Fxpunto=(Fxnew-Fx)/dt; Fx=Fxnew;
}
void Cuerpo::Muevase(double dt, Crandom & ran){
double sigma=sqrt(2*D*dt);
 double xnew =  x +(Fx*dt+Fxpunto*dt*dt/2)/(m*Gamma)+ ran.gauss(0, sigma); //Calcular la nueva posición
 Vx=(xnew-x)/dt; //Calcular la velocidad
 x=xnew;  if (x>L) x-=L; if(x<0) x+=L; // condiciones de frontera periodicas
  
}
//-------------------------------------
void InicieAnimacion(void){
//  cout<<"set terminal gif animate"<<endl; 
 // cout<<"set output 'MiBalon.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange [-10:110]"<<endl;
  cout<<"set yrange [-10:10]"<<endl;
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
  Cuerpo Ion[N];
  Crandom ran(29);
  double t,Deltat=0.1,tmax=100;
  double E;
  double xprom, sigma2,jprom;
  int i;


//  InicieAnimacion();

//  void Cuerpo::Inicie(double x0,double Vx0,double Fx0,double m0,double q0,double R0)
for(i=0;i<N;i++){
  
  //Varias Partículas
  // Ion[i].Inicie(i*L/N ,0 ,0 ,Masa ,1,2);
  // Una Partícula
  Ion[i].Inicie(L/2 ,0 ,0 ,Masa ,1,2);
 }
 for(E=-1000;E<1000;E+=1){
  for(t=0;t<tmax;t+=Deltat){
      //  InicieCuadro();
   // Ion.Dibujese();
   // TermineCuadro();
   //EVOLUCIONE
 for(i=0;i<N;i++)   Ion[i].CalculeFuerza(E,Deltat);
  for(i=0;i<N;i++)  Ion[i].Muevase(Deltat,ran);

// Calcular xprom
/*
  for(xprom=0,i=0;i<N;i++){
xprom+=Ion[i].Getx();
xprom /= N;
}
*/
  //  cout<<t<<" "<<xprom<<endl;
//Calcular jprom

  for(jprom=0,i=0;i<N;i++){
jprom+=Ion[i].Getj();
}
jprom /= N;  
//Calcular Sigma2
/*
for(sigma2=0,i=0;i<N;i++)
{sigma2+= (Ion[i].Getx()-xprom)*(Ion[i].Getx()-xprom);
sigma2/=(N-1);
}
*/
//cout<<t<<" "<<sigma2<<endl;
// cout<<t<<" "<<Ion[0].Getx()<<endl;
  }
  
  cout<<E<<" "<<jprom<<endl;
 }
  return 0;
}

/* Tarea;

Set E=0, the print x vs t;
Hacer muchas partículas N= 10000
Como el histograma crece con el numero de partículas
Mostrar como Sigma2 crece linealmente con el tiempo
Mostrar cómo varia la conductividad con el campo electrico (ley de ohm)
*/

 
