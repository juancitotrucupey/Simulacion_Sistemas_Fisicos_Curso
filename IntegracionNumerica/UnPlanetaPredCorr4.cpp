//Un Planeta Predictor Corrector de Orden 4
#include <iostream>
#include <cmath>
#include "Vector.h" 
using namespace std;

const double GM=1.0;
const double Deltat=0.1;

//------------------- Clase Cuerpo ----------------------------------------
class Cuerpo{
private:
  vector3D r,V,F,Fold1,Fold2;
  double m,R;
 
public:
  void Inicie(double m0,double R0,double x0,double y0, double z0,double Vx0,double Vy0, double Vz0);
  void CalculeFuerza(void);
  void Arranque(double dt);
  void Muevase(double dt);
  void Dibujese(void);
  double getx(void){return r.x();}; //Inline
  double gety(void){return r.y();}; //Inline
};
void Cuerpo::Inicie(double m0,double R0,double x0,double y0, double z0 ,double Vx0,double Vy0, double Vz0){
  r.cargue(x0,y0,z0); V.cargue(Vx0,Vy0,Vz0);  m=m0; R=R0;
}
void Cuerpo::CalculeFuerza(void){
  //Fuerza de gravedad
  double aux=GM*m*pow(norma2(r),-1.5);
  F=r*aux;
}
void Cuerpo::Arranque(double dt){
    float xa,ya,Vxa,Vya,Fxa,Fya;
    vector3D ra,Va,Fa;
    float umass = 1.0/m;
    CalculeFuerza();
    ra.cargue(r.x(),rxa=x; ya=y; Vxa=Vx; Vya=Vy; Fxa=Fx; Fya=Fy;
    //Predictor
    x = x - Vx*dt + 0.5*umass*dt*dt*Fx;
    y=y-Vy*dt+0.5*umass*dt*dt*Fy;
    Vx = Vx + dt*umass*Fx;
    Vy=Vy+dt*umass*Fy;
    CalculeFuerza();
    Fxold1=Fx; Fyold1=Fy;

    Fx=Fxa; Fy=Fya;
    //Corrector
    x = x - Vx*2*dt + 2*dt*dt*Fx;
    y=y-Vy*2*dt+2*umass*dt*dt*Fy;
    Vx = Vx + 2*umass*dt*Fx;
    Vy=Vy+2*dt*umass*Fy;
    CalculeFuerza();
    Fxold2=Fx; Fyold2=Fy;

    x=xa; y=ya; Vx=Vxa; Vy=Vya; Fx=Fxa; Fy=Fya;

}
void Cuerpo::Muevase(double dt){
    float xp,yp,Vxp,Vyp;
    float xa,ya,Vxa,Vya;
    float umass = 1./m;
    CalculeFuerza();
    //Predictor
    xp = x + Vx*dt + umass*dt*dt*(19*Fx - 10*Fxold1 + 3*Fxold2)/24;
    yp = y + dt*Vy + umass*dt*dt*(19*Fy - 10*Fyold1 + 3*Fyold2)/24;
    Vxp = (xp - x)/dt + umass*dt*(27*Fx - 22*Fxold1 + 7*Fxold2)/24;
    Vyp = (yp - y)/dt + umass*dt*(27*Fy - 22*Fyold1 + 7*Fyold2)/24;
    //Calculo de nueva Fuerza
    xa=x; ya=y; 
    Vxa=Vx; Vya=Vy;
    x=xp; y=yp; 
    Vx=Vxp; Vy=Vyp;
    Fxold2=Fxold1; 
    Fyold2=Fyold1; 
    Fxold1=Fx;
    Fyold1=Fy;
    CalculeFuerza();
    //Corrector
    x = xa + dt*Vxa + umass*dt*dt*(3*Fx + 10*Fxold1 - 1*Fxold2)/24;
    y = ya + dt*Vya + umass*dt*dt*(3*Fy + 10*Fyold1 - 1*Fyold2)/24;
    Vx = (x - xa)/dt + umass*dt*(7*Fx + 6*Fxold1 - 1*Fxold2)/24;
    Vy = (y - ya)/dt + umass*dt*(7*Fy + 6*Fyold1 - 1*Fyold2)/24;
}

void Cuerpo::Dibujese(void){
    cout<<" , "<<x<<"+"<<R<<"*cos(t),"<<y<<"+"<<R<<"*sin(t)";
}

//------------------- Funciones Globales de Animacion ------------------
void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'pelicula.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-12:12]"<<endl;
  cout<<"set yrange[-12:12]"<<endl;
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

int main(){
  Cuerpo Planeta;
  double t,tdibujar;
  double r=10, m=1.0;
  double omega,V,T;

  omega=sqrt(GM/(r*r*r)); V=omega*r; T=2*M_PI/omega;
  
  //------------(m0,R0, x0,y0,Vx0,Vy0);
  Planeta.Inicie(m,0.5,   r, 0, 0,0.5*V);
  Planeta.CalculeFuerza();
  Planeta.Arranque(Deltat);

  //  InicieAnimacion();
  for(t=tdibujar=0;t<1.1*T;t+=Deltat,tdibujar+=Deltat){
    if(tdibujar>=T/100){
      /*
      InicieCuadro();
      Planeta.Dibujese();
      TermineCuadro();
      */   
      cout<<Planeta.getx()<<" "<<Planeta.gety()<<endl;
      tdibujar=0;
    }
    Planeta.CalculeFuerza();
    Planeta.Muevase(Deltat);
  } 
  
  return 0;
}
