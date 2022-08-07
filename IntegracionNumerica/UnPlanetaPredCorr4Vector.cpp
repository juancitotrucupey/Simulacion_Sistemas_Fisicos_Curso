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
   double getz(void){return r.z();}; //Inline
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
  // float xa,ya,Vxa,Vya,Fxa,Fya;
    vector3D ra,Va,Fa;
    float umass = 1.0/m;
    CalculeFuerza();
    ra=r; Va=V; Fa=F;
    //Predictor
    r +=(0.5*umass*dt*dt*F - V*dt);
    // x = x - Vx*dt + 0.5*umass*dt*dt*Fx;
    // y=y-Vy*dt+0.5*umass*dt*dt*Fy;
    V+= (dt*umass*F);
    // Vx = Vx + dt*umass*Fx;
    // Vy=Vy+dt*umass*Fy;
    CalculeFuerza();
    Fold1=F;
    // Fxold1=Fx; Fyold1=Fy;

    F=Fa;
    // Fx=Fxa; Fy=Fya;
    //Corrector
    r+= ( 2*dt*dt*F - V*2*dt);
    //    x = x - Vx*2*dt + 2*dt*dt*Fx;
    //   y=y-Vy*2*dt+2*umass*dt*dt*Fy;

    V+=(2*umass*dt*F);
    //   Vx = Vx + 2*umass*dt*Fx;
    //  Vy=Vy+2*dt*umass*Fy;
    CalculeFuerza();
    Fold2=F;
    //    Fxold2=Fx; Fyold2=Fy;
    r=ra; V=Va;F=Fa;
    //    x=xa; y=ya; Vx=Vxa; Vy=Vya; Fx=Fxa; Fy=Fya;

}
void Cuerpo::Muevase(double dt){
  vector3D rp,Vp,ra,Va;
  //    float xp,yp,Vxp,Vyp;
  //   float xa,ya,Vxa,Vya;
    float umass = 1./m;
    CalculeFuerza();
    //Predictor

   rp = r + (V*dt + umass*dt*dt*(19*F-10*Fold1+3*Fold2)/24);
   //  xp = x + Vx*dt + umass*dt*dt*(19*Fx - 10*Fxold1 + 3*Fxold2)/24;
   //  yp = y + dt*Vy + umass*dt*dt*(19*Fy - 10*Fyold1 + 3*Fyold2)/24;
   Vp= (rp-r)/dt + umass*dt*(27*F-22*Fold1+7*Fold2)/24;

   //   Vxp = (xp - x)/dt + umass*dt*(27*Fx - 22*Fxold1 + 7*Fxold2)/24;
   //   Vyp = (yp - y)/dt + umass*dt*(27*Fy - 22*Fyold1 + 7*Fyold2)/24;
    //Calculo de nueva Fuerza
   ra =r;

   // xa=x; ya=y; 
   Va=V;
   // Vxa=Vx; Vya=Vy;
   r=rp;
   //  x=xp; y=yp; 
   V =Vp;
   //   Vx=Vxp; Vy=Vyp;
   Fold2 = Fold1;
   //  Fxold2=Fxold1; 
   //  Fyold2=Fyold1; 
   Fold1 =F;
   //   Fxold1=Fx;
   //   Fyold1=Fy;
    CalculeFuerza();
    //Corrector
    r= ra + dt*Va + umass*dt*dt*(3*F+10*Fold1-1*Fold2)/24;
    //    x = xa + dt*Vxa + umass*dt*dt*(3*Fx + 10*Fxold1 - 1*Fxold2)/24;
    //   y = ya + dt*Vya + umass*dt*dt*(3*Fy + 10*Fyold1 - 1*Fyold2)/24;

    V = (r-ra)/dt + umass*dt*(7*F + 6*Fold1 -1*Fold2)/24;
    //    Vx = (x - xa)/dt + umass*dt*(7*Fx + 6*Fxold1 - 1*Fxold2)/24;
    //   Vy = (y - ya)/dt + umass*dt*(7*Fy + 6*Fyold1 - 1*Fyold2)/24;
}

void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
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
  double t,tdibujo,tmax;
  double r=10, m=1.0;
  double Omega,V,T;
  int ndibujos=1000;

  

  Omega=sqrt(GM*pow(r,-3)); V=Omega*r; T=2*M_PI/Omega;
  tmax=1.1*T; 

  // InicieAnimacion();
  
  //------------(m0,R0, x0, y0,z0,Vx0,Vy0,  Vz0);
  Planeta.Inicie(m,1.0,   r, 0,0, 0  ,0.0001*V,  0);
  Planeta.CalculeFuerza();
  Planeta.Arranque(Deltat);
  //   cout<<Planeta.getx()<<" "<<Planeta.gety()<<endl;

  for(t=0,tdibujo=0;t<tmax;t+=Deltat,tdibujo+=Deltat){
     cout<<Planeta.getx()<<" "<<Planeta.gety()<<endl;
    //Animacion
    /*
    if(tdibujo>tmax/ndibujos){
 

      InicieCuadro();
      Planeta.Dibujese();
      TermineCuadro();
      
      // cout<<Planeta.getx()<<" "<<Planeta.gety()<<endl;
      tdibujo=0;
    }
    */
    Planeta.CalculeFuerza();
    Planeta.Muevase(Deltat);
  } 
  
  return 0;
}
