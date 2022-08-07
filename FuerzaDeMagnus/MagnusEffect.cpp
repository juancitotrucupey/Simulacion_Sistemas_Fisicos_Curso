// Efecto Magnus
#include <iostream>
#include <cmath>
#include "Vector.h"
using namespace std;

// Se manejan unidades SI, i.e. Kg, m, s
const double g = 9.8;
const double Ro = 1.224;
const double Alpha = 0.001;
const double Mu = 1.83e-5;

// constantes del integrador
const double CHI   =0.1786178958448091;
const double LAMBDA=-0.2123418310626054;
const double XI    =-0.06626458266981849;
const double uno_m2LAMBDA_S2=(1-2*LAMBDA)/2;
const double uno_m2_XIplusCHI=(1-2*(XI+CHI));



//Enunciando las clases
class Cuerpo;

//............................
class Cuerpo{
private:
  vector3D r, V, F, W, T;
  double m,R,Area,I;
public:
  void Inicie(double x0,double y0, double z0, double Vx0,double Vy0, double Vz0, double Wx0, double Wy0, double Wz0, double m0,double R0);
  void CalculeFuerza(void);
    void Mueva_r(double dt,double CTE);
  void Mueva_V(double dt,double CTE);
   void Dibujese(void);
   double Getx(void){return r.x();}; //Funcion inline (macro)
  double Gety(void){return r.y();}; //Funcion inline (macro)
  double Getz(void){return r.z();}; //Funcion inline (macro)
  
  
};
void Cuerpo::Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0, double Wx0, double Wy0, double Wz0, double m0,double R0){
  r.cargue(x0,y0,z0); V.cargue(Vx0,Vy0,Vz0); W.cargue (Wx0,Wy0,Wz0); m=m0;  R=R0; I=2.0*m*R*R/5.0; Area = M_PI*R*R;
}

void Cuerpo::CalculeFuerza(void){
  double normaV= norma(V);  vector3D cruz = (W^V); double Re,Ca;
 //__--------------------------------------------------------------------------------------------------------------------- 
//Calculo del coeficiente arrastre

Re= 2*Ro*R*normaV/Mu; 

int i,j,N=4;

double XData [N];
double YData [N];
XData[0]=pow(10,5)*(pow(2,0.3)-pow(2,0.575)); YData[0]=0.22;
XData[1]=0; YData[1]=0;
XData[2]=pow(10,5)*(pow(2,0.85)-pow(2,0.575)); YData[2]=-0.22;
XData[3]=pow(10,5)*(pow(2,0.575)-pow(2,0.85)); YData[3]=0.22;
//Parametros del polinomio
double NRe;
double Pol=0;
double Mul=1;

//Parametros del ajuste exponencial

double A, Alpha,SRe;


if(Re < 100000*pow(2,0.3)){
Ca=0.5;
}

else if((100000*pow(2,0.3)<= Re)&(Re <=100000*pow(2,0.85))){
NRe=Re-pow(2,0.575)*pow(10,5);


for(i=0;i<4;i++){
Mul=Mul*YData[i];
for(j=0;j<4;j++){
if(i!=j){
Mul=Mul*(NRe-XData[j])/(XData[i]-XData[j]);
}
}
Pol+=Mul;
Mul=1;
}
Ca =Pol+0.28;
Pol=0;
}

else if(Re> 100000*pow(2,0.85)){

A=0.375-0.06;
Alpha=log(A/(0.375-0.372))/( pow(10,5)*(10-pow(2,0.85)));
SRe=Re-100000*pow(2,0.85);
Ca=0.375-A*exp((-1)*Alpha*SRe);

}



//______________________________________________________________________________________________________________________________________
//Calculo de la fuerza total

  F.cargue(0,0,-m*g); // Fuerza de Gravedad

   F+=(-0.5)*(Ca*Ro*Area*normaV)*V; // Fuerza de Arrastre

  //F_{Magnus]=Cs*0.5*Area*R*Ro*cruz, con Cs=7*Ca y se utilizó un multiplo de Ca para dar cuenta del cambio en la interacción balón aire para grandes velocidades (Esta constante fue fijada numericamente para que la simulación generara los valores esperados)


  F+= 0.5*7*Ca*Area*R*Ro*cruz;

// Esta fuerza no da cuenta de la cantidad de momento angular del balon
/*
double normacruz=norma(cruz);
  F+= cruz*0.5*Ro*0.9*Ca*Area*normaV*normaV/normacruz; //Fuerza de Magnus
*/

 T=(-1)*0.001*W;// Torque 1
// T=(-1)*0.0001*normaV*W; // Torque 2
//T.cargue(0,0,0);
    }

void Cuerpo::Mueva_r(double dt,double CTE){
  r+=V*CTE*dt;
}
void Cuerpo::Mueva_V(double dt,double CTE){
  V+=F*(CTE*dt/m);
    W+=T*(CTE*dt/I);
/*
double normaW=norma(W);
cout<<normaW<<endl;
*/
}
void Cuerpo::Dibujese(void){
  cout<<", "<<r.x()<<"+"<<R<<"*cos(t),"<<r.z()<<"+"<<R<<"*sin(t)";
}
//-------------------------------------
void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'MiBalon.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange [0:5]"<<endl;
  cout<<"set yrange [0:4]"<<endl;
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
  Cuerpo Balon;
  double t,dt=0.000001;
  double r=0.11, m0=0.43;
  double tmax;
  double tdibujo, ndibujos=10000;
  double Alpha1 = 13.3*M_PI/180;
  double Beta = 17*M_PI/180;
  double V = 30.4;
  double Balony, Xrotado, Yrotado;
    
  tmax=10;
  
   InicieAnimacion();


  //  void Inicie(double x0,double y0, double z0, double Vx0,                  double Vy0,       double Vz0, double Wx0, double Wy0, double Wz0, double m0,double R0); 
 Balon.Inicie(2,       0,            0,     V*cos(Beta)*sin(Alpha1), V*cos(Beta)*cos(Alpha1), V*sin(Beta),      0,   -2, 30,      m0, r);
 
  
// La animación termina cuando el balon llega a la linea final de la cacha, i.e. Balony=32
Balony=Balon.Gety();   

   for(t=0,tdibujo=0;Balony<32 ;t+=dt,tdibujo+=dt){
//Rotar para que la velocidad inicial sea paralela al eje y

Xrotado=cos(Alpha1)*Balon.Getx()-sin(Alpha1)*Balon.Gety();
Yrotado=sin(Alpha1)*Balon.Getx()+cos(Alpha1)*Balon.Gety();

     
     // Salida de datos
  //  cout<<Balon.Getx()<<" "<<Balon.Getz()<<endl;
    //cout<<Xrotado<<" "<<Yrotado<<endl;
     //Animacion
     
    if(tdibujo>tmax/ndibujos){
      InicieCuadro();
      Balon.Dibujese();
      TermineCuadro();
      tdibujo=0;
    }
    
    Balon.Mueva_r(dt,CHI);
    Balon.CalculeFuerza(); Balon.Mueva_V(dt,uno_m2LAMBDA_S2);
    Balon.Mueva_r(dt,XI);
    Balon.CalculeFuerza(); Balon.Mueva_V(dt,LAMBDA); 
    Balon.Mueva_r(dt,uno_m2_XIplusCHI);
    Balon.CalculeFuerza(); Balon.Mueva_V(dt,LAMBDA); 
    Balon.Mueva_r(dt,XI);
    Balon.CalculeFuerza(); Balon.Mueva_V(dt,uno_m2LAMBDA_S2);
    Balon.Mueva_r(dt,CHI);
Balony=Balon.Gety();
if(Balony>31.9){t=0; tdibujo=0;  Balon.Inicie(2,       0,            0,     V*cos(Beta)*sin(Alpha1), V*cos(Beta)*cos(Alpha1), V*sin(Beta),      0,   -2, 30,      m0, r);}
//Balony=Yrotado;
  }
//Tiempo que tarda en llegar a la cancha
//cout<<t<<endl;
   
  return 0;
}
