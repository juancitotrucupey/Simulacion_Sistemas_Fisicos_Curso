// Osciladores Acoplados (#=N-2)
#include <iostream>
#include <cmath>
using namespace std;

//------ConstantesDelIntegrador-------------------------------
const double CHI   =0.1786178958448091;
const double LAMBDA=-0.2123418310626054;
const double XI    =-0.06626458266981849;
const double uno_m2LAMBDA_S2=(1-2*LAMBDA)/2;
const double uno_m2_XIplusCHI=(1-2*(XI+CHI));
//--------------ConstantesDelSistema------------------------------
const double k=1;
const int N=5;

class Cuerpo;
class Colisionador;

//-------------------Clase cuerpo------------------------------
class Cuerpo{
private:
  double Xinit,M,V,R,F,DeltaX;
public:
  void Inicie(double M0,double X0,double V0,double R0,double DeltaX0);
  void Mueva_r(double dt,double CTE);
  void Mueva_V(double dt,double CTE);
  void AdicioneFuerza(double F0);
  void BorreFuerza(void){F=0;};
  void Dibujese(void);
  double GetDeltaX(void){return DeltaX;}//Desplazamiento; se va a uilizar para calcular la fuerza
  double Getx(void){return Xinit + DeltaX;}; //Posición de reposos más el desplazamiento
  double Gety(void){return 0;}//No hay desplazamiento en la dirección de y
friend class Colisionador;
};

void Cuerpo::Inicie(double M0,double X0,double V0,double R0,double DeltaX0){
  Xinit=X0;  M=M0; V=V0; R=R0;  DeltaX=DeltaX0;
}
void Cuerpo::Mueva_r(double dt,double CTE){
  DeltaX+=V*(CTE*dt);
}
void Cuerpo::Mueva_V(double dt,double CTE){
  V+=F*(CTE*dt/M);
}
void Cuerpo::AdicioneFuerza(double F0){
F+=F0;
}
void Cuerpo::Dibujese(void){
  cout<<", "<<Getx()<<"+"<<R<<"*cos(t),"<<Gety()<<"+"<<R<<"*sin(t)";
}


//-------------Clase colisionador--------------
class Colisionador{
private:
public:  
    void CalculeTodasLasFuerzas(Cuerpo*Oscilador);
    void CalculeFuerzaEntre(Cuerpo & Oscilador1, Cuerpo & Oscilador2);  
};

void Colisionador::CalculeTodasLasFuerzas(Cuerpo*Oscilador){
	int i;
        //Se reinicia la fuerza para cada oscilador
        for(i=0;i<N;i++) Oscilador[i].BorreFuerza();
	//Para cada oscilador quitando el ultimo se calcula la fuerza realizada
	for(i=0;i<N-1;i++)
			CalculeFuerzaEntre(Oscilador[i],Oscilador[i+1]);
}
	
	
void Colisionador::CalculeFuerzaEntre(Cuerpo & Oscilador1, Cuerpo & Oscilador2){ 
        //Calcula la fuerza elastica 
	double dif,x1,x2,Fint;
	x1=Oscilador1.GetDeltaX();   x2=Oscilador2.GetDeltaX();
	dif=x2-x1;
	Fint=k*dif;
		Oscilador2.AdicioneFuerza(-Fint);    Oscilador1.AdicioneFuerza(Fint);
}



//-------------------------------------
void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'MiBalon.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange [-3:43]"<<endl;
  cout<<"set yrange [-4:3]"<<endl;
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
  Cuerpo Oscilador[N];   int i;
  Colisionador Elasticidad;
  double t,Deltat=0.005;
  double m0=1,R=1; 
  double T,tmax,tdibujo;
  int ndibujos=10000;
  
  T=20;
  
  tmax=T;



 InicieAnimacion();

 //Inicie(            double M0,double X0,double V0,double R0,double DeltaX0)
  Oscilador[0].Inicie(m0,           0,         0,         R,           0);
  Oscilador[1].Inicie(m0,          10,         0,         R,           3);
  Oscilador[2].Inicie(m0,          20,         0,         R,           -2*sqrt(2));
  Oscilador[3].Inicie(m0,          30,         0,         R,           0);
  Oscilador[4].Inicie(m0,          40,         0,         R,           0);
  
  for(t=0,tdibujo=0;t<tmax;t+=Deltat,tdibujo+=Deltat){
//cout<<t<<" "<<Oscilador[1].GetDeltaX()<<endl;    
//Animacion
  
    if(tdibujo>tmax/ndibujos){
      InicieCuadro();
      for(i=0;i<N;i++)  Oscilador[i].Dibujese();
      TermineCuadro();
      tdibujo=0;
    }
  
    for(i=1;i<N-1;i++)  Oscilador[i].Mueva_r(Deltat,CHI);
    Elasticidad.CalculeTodasLasFuerzas(Oscilador); for(i=1;i<N-1;i++) Oscilador[i].Mueva_V(Deltat,uno_m2LAMBDA_S2);
    for(i=1;i<N-1;i++) Oscilador[i].Mueva_r(Deltat,XI);
    Elasticidad.CalculeTodasLasFuerzas(Oscilador); for(i=1;i<N-1;i++) Oscilador[i].Mueva_V(Deltat,LAMBDA); 
    for(i=1;i<N-1;i++) Oscilador[i].Mueva_r(Deltat,uno_m2_XIplusCHI);
    Elasticidad.CalculeTodasLasFuerzas(Oscilador); for(i=1;i<N-1;i++) Oscilador[i].Mueva_V(Deltat,LAMBDA); 
    for(i=1;i<N-1;i++) Oscilador[i].Mueva_r(Deltat,XI);
    Elasticidad.CalculeTodasLasFuerzas(Oscilador); for(i=1;i<N-1;i++) Oscilador[i].Mueva_V(Deltat,uno_m2LAMBDA_S2);
    for(i=1;i<N-1;i++) Oscilador[i].Mueva_r(Deltat,CHI);

  }

  return 0;
}
