// Integracion Numerica
#include <iostream>
#include <cmath>

using namespace std;

const double GM=1;

  
class Cuerpo{
private:
  double x,y,Vx,Vy,Fx,Fy,m,R;
public:
 void  Inicie( double x0, double y0, double Vx0, double Vy0, double m0, double R0);
  void CalculeFuerza(void); 
  void Muevase (double dt);
  void Dibujese(void);
  double Getx(void){    return x;  }; //Funcion inline (macro) 
 double Gety(void){    return y;  }; //Funcion inline (macro) ; 

};

void Cuerpo::Dibujese(void){
  cout<<", "<<x<<"+"<<R<<"*cos(t),"<<y<<"+"<<R<<"*sin(t)";
}

void Cuerpo::Inicie( double x0, double y0, double Vx0, double Vy0, double m0, double R0){
  x=x0;  y=y0;  Vx=Vx0;  Vy=Vy0;  m=m0;  R=R0;
}

void Cuerpo::CalculeFuerza(void){
  double aux=-GM*m*pow(x*x + y*y,-1.5);
  
  Fx=x*aux; Fy=y*aux;
}

void Cuerpo::Muevase (double dt){
  x+=Vx*dt;     y+=Vy*dt;
  Vx+=Fx*dt/m;  Vy+=Fy*dt/m; 
}

void InicieAnimacion(void){
  /*
    cout<<"set terminal gif animate"<<endl;
    cout<<"set output 'MiPlaneta.gif'"<<endl;
  */  
cout<<"unset key"<<endl;
  cout<<"set xrange [-12:12]"<<endl;
  cout<<"set yrange [-12:12]"<<endl;
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

//double Cuerpo::Getx(void){
// return x;
//}

int main(void){

  Cuerpo Planeta;
  double t, dt=0.0001;
  double r=10; double Omega,V0y;
  double T,tmax,tdibujo;
  int ndibujos=1000;
  Omega = sqrt(GM*pow(r,-3)); T=2*M_PI/Omega;
  V0y=Omega*r;
  tmax= 1.1*T;

  InicieAnimacion();

 // ----- (  x0, y0, Vx0, Vy0, m0,  R0);
  Planeta.Inicie(r,0,0,0.5*V0y,1,1);

for(t=0,tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
    //    cout<<Planeta.Getx()<<" "<<Planeta.Gety()<<endl;
    //Animacion
    if(tdibujo>tmax/ndibujos){
      InicieCuadro();
      Planeta.Dibujese();
      TermineCuadro();
      tdibujo=0;
    }
    Planeta.CalculeFuerza();
    Planeta.Muevase(dt);
    /*
  for(t=0;t<10;t+=dt)
      {
	//	cout<<Balon.Getx()<<" "<<Balon.Gety()<<endl;
	InicieCuadro();
	Balon.Dibujese();
	TermineCuadro();
	Balon.CalculeFuerza();
	Balon.Muevase(dt);
    */      
}
  
  return 0;
}


