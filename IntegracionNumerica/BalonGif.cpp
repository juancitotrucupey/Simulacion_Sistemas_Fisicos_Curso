// Integracion Numerica
#include <iostream>
#include <cmath>

using namespace std;

const double g=9.795;
  
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
  Fx=0; Fy=-m*g;
}

void Cuerpo::Muevase (double dt){
  x+=Vx*dt;     y+=Vy*dt;
  Vx+=Fx*dt/m;  Vy+=Fy*dt/m; 
}

void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl;
  cout<<"set output 'MiBalon.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange [0:300]"<<endl;
  cout<<"set yrange [-100:100]"<<endl;
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

  Cuerpo Balon;
  double t, dt=0.01;
 

  InicieAnimacion();

 // ----- (  x0, y0, Vx0, Vy0, m0,  R0);
  Balon.Inicie(0,0,30,40,0.453,10);

  for(t=0;t<10;t+=dt)
      {
	//	cout<<Balon.Getx()<<" "<<Balon.Gety()<<endl;
	InicieCuadro();
	Balon.Dibujese();
	TermineCuadro();
	Balon.CalculeFuerza();
	Balon.Muevase(dt);
      }
  
  return 0;
}


