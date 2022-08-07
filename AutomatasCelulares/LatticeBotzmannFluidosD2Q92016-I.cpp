#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=128;
const int Ly=64;

const double C=0.5;     // OJO:  C< sqrt(1/2)=0.707
const double C2=C*C;
const double UM2C2=1-2*C2;
const double Tau=0.53; 
const double RHO0=1.0, UX0=0.06, UY0=0;
const double A=10;
const double Lambda=10;
const double Omega=2*M_PI*C/Lambda;

enum TipoCelda{aire,obstaculo,ventilador};

class LatticeBoltzmann{
private:
  double w[9];
  int V[2][9];  //V[x=0,y=1][i]  , i=el vector velocidad
  double f[Lx][Ly][9],fnew[Lx][Ly][9]; //f[ix][iy][i]  ix,iy=la celda  ,  i=el vector velocidad
  TipoCelda Celda[Lx][Ly];
public:
  LatticeBoltzmann(void);
  void Inicie(void);
  void ConstruyaLaGeometria(voi);
  double Rho(int ix,int iy,int t,bool CalculeConNew);
  double Ux(int ix,int iy,int t,bool CalculeConNew);
  double Uy(int ix,int iy,int t,bool CalculeConNew);
  double feq(double Rho0,double Ux0,double Uy0,int i);
  void Colisione(int t);
  void Adveccione(void);
  void Imprimase(char * NombreArchivo,int t);
 };
LatticeBoltzmann::LatticeBoltzmann(void){
  //D2Q5
  //Cargar los pesos
  w[0]=4.0/9; w[1]=w[2]=w[3]=w[4]=1.0/9;  w[5]=w[6]=w[7]=w[8]=1.0/36;
  //Cargar los vectores velicidad
  V[0][0]=0; 
  V[1][0]=0;

  V[0][1]=1;   V[0][2]=0;   V[0][3]=-1;  V[0][4]=0; 
  V[1][1]=0;   V[1][2]=1;   V[1][3]=0;   V[1][4]=-1;

   V[0][5]=1;   V[0][6]=-1;   V[1][7]=-1;  V[0][8]=1;
    V[1][1]=1;   V[1][6]=1;   V[1][7]=-1;  V[1][8]=-1; 
}

void LatticeBoltzman::ConstruyaLaGeometria(void){
  int ixc=Lx/4,iyc=Ly/2, R=Ly/5;
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      if(ix==0)
	Celda=ventilador;
      else if((ix-ixc)*(ix-ixc)+(iy-iyc)*(iy-iyc)<=R*R)
	Celda=obstaculo
	else
	  Celda=aire;
}
void LatticeBoltzmann::Inicie(void){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<9;i++)
	if(Celda[ix][iy]==obstaculo)
	  f[ix][iy][i]=fnew[ix][iy][i]=feq(RHO0,0,0,i);
	else
	   f[ix][iy][i]=fnew[ix][iy][i]=feq(RHO0,UX0,UY0,i);
}
double LatticeBoltzmann::Rho(int ix,int iy,int t,bool CalculeConNew){
  
    double suma=0;
    for(int i=0;i<9;i++)
      if(CalculeConNew)
	suma+=fnew[ix][iy][i];
      else
	suma+=f[ix][iy][i];
    return suma;
  
}
double LatticeBoltzmann::Ux(int ix,int iy,int t,bool CalculeConNew){
  if(Celda[ix][iy]==ventilador)
    return UX0;
  else if(Celda[ix][iy]==obstaculo)
    return 0;
  else{
  double suma=0;
  for(int i=0;i<9;i++)
    if(CalculeConNew)
    suma+=f[ix][iy][i]*V[0][i];
    else
       suma+=f[ix][iy][i]*V[0][i];
  return suma/Rho(ix,iy,t,CalculeConNew);
}
}
double LatticeBoltzmann::Uy(int ix,int iy,int t,bool CalculeConNew){
  if(Celda[ix][iy]==obstaculo)
    return 0;
  else if(Celda[ix][iy]==ventilador)
    return UY0;
  else{
  double suma=0;
  for(int i=0;i<9;i++)
    if(CalculeConNew)
    suma+=f[ix][iy][i]*V[1][i];
    else
       suma+=f[ix][iy][i]*V[1][i];
    return suma/Rho(ix,iy,t,CalculeConNew);
}
}
double LatticeBoltzmann::feq(double Rho0,double Ux0,double Uy0,int i){
  double UdotVi=(Ux0*V[0][i]+Uy0*V[1][i]);
  double U2=Ux0*Ux0+Uy0*Uy0;
  return w[i]*Rho0*(1+3*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
}
void LatticeBoltzmann::Colisione(int t){
  int ix,iy,i; double Rho0,Ux0,Uy0;
  for(ix=0;ix<Lx;ix++) //Para Cada Celda
    for(iy=0;iy<Ly;iy++){
      Rho0=Rho(ix,iy,t,false); Ux0=Ux(ix,iy);  Uy0=Uy(ix,iy);
      for(i=0;i<5;i++) //en cada direccion      
	if(Celda[ix][iy]==ventilador)
	  fnew[ix][iy][i]=feq(Rho0,UY0,UY0,i);
	else if(Celda[ix][iy]==obstaculo)
	  fnew[ix][iy][i]=feq(Rho0,0,0,i);
      else
	  fnew[ix][iy][i]=f[ix][iy][i]-1.0/Tau*(f[ix][iy][i]-feq(Rho0,Ux0,Uy0,i));
    }
}
void LatticeBoltzmann::Adveccione(void){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<5;i++)
	f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
}
void LatticeBoltzmann::Imprimase(char * NombreArchivo,int t){
  ofstream MiArchivo(NombreArchivo);
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++)
      MiArchivo<<ix<<" "<<iy<<" "<<Ux(ix,iy,t,true)<<" "<<Uy(ix,iy,t,true)<<endl;
    MiArchivo<<endl;
  }
  MiArchivo.close();
}

//------------------------ Funciones Globales ---------------


int main(void){
  LatticeBoltzmann Ondas;
  int t,tmax=100;
  char NombreArchivo[]="Ondas.dat";
  char NombreArchivoLinea[]="Linea.dat";

  Ondas.Inicie();
  for(t=0;t<tmax;t++){
    Ondas.Adveccione();
    Ondas.Colisione(t);
  }
  Ondas.Imprimase(NombreArchivo,t);
  Ondas.ImprimaUnaLinea(NombreArchivoLinea,t);

  return 0;
}
