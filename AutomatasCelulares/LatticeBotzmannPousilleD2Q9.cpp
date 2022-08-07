#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=200;
const int Ly=200;

const double Tau=1.5;
const double RHO0=1.0, UX0=0.00, UY0=0.00;
const double g=0.003;


enum TipoCelda{aire,obstaculo,ventilador};

class LatticeBoltzmann{
private:
  double w[9];
  int V[2][9];  //V[x=0,y=1][i]  , i=el vector velocidad
  double f[Lx][Ly][9],fnew[Lx][Ly][9]; //f[ix][iy][i]  ix,iy=la celda  ,  i=el vector velocidad
  TipoCelda Celda[Lx][Ly];
public:
  LatticeBoltzmann(void);
  void ConstruyaLaGeometria(void);
  void Inicie(void);
  double Rho(int ix,int iy,int t,bool CalculeConNew);
  double Ux(int ix,int iy,int t,bool CalculeConNew);
  double Uy(int ix,int iy,int t,bool CalculeConNew);
  double feq(double Rho0,double Ux0,double Uy0,int i);
  double Fi(double Ux0, double Uy0,double Rho0, int i); 
  void Colisione(int t);
  void Adveccione(bool BounceBack);
  void Grafique(char const * NombreArchivo,int t);
  void Visualizacion_2D_Evolucion(int t, char const * Campo );
  int Opuesto(int ii);
};
LatticeBoltzmann::LatticeBoltzmann(void){
  //D2Q9
  //Cargar los pesos
  w[0]=4.0/9; w[1]=w[2]=w[3]=w[4]=1.0/9; w[5]=w[6]=w[7]=w[8]=1.0/36;
  //Cargar los vectores velicidad
  V[0][0]=0; 
  V[1][0]=0;

  V[0][1]=1;   V[0][2]=0;   V[0][3]=-1;  V[0][4]=0; 
  V[1][1]=0;   V[1][2]=1;   V[1][3]=0;   V[1][4]=-1;

  V[0][5]=1;   V[0][6]=-1;  V[0][7]=-1;  V[0][8]=1; 
  V[1][5]=1;   V[1][6]=1;   V[1][7]=-1;  V[1][8]=-1;
}
void LatticeBoltzmann::ConstruyaLaGeometria(void){
   for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      if((ix==0)||(ix==Lx-1)||(iy==0)||(iy==Ly-1)){
	      Celda[ix][iy]=obstaculo;
      }
      else{
	      Celda[ix][iy]=aire;
      }
    }
  }
}
void LatticeBoltzmann::Inicie(void){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<9;i++)
	      if(Celda[ix][iy]==obstaculo)
	        f[ix][iy][i]=fnew[ix][iy][i]=feq(RHO0,0.0,0.0,i);
	      else
          f[ix][iy][i]=fnew[ix][iy][i]=feq(RHO0,0.0,0.0,i);
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
  double Fext[2]; Fext[1]=Rho(ix,iy,t,CalculeConNew)*g; Fext[0]=0;
  if(Celda[ix][iy]==ventilador)
    return UX0;
  else if(Celda[ix][iy]==obstaculo)
    return 0;
  else{
    double suma=0;
    for(int i=0;i<9;i++)
      if(CalculeConNew)
	      suma+=fnew[ix][iy][i]*V[0][i];
      else
	      suma+=f[ix][iy][i]*V[0][i];
    return (suma+0.5*Fext[0])/Rho(ix,iy,t,CalculeConNew);
  }
}
double LatticeBoltzmann::Uy(int ix,int iy,int t,bool CalculeConNew){
  double Fext[2]; Fext[1]=Rho(ix,iy,t,CalculeConNew)*g; Fext[0]=0;
  if(Celda[ix][iy]==ventilador)
    return UY0;
  else if(Celda[ix][iy]==obstaculo)
    return 0;
  else{
    double suma=0;
    for(int i=0;i<9;i++)
      if(CalculeConNew)
	      suma+=fnew[ix][iy][i]*V[1][i];
      else
	      suma+=f[ix][iy][i]*V[1][i];
    return (suma+0.5*Fext[1])/Rho(ix,iy,t,CalculeConNew);
  }
}
double LatticeBoltzmann::feq(double Rho0,double Ux0,double Uy0,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i];
  double U2=Ux0*Ux0+Uy0*Uy0;
  double Normal = w[i]*(3.0*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
  
  if(i==0){//Fluido 
    return Rho0*(0.1-U2*(2.0/3.0)); //v0
  }
  else if((1<=i)&&(i<=4)){
    return Rho0*(((1.0-0.1)/5.0) + Normal);//v1,...,v4
  }
  else{
    return Rho0*(((1.0-0.1)/20.0) + Normal);//v5,...,v8
  }
  
  //return w[i]*Rho0*(1+3*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
}

double LatticeBoltzmann::Fi(double Ux0, double Uy0,double Rho0, int i){
  double Fext[2]; Fext[1]=Rho0*g; Fext[0]=0;

  double VidotFext, VidotU, UdotFext;

  VidotFext=Fext[0]*V[0][i]+Fext[1]*V[1][i];
  VidotU= V[0][i]*Ux0+V[1][i]*Uy0;
  UdotFext=Ux0*Fext[0]+Uy0*Fext[1];

  return (1-1/(2*Tau))*w[i]*(3*(VidotFext-UdotFext)+9*VidotU*VidotFext);
}

void LatticeBoltzmann::Colisione(int t){
  int ix,iy,i; double Rho0,Ux0,Uy0;
  for(ix=0;ix<Lx;ix++) //Para Cada Celda
    for(iy=0;iy<Ly;iy++){
      Rho0=Rho(ix,iy,t,false); Ux0=Ux(ix,iy,t,false);  Uy0=Uy(ix,iy,t,false);
      for(i=0;i<9;i++) //en cada direccion      
	      if(Celda[ix][iy]==ventilador)
	        fnew[ix][iy][i]=feq(Rho0,UX0,UY0,i);
	      else if(Celda[ix][iy]==obstaculo)
	        fnew[ix][iy][i]=feq(Rho0,0,0,i);
	      else
	        fnew[ix][iy][i]=f[ix][iy][i]-1.0/Tau*(f[ix][iy][i]-feq(Rho0,Ux0,Uy0,i))+Fi(Ux0,Uy0,Rho0,i);
    }
}
void LatticeBoltzmann::Adveccione(bool BounceBack){
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      if(BounceBack){
        if(aire==Celda[ix][iy]){
          for(int i=0;i<9;i++){
            if(obstaculo==Celda[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly]){//Si el nodo se encuentra junto a un pared
              f[ix][iy][Opuesto(i)]=fnew[ix][iy][i];              //las funciones de densidad son reflejadas 
            }
            else{
              f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
            }  
          }
        }
      }
      else{
        for(int i=0;i<9;i++){
          f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];    
        }
      }
    }
  }
}

void LatticeBoltzmann::Grafique(char const * NombreArchivo,int t){
  ofstream MiArchivo(NombreArchivo);
  int iy=0;
  for(int ix=0;ix<Lx;ix++){
    MiArchivo<<ix<<" "<<Uy(ix,iy,t,true)<<endl;
  }
  MiArchivo.close();
}

void LatticeBoltzmann::Visualizacion_2D_Evolucion(int t, char const * Campo ){
  //cout<<"set xrange [0:"<<Lx<<"]"<<endl;
  //cout<<"set yrange [0:"<<Ly<<"]"<<endl;
  cout<<"set size ratio 1"<<endl;
  //cout<<"set palette defined (1 \"black\",2 \"red\")"<<endl;
  cout<<"set view map"<<endl;
  cout<<"set title \" "<<t<<"\" "<<endl;
  //cout<<"set countour base"<<endl;
  //cout<<"splot '-' using 2:1:3 with pm3d title \""<<t<<"\""<<endl;
  int ix,iy;
  double Phi;

  //Se decide que variable va a ser graficada a partir del valor de la variable Campo
  
  //La magnitud del campo de velocidades
  if(Campo=="MagnitudVelocidad"){
    //cout<<"set palette defined ("<<0<< " \"blue\")"<<endl;//,"<<max(OmegaB,OmegaR)<<" \"red\")"<<endl;
    cout<<"plot '-' using 2:1:3 with image pixels title \""<<t<<"\""<<endl;
    for(ix=0;ix<Lx;ix++){
      for(iy=0;iy<Ly;iy++){
        //Graficar la magnitud de la velocidad
        cout<<ix+0.5<<" "<<iy+0.5<<" "<<sqrt(Ux(ix,iy,t,false)*Ux(ix,iy,t,false) + Uy(ix,iy,t,false)*Uy(ix,iy,t,false))<<endl;
        //cout<<ix<<" "<<iy<<" "<<Rho(ix,iy,t,false)<<endl;//Graficar indicador
      }
      if(Lx-1==ix){
        cout<<"EOF"<<endl;
      }
      else{
        cout<<endl; 
      }
    } 
  }
  //Tipo de celda 
  if(Campo=="TipoCelda"){
    cout<<"set palette defined ("<<"-1"<< " \"blue\","<<"1"<<" \"red\")"<<endl;
    cout<<"plot '-' using 2:1:3 with image pixels title \""<<t<<"\""<<endl;
    for(ix=0;ix<Lx;ix++){
      for(iy=0;iy<Ly;iy++){
        if(obstaculo==Celda[ix][iy]){
          cout<<ix+0.5<<" "<<iy+0.5<<" "<<"-1"<<endl;
        }
        else if (aire==Celda[ix][iy]){
          cout<<ix+0.5<<" "<<iy+0.5<<" "<<"1"<<endl;
        }
        
      }
      if(Lx-1==ix){
        cout<<"EOF"<<endl;
      }
      else{
        cout<<endl; 
      }
    } 
  }
  //La densidad en cada celda
  if(Campo=="Densidad"){
    //cout<<"set palette defined ("<<0<< " \"blue\","<<RHOB<<" \"red\")"<<endl;
    cout<<"plot '-' using 2:1:3 with image pixels title \""<<t<<"\""<<endl;
    for(ix=0;ix<Lx;ix++){
      for(iy=0;iy<Ly;iy++){
        cout<<ix+0.5<<" "<<iy+0.5<<" "<<Rho(ix,iy,t,false)<<endl;//Graficar indicador
      }
      if(Lx-1==ix){
        cout<<"EOF"<<endl;
      }
      else{
        cout<<endl; 
      }
    } 
  }
}

//Se implementa una función que calcula el vector opuesto de unos de los vectores de transporte
int LatticeBoltzmann::Opuesto(int ii){
    if(ii==2){return 4;}
    if(ii==4){return 2;}
    if(ii==1){return 3;}
    if(ii==3){return 1;}
    if(ii==5){return 7;}
    if(ii==7){return 5;}
    if(ii==6){return 8;}
    if(ii==8){return 6;}
}




//------------------------ Funciones Globales ---------------


int main(void){
  LatticeBoltzmann Aire;
  int t,tmax=20000;

  Aire.ConstruyaLaGeometria();
  Aire.Inicie();
  for(t=0;t<tmax;t++){
    if((t%100)==0){
      Aire.Visualizacion_2D_Evolucion(t,"MagnitudVelocidad");
      //cout<<"pause 0.5"<<endl;
    }
    //Aire.Adveccione();
    Aire.Colisione(t);
    Aire.Adveccione(true);
    
  }
  //Aire.Grafique("Rio_JD.dat",t);

  return 0;
}
