#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;


const int Lx=200;
const int Ly=200;


const double C=0.5;
const double C2=C*C;
const double UNO_2C2=1-2*C2;
const double Tau=05;

const double A=10;
const double Lambda=10;
const double Omega=2*M_PI*C/Lambda;




class LatticeBoltzman{
	private:
  double w[5];
  int V[2][5];  // V[x=0][i]=i, i vector velocidad 
		double f[Lx][Ly][5], fnew[Lx][Ly][5]; //int n[ix,i], nnew[ix,i]; ix=la posicion de la celda, i=la direccion . Cambiar n por f, y nnew por fnew
	public:
	         LatticeBoltzman(void);//constructor
                void Inicie(void);
  double Rho(int ix,int iy,int t,bool New);
  double Jx(int ix,int iy);
  double Jy(int ix,int iy);
  double  feq(double Rho0, double Jx0, double Jy0, int i);
  void Imprimase(char * NombreArchivo, int t);
  void ImprimaUnaLinea(char * NombreArchivo, int t);
                void Colisione(int t);
		void Adveccione(void);
               
};

LatticeBoltzman::LatticeBoltzman(void){
  //Cargar los vectores velocidad
  V[0][0]=0; V[1][0]=0;
		V[0][1]=1; V[0][2]=0; 	V[0][3]=-1; V[0][4]=0;
		V[1][1]=0; V[1][2]=1;	V[1][3]=0;  V[1][4]=-1;
  // Cargar los pesos				
  w[0]=1/3; w[1]=w[2]=w[3]=w[4]=1/6;
  }
void LatticeBoltzman::Inicie(void){
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      for(int i=0; i<5; i++)
	f[ix][iy][i]=fnew[ix][iy][i]=0;
  }
  }
  }


double LatticeBoltzman::Rho(int ix,int iy,int t,bool New){
  if(ix==Lx/2 && iy==Ly/2) return sin(Omega*t);
  else{
  double suma=0;
  for(int i=0;i<5;i++)
    if(New)
    suma+=fnew[ix][iy][i];
    else
       suma+=f[ix][iy][i];
  return suma;
}
}

double LatticeBoltzman::Jx(int ix,int iy){
  double suma=0;
  for(int i=0;i<5;i++)
    suma+=f[ix][iy][i]*V[0][i];
  return suma;
}


double LatticeBoltzman::Jy(int ix,int iy){
  double suma=0;
  for(int i=0;i<5;i++)
    suma+=f[ix][iy][i]*V[1][i];
  return suma;
}



double LatticeBoltzman::feq(double Rho0, double Jx0, double Jy0, int i){
  if(i!=0)
    return 3*w[i]*(C2*Rho0 + (V[0][i]*Jx0 + V[1][i]*Jy0));
  else
    Rho0*(1-3*C2*(1-w[0]));
}



void LatticeBoltzman::Imprimase(char * NombreArchivo, int t){
  ofstream MiArchivo(NombreArchivo);
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0; iy<Ly;iy++)
      MiArchivo<<ix<<" "<<iy<<" "<<Rho(ix,iy,t,true)<<endl;    
    MiArchivo<<endl;
  }
  MiArchivo.close();
}

void LatticeBoltzman::ImprimaUnaLinea(char * NombreArchivo, int t){
  ofstream MiArchivo(NombreArchivo);
  int ix=Lx/2;
    for(int iy; iy<Ly;iy++)
      MiArchivo<<iy<<" "<<Rho(ix,iy,t,true)<<endl;
     MiArchivo.close();
		}

void LatticeBoltzman::Colisione(int t){
  int ix,iy,i; double Rho0,Jx0,Jy0;

  for ( ix=0; ix<Lx; ix++){ // Para cada celda
    for ( iy=0; iy<Ly; iy++){
      Rho0=Rho(ix,iy,t,false); Jx0=Jx(ix,iy); Jy0=Jy(ix,iy);
       
       for ( i=0; i<5; i++){  //En cada dirección
	 if(ix==Lx/2 && iy==Ly/2)
	   fnew[ix][iy][i]=feq(Rho0,Jx0,Jy0,i);
	 else
	 fnew[ix][iy][i]=f[ix][iy][i]-1.0/Tau*(f[ix][iy][i]-feq(Rho0,Jx0,Jy0,i));
			
		}
     }
  }
	}

void LatticeBoltzman::Adveccione(void){
	for(int ix=0;ix<Lx;ix++){
	  	for(int iy=0;iy<Ly;iy++){
		  	for(int i=0;i<5;i++){
				f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
		}
	}
}
	}
  //______________________________Funciones Globales_________________________________________________
  

int main(void){
	LatticeBoltzman Ondas;
	int t,tmax=400;
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
