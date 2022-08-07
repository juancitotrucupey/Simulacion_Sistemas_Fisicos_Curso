//Lattice Boltzman Electrodinamica
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=1, Ly=1, Lz=200;

const double Tau=0.5;

const double mur=1, sigma=0;
const double mu0=2, epsilon0=1;
const double mu=mu0*mur;
const int tmax=140;


//----------Funciones Medio----------
double epsilonr(int ix, int iy, int iz){
  return 1.75+0.75*tanh(iz-Lz/2);
}

double epsilon(double epsilonr0){
  return epsilon0*epsilonr0;
}

double prefactor(double epsilonr0){
  return 1.0/(1+sigma*mu0/(4*epsilonr0));
}


//----------Clase Lattice-Boltzmann----------
class LatticeBoltzmann{
private:
  int V[3][3][4], V0[3];//V[xyz][p][i]
  double e[3][3][4][2], e0[3]; //e[xyz][p][i][j]
  double b[3][3][4][2], b0[3];//b[xyz][p][i][j]
  double f[Lx][Ly][Lz][2][3][4][2],fnew[Lx][Ly][Lz][2][3][4][2];//f[ix][iy][iz][r][p][i][j]
  double f0[Lx][Ly][Lz][2],f0new[Lx][Ly][Lz][2];//f0[ix][iy][iz][r]
public:
  LatticeBoltzmann(void);
  void Inicie(void);
  void Colisione(double t);
  void Adveccione(void);
  double feq(double Bx0,double By0,double Bz0,double Eprimax0,double Eprimay0,double Eprimaz0,int r,int p,int i,int j,double epsilonr0);
  double feq0(double rho0);
  double rho(int ix,int iy,int iz);
  double Bx(int ix,int iy,int iz);
  double By(int ix,int iy,int iz);
  double Bz(int ix,int iy,int iz);
  double Eprimax(int ix,int iy,int iz);
  double Eprimay(int ix,int iy,int iz);
  double Eprimaz(int ix,int iy,int iz);
  double Jprimax(int ix,int iy,int iz){return sigma*Eprimax(ix,iy,iz);};
  double Jprimay(int ix,int iy,int iz){return sigma*Eprimay(ix,iy,iz);};
  double Jprimaz(int ix,int iy,int iz){return sigma*Eprimaz(ix,iy,iz);};
  void Muestre(void);
};

LatticeBoltzmann::LatticeBoltzmann(void){
  int ix,iy,iz,alpha,r,p,i,j;
  //Los vectores Velocidad V[alpha][p][i]=V^p_i_(alpha)
  V0[0]=V0[1]=V0[2]=0;

  V[0][0][0]=V[0][1][0]=V[1][2][0]=1;
  V[1][0][0]=V[2][1][0]=V[2][2][0]=1;
  V[2][0][0]=V[1][1][0]=V[0][2][0]=0;

  V[0][0][1]=V[0][1][1]=V[1][2][1]=-1;
  V[1][0][1]=V[2][1][1]=V[2][2][1]=1;
  V[2][0][1]=V[1][1][1]=V[0][2][1]=0;

  V[0][0][2]=V[0][1][2]=V[1][2][2]=-1;
  V[1][0][2]=V[2][1][2]=V[2][2][2]=-1;
  V[2][0][2]=V[1][1][2]=V[0][2][2]=0;

  V[0][0][3]=V[0][1][3]=V[1][2][3]=1;
  V[1][0][3]=V[2][1][3]=V[2][2][3]=-1;
  V[2][0][3]=V[1][1][3]=V[0][2][3]=0;
  //Los vectores Electricos e[alpha][p][i][j]=e^p_{ij}_(alpha)
  e0[0]=e0[1]=e0[2]=0;
  for(alpha=0;alpha<3;alpha++)
    for(p=0;p<3;p++)
      for(i=0;i<4;i++){
	e[alpha][p][i][0]=V[alpha][p][(i+1)%4]*0.5;
	e[alpha][p][i][1]=V[alpha][p][(i+3)%4]*0.5;
  }  
  //Los vectores Magneticos e[alpha][p][i][j]=b^p_{ij}_(alpha)
  b0[0]=b0[1]=b0[2]=0;  
  for(p=0;p<3;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++){
	b[0][p][i][j]=V[1][p][i]*e[2][p][i][j]-V[2][p][i]*e[1][p][i][j];
	b[1][p][i][j]=V[2][p][i]*e[0][p][i][j]-V[0][p][i]*e[2][p][i][j];
	b[2][p][i][j]=V[0][p][i]*e[1][p][i][j]-V[1][p][i]*e[0][p][i][j];
  }
}


void LatticeBoltzmann::Inicie(void){
  int ix,iy,iz,r,p,i,j; double E0=0.001,alpha=0.01,z0=40,c=1/sqrt(2); double B0=E0/c;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++)
	for(r=0;r<2;r++)
	  for(p=0;p<3;p++)
	    for(i=0;i<4;i++)
	      for(j=0;j<2;j++)
		f[ix][iy][iz][r][p][i][j]=fnew[ix][iy][iz][r][p][i][j]=
		  feq(0,B0*exp(-alpha*(iz-z0)*(iz-z0)),0,E0*exp(-alpha*(iz-z0)*(iz-z0)),0,0,r,p,i,j,epsilonr(ix,iy,iz));
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++)
	for(r=0;r<2;r++)
	  f0[ix][iy][iz][r]=f0new[ix][iy][iz][r]=feq0(-E0*2*alpha*(iz-z0)*exp(-alpha*(iz-z0)*(iz-z0)));
}


void LatticeBoltzmann::Colisione(double t){
  double rho0,Bx0,By0,Bz0,Eprimax0,Eprimay0,Eprimaz0,epsilonr0;
  int ix,iy,iz,r,p,i,j;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
	//Calculo las variables macroscopicas
	rho0=rho(ix,iy,iz);  Bx0=Bx(ix,iy,iz);  By0=By(ix,iy,iz);  Bz0=Bz(ix,iy,iz);
	Eprimax0=Eprimax(ix,iy,iz);  Eprimay0=Eprimay(ix,iy,iz);  Eprimaz0=Eprimaz(ix,iy,iz);
	epsilonr0=epsilonr(ix,iy,iz);
	//Hago la evolucion de la ec. de Boltzmann BGK
	for(r=0;r<2;r++)
	  for(p=0;p<3;p++)
	    for(i=0;i<4;i++)
	      for(j=0;j<2;j++)
		fnew[ix][iy][iz][r][p][i][j]=f[ix][iy][iz][r][p][i][j]-1.0/Tau*(f[ix][iy][iz][r][p][i][j]
										-feq(Bx0,By0,Bz0,Eprimax0,Eprimay0,Eprimaz0,r,p,i,j,epsilonr0));
	for(r=0;r<2;r++)
	  f0new[ix][iy][iz][r]=f0[ix][iy][iz][r]-1.0/Tau*(f0[ix][iy][iz][r]-feq0(rho0));
      }
}

void LatticeBoltzmann::Adveccione(void){
  int ix,iy,iz,r,p,i,j;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
	for(r=0;r<2;r++)
	  for(p=0;p<3;p++)
	    for(i=0;i<4;i++)
	      for(j=0;j<2;j++)
		f[(ix+V[0][p][i]+Lx)%Lx][(iy+V[1][p][i]+Ly)%Ly][(iz+V[2][p][i]+Lz)%Lz][r][p][i][j]=
		  fnew[ix][iy][iz][r][p][i][j];
	for(r=0;r<2;r++)
	  f0new[(ix+V[0][p][i]+Lx)%Lx][(iy+V[1][p][i]+Ly)%Ly][(iz+V[2][p][i]+Lz)%Lz][r]=
	    f0[ix][iy][iz][r];
      }
}

double LatticeBoltzmann::feq(double Bx0,double By0,double Bz0,double Eprimax0,double Eprimay0,double Eprimaz0,int r,int p,int i,int j,double epsilonr0){
  double VdotJ,Edote,Bdotb;
  double Jprimax0=sigma*Eprimax0,Jprimay0=sigma*Eprimay0,Jprimaz0=sigma*Eprimaz0;
  VdotJ=V[0][p][i]*Jprimax0+V[1][p][i]*Jprimay0+V[2][p][i]*Jprimaz0;
  Edote=e[0][p][i][j]*Eprimax0+e[1][p][i][j]*Eprimay0+e[2][p][i][j]*Eprimaz0;
  Bdotb=b[0][p][i][j]*Bx0+b[1][p][i][j]*By0+b[2][p][i][j]*Bz0;
  if(p==0)
    return 0.25*(0.25*VdotJ+epsilon(epsilonr0)*Edote+0.5/mu*Bdotb);
  if(p==1)
    return 0.25*(0.25*VdotJ+Edote+0.5*Bdotb);
}

double LatticeBoltzmann::feq0(double rho0){
  return rho0;
}


double LatticeBoltzmann::rho(int ix,int iy,int iz){
  int p,i,j; double suma=f0[ix][iy][iz][0];
  for(p=0;p<3;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++)
	suma+=f[ix][iy][iz][0][p][i][j];
  return suma;
}

double LatticeBoltzmann::Bx(int ix,int iy,int iz){
  int p,i,j; double suma=0;
  for(p=0;p<3;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++)
	suma+=f[ix][iy][iz][1][p][i][j]*b[0][p][i][j];
  return suma;
}

double LatticeBoltzmann::By(int ix,int iy,int iz){
  int p,i,j; double suma=0;
  for(p=0;p<3;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++)
	suma+=f[ix][iy][iz][1][p][i][j]*b[1][p][i][j];
  return suma;
}

double LatticeBoltzmann::Bz(int ix,int iy,int iz){
  int p,i,j; double suma=0;
  for(p=0;p<3;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++)
	suma+=f[ix][iy][iz][1][p][i][j]*b[2][p][i][j];
  return suma;
}

double LatticeBoltzmann::Eprimax(int ix,int iy,int iz){
  int p,i,j; double suma=0;
  for(p=0;p<3;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++)
	suma+=f[ix][iy][iz][0][p][i][j]*e[0][p][i][j];
  return prefactor(epsilonr(ix,iy,iz))*suma/epsilonr(ix,iy,iz);
}

double LatticeBoltzmann::Eprimay(int ix,int iy,int iz){
  int p,i,j; double suma=0;
  for(p=0;p<3;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++)
	suma+=f[ix][iy][iz][0][p][i][j]*e[1][p][i][j];
  return prefactor(epsilonr(ix,iy,iz))*suma/epsilonr(ix,iy,iz);
}

double LatticeBoltzmann::Eprimaz(int ix,int iy,int iz){
  int p,i,j; double suma=0;
  for(p=0;p<3;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++)
	suma+=f[ix][iy][iz][0][p][i][j]*e[2][p][i][j];
  return prefactor(epsilonr(ix,iy,iz))*suma/epsilonr(ix,iy,iz);
}

void LatticeBoltzmann::Muestre(void){
  int iz; double E2,B2;
  for(iz=0;iz<Lz;iz++){
    E2=pow(Eprimax(0,0,iz),2)+pow(Eprimay(0,0,iz),2)+pow(Eprimaz(0,0,iz),2);
    B2=pow(Bx(0,0,iz),2)+pow(By(0,0,iz),2)+pow(Bz(0,0,iz),2);
    cout<<iz<<" "<<0.5*(epsilonr(0,0,iz)*E2+B2/mur)<<endl;
  }
}

int main(){
  LatticeBoltzmann Pulso;
  int t, tmax=180;
  
  Pulso.Inicie();
  for(t=0;t<tmax;t++){
    Pulso.Colisione(t);
    Pulso.Adveccione();
  }
  Pulso.Muestre();

  return 0;
}
