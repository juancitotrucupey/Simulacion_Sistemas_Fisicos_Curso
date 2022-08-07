#include <iostream>
#include <cmath>

using namespace std;

const int Lx=200;
const double p=0.5;


class LatticeGas{
	private:
		int V[2];
		double n[Lx][2], nnew[Lx][2]; //int n[ix,i], nnew[ix,i]; ix=la posicion de la celda, i=la direccion
	public:
	         LatticeGas(void);//constructor
                void Inicie(double mu, double sigma);
	        void Show(void);
                void Colisione(void);
		void Adveccione(void);
                double Getsigma2(void);
};

LatticeGas::LatticeGas(void){
	V[0]=1; V[1]=-1;
}
void LatticeGas::Inicie(double mu, double sigma){
	for(int ix=0;ix<Lx;ix++)
	n[ix][0]=n[ix][1]=0.5/(sigma*sqrt(2*M_PI))*exp(-0.5*pow((ix-mu)/sigma,2));
}
  void LatticeGas::Show(void){
		for(int ix=0;ix<Lx;ix++)
		  cout<<ix<<" "<<n[ix][0]+n[ix][1]<<endl;
		}
void LatticeGas::Colisione(void){
	for (int ix=0; ix<Lx; ix++){
	  nnew[ix][0]=p*n[ix][0]+(1-p)*n[ix][1];
	  nnew[ix][1]=p*n[ix][1]+(1-p)*n[ix][0];		
		}
	}

void LatticeGas::Adveccione(void){
	for(int ix=0;ix<Lx;ix++){
		for(int i=0;i<2;i++){
			n[(ix+V[i]+Lx)%Lx][i]=nnew[ix][i];
		}
	}
}
  double LatticeGas::Getsigma2(void){
    int ix;	double Norma, Xprom, s2;
        //Calcule Norma
    for(Norma=0,ix=0;ix<Lx;ix++)  
      {
	Norma+=(n[ix][0]+n[ix][1]);
      }
      //Calcular Xprom
	for(Xprom=0, ix=0; ix<Lx; ix++)
	  Xprom+=ix*(n[ix][0]+n[ix][1]);
	Xprom/=Norma;
	//Calcular s^2
	for(s2=0, ix=0; ix<Lx; ix++)
	  s2+=pow(ix-Xprom,2)*(n[ix][0]+n[ix][1]);
	s2/=Norma;
	//Devuelva el resultado
	return s2;

  }
  //______________________________Funciones Globales_________________________________________________
  

int main(void){
	LatticeGas Difusion;
	int t,tmax=400;
	
	Difusion.Inicie(Lx/2,5.0);
	for(t=0;t<tmax;t++){
	Difusion.Colisione();
	Difusion.Adveccione();
	cout<<t<<" "<<Difusion.Getsigma2()<<endl;
	}
	return 0;
}
