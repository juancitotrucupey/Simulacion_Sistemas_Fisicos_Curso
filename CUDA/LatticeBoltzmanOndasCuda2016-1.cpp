#include<iostream>
#include <fstream>
#include <cmath>
#include <GL/glew.h>
#include <GL/glut.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>
using namespace std;

const double C=0.5;
const double C2=C*C;
const double UNO_2C2=1-2*C2;
const double Tau=05;

const double A=10;
const double Lambda=10;
const double Omega=2*M_PI*C/Lambda;

#define Lx 256
#define Ly 256
#define Nx 8
#define Ny 8
const int Mx=(Lx+Nx-1)/Nx;
const int My=(Ly+Ny-1)/Ny;

//--------------------KERNELS----------------

__constant__ float d_w[5];
__constant__ float v_x[5];
__constant__ float v_y[5];

__device__ float SumeleUno(float x){
  return x+1;
}
__global__ void SumarConstante(float * d_a,size_t pitcha){
  int ix,iy; float *aux;
  ix=blockIdx.x*blockDim.x+threadIdx.x;  iy=blockIdx.y*blockDim.y+threadIdx.y;
  
  aux=d_a+(ix*pitcha)/sizeof(float)+iy; // aux es &(d_a[ix][iy])
  
  (*aux)=SumeleUno((*aux)); //  (*aux) es d_a[ix][iy]
}

//------------------- CLASES ----------------
class LatticeBoltzmann{
private:
  //Los pesos
  float h_w[5];
  //Los vectores
  int h_Vx[5],h_Vy[5];
  //Las funciones f
  float h_f0[Lx][Ly]; float *d_f0; size_t pitchf0;
   float h_f1[Lx][Ly]; float *d_f0; size_t pitchf1;
    float h_f2[Lx][Ly]; float *d_f0; size_t pitchf2;
    float h_f3[Lx][Ly]; float *d_f0; size_t pitchf3;
    float h_f4[Lx][Ly]; float *d_f0; size_t pitchf4;
  //Las funciones fnew
    float h_f0new[Lx][Ly]; float *d_f0new; size_t pitchf0new;
    float h_f1new[Lx][Ly]; float *d_f1new; size_t pitchf1new;
    float h_f2new[Lx][Ly]; float *d_f2new; size_t pitchf2new;
    float h_f3new[Lx][Ly]; float *d_f3new; size_t pitchf3new;
    float h_f4new[Lx][Ly]; float *d_f4new; size_t pitchf4new;

}
 
public:
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  void Inicie(void);
  void Incremente(void);
  void Muestre(void);
void Colisinone(int t);
};
LatticeBoltzmann::LatticeBoltzmann(void){
  cudaMallocPitch((void**) &d_a,&pitcha,Ly*sizeof(float),Lx);
 //Cargar los vectores velocidad
  h_Vx[0]=0; h_Vy[0]=0;
		h_Vx[1]=1; h_Vy[2]=0; 	h_Vx[3]=-1; h_Vy4]=0;
		h_Vx[1]=0; h_Vy[2]=1;	h_Vx[3]=0;  h_Vy[4]=-1;
  // Cargar los pesos				
  h_w[0]=1/3; h_w[1]=h_w[2]=h_w[3]=h_w[4]=1/6;
  //Enviar al device
  cudaMencpy(d_w,h_w,5*sizeof(float),cudaMencpyHostToDevice);
  cudaMencpy(d_Vx,h_Vx,5*sizeof(float),cudaMencpyHostToDevice);
  cudaMencpy(d_Vy,h_Vy,5*sizeof(float),cudaMencpyHostToDevice);


}
LatticeBoltzmann::~LatticeBoltzmann(void){
  cudaFree(d_f0); cudaFree(d_f1); cudaFree(d_f2); cudaFree(d_f3); cudaFree(d_f4);
   cudaFree(d_f0new); cudaFree(d_f1new); cudaFree(d_f2new); cudaFree(d_f3new); cudaFree(d_f4new);
}
void LatticeBoltzmann::Inicie(void){
  
}
void LatticeBoltzmann::Incremente(void){
  dim3 ThreadsPerBlock(Nx,Ny,1); //x,y,z
  dim3 BlocksPerGrid(Mx,My,1);
  SumarConstante<<<BlocksPerGrid,ThreadsPerBlock>>>(d_a,pitcha);
}
void LatticeBoltzmann::Muestre(void){
  //Devolverlos al Host
  cudaMemcpy2D(h_a,Ly*sizeof(float),d_a,pitcha,Ly*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  //Imprimirlos
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++)
      cout<<h_a[ix][iy]<<" ";
    cout<<endl;
  }
  cout<<endl;
}
void LatticeBoltzmann::Colisione(int t){
  // Procesar en el device
  dim3 ThreadsPerBolck(Nx,Ny,1);
  dim3 BlocksPerGrid(Mx,My,1);
  d_Colisione<<<BlocksPerGrid


}
//---------------- FUNCIONES GLOBALES Y PROGRAMA PRINCIPAL -------
int main(void){
  LatticeBoltzmann Ondas;

  Ondas.Inicie();
  Ondas.Muestre();
  Ondas.Incremente();
  Ondas.Muestre();

  return 0;
}
