// Suma de dos vectores c=a+b en CUDA
#include<iostream>
#include <fstream>
#include <cmath>
#include <GL/glew.h>
#include <GL/glut.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>
using namespace std;

#define Lx 32
#define Nx 8
const int Mx=(Lx+Nx-1)/Nx;

//--------------- KERNELS ----------------
__global__ void SumaDeATres(float * d_a){
  __shared__ float Shared[Nx+2];
  int ix; ix=blockIdx.x*blockDim.x+threadIdx.x;
  int is; is=threadIdx.x+1;

  //Copiar los datos a la memoria compartida
  switch(is){
  case 1:    Shared[is]=d_a[ix]; Shared[is-1]=d_a[(ix-1+Lx)%Lx]; break; 
  case Nx-1: Shared[is]=d_a[ix]; Shared[is+1]=d_a[(ix+1)%Lx]; break; 
  default:   Shared[is]=d_a[ix];
  }
  __syncthreads();

  //Procesamos
  float suma=Shared[is-1]+Shared[is]+Shared[is+1];

  __syncthreads();

  //Colocar los resultados en la memoria general
  d_a[ix]=suma;
}

int main(){
  int ix;
  //DECLARAR LAS MATRICES
  //En el Host
  float h_a[Lx];
  //En el Device
  float *d_a;  cudaMalloc((void**) &d_a,Lx*sizeof(float));

  //INICIALIZAR LOS DATOS
  //Cargarlos en el Host
  for(ix=0;ix<Lx;ix++) h_a[ix]=ix;
  //Enviarlos al Device
  cudaMemcpy(d_a,h_a,Lx*sizeof(float),cudaMemcpyHostToDevice);
  //Imprimirlos
  for(ix=0;ix<Lx;ix++) cout<<h_a[ix]<<" "; cout <<endl;
  
  //PROCESAR EN LA TARJETA GRAFICA
  dim3 ThreadsPerBlock(Nx,1,1);
  dim3 BlocksPerGrid(Mx,1,1);
  SumaDeATres<<<BlocksPerGrid,ThreadsPerBlock>>>(d_a);

  //IMPRIMIR LOS DATOS
  //Devolverlos al Host
  cudaMemcpy(h_a,d_a,Lx*sizeof(float),cudaMemcpyDeviceToHost);
  //Imprimirlos
  for(ix=0;ix<Lx;ix++) cout<<h_a[ix]<<" "; cout <<endl;

  //LIBERAR MEMORIA
  cudaFree(d_a);

  return 0;
}
