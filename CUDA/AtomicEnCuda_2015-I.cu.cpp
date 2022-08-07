#include<iostream>
#include <fstream>
#include <cmath>
#include <GL/glew.h>
#include <GL/glut.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>
using namespace std;

#define Lx 128
#define Nx 64
const int Mx=(Lx+Nx-1)/Nx;


//--------------- KERNELS ----------------
__global__ void SumarNuestrosIndices(float * d_test){
  int ix; ix=blockIdx.x*blockDim.x+threadIdx.x;
  //  d_test[0]+=ix;
  atomicAdd(&d_test[0],ix);
}

int main(){
  //DECLARAR LAS MATRICES
  float h_test[1];
  float *d_test;  cudaMalloc((void**) &d_test,sizeof(float));

  //INICIALIZAR LOS DATOS
  //Cargarlos en el Host
  h_test[0]=0;
  //Enviarlos al Device
  cudaMemcpy(d_test,h_test,sizeof(float),cudaMemcpyHostToDevice);
  

  //PROCESAR EN LA TARJETA GRAFICA
  dim3 ThreadsPerBlock(Nx,1,1);
  dim3 BlocksPerGrid(Mx,1,1);
  SumarNuestrosIndices<<<BlocksPerGrid,ThreadsPerBlock>>>(d_test);

  //IMPRIMIR LOS DATOS
  //Devolverlos al Host
  cudaMemcpy(h_test,d_test,sizeof(float),cudaMemcpyDeviceToHost);
  //Imprimirlos
  cout<<h_test[0]<<endl;

  //LIBERAR MEMORIA
  cudaFree(d_test);

  return 0;
}
