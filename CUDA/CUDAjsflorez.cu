#include <iostream>
#include <fstream>
#include <cmath>
#include <GL/glew.h>
#include <Gl/glut.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>
using namespace std;


__global__ void Test(float *d_test){
  (*(d_test+0))=5.0;

}

int main(void){
  //Declarar las matrices

  float h_test[1];
  float *d_test; cudaMalloc((void**)&d_test,sizeof(float));

  //Inicializar los datos
  //cargar los datos en el host
  h_test[0]=2.5;
  //Enviar al device
  cudaMencpy(d_test,h_test,sizeof(float),cudaMencpyHostToDevice);


  //Procesar en la tarjeta grafica
  dim3 ThreadsPerBlock(1,1,1);
  dim3 BlocksPerGrid(1,1,1);
  Test<<<BlocksPerGrid,ThreadsPerGrid>>>(d_test);

  //Devolver e imprimir
  //Devolverlos al Host
   cudaMencpy(d_test,h_test,sizeof(float),cudaMencpyDeviceToHost);
   //Imprimirlos
   cout<<h_test[0]<<endl;

  return 0;

}
