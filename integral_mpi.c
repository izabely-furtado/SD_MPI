// Author: Izabely Furtado
//
// Program that computes the integral of an array of elements in parallel using
// MPI_Reduce.
//
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/sysinfo.h>

//double integral = 0.0;


 //função que calcula a area de um trapézio
double AreaTrapezio(double dx, double h1, double h2){
	double area;
	area = dx * (h1+h2)/2;
	return area;
}

 //função qualquer f(x)
double f(double x){
	return (4*sqrt(1-x*x));
}

/*

 //criando limite
typedef struct limite {
	double inferior;
	double superior;
	long n;
}Limite;

 //construtor do struct
Limite *CriaLim(double inferior, double superior, long numero){
	Limite *l = (Limite *) malloc(sizeof(Limite));
	l->inferior = inferior;
	l->superior = superior;
	l->n = numero;
	return l;
}
*/

 //calcula a area da integral da função f(x)
double CalculaArea(double a, double b, int N){
	int i;
	double area, x1, x2, f1, f2 = 0.0;
	double dx  = (b-a)/N;
	for(i = 0; i<N; i++){
		x1 = a + dx * i;
		x2 = a + dx * (i + 1);
		if(x1>1 || x2>1){
			printf("%lf, %lf \n", x1, x2);
		}
		f1 = f(x1);
		f2 = f(x2);
		area += AreaTrapezio(dx, f1, f2);
	}
	return area;
}


//ponteiro pra função ...
//dividindo f(x) para ser calculada por partes... (estilo Jack)
double RankCalculaArea(double A, double B, long int nDivi, int rank, int size){
	//long nDivi;  		//indica o numero de divisões na integral
	//int A, B;    		//limite inferior e superior 
	int nLocal;  		//indica o numero de divisões feitas por bloco
	double larguraLocal;   //indica a largura processada por cada bloco
	double aLocal;         //indica o limite inferior processado por cada bloco
	double bLocal;         //indica o limite superior processado por cada bloco
	
	larguraLocal = (B - A) / size;
	aLocal = A + rank * larguraLocal;
	bLocal = aLocal + larguraLocal;
	nLocal = nDivi/size;
	
	double integralRank = CalculaArea(aLocal, bLocal, nLocal);
	
	return integralRank;
}


int main(int argc, char **argv) {
	MPI_Init(&argc, &argv); //????
		
	
	int size, rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	double global_sum;
	double local_sum;
		
	
	//local_sum = RankCalculaArea(A, B, nDivi, rank, size);
	local_sum = RankCalculaArea(0.0, 1.0, 2000, rank, size);
	
	MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	
	//Imprime o resultado
	if (rank == 0) {
		printf("Total sum = %f \n", global_sum);
	}

	MPI_Finalize();
	return 0;

}
