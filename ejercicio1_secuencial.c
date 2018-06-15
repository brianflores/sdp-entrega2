#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

double dwalltime()
{
  double sec;
  struct timeval tv;

  gettimeofday(&tv,NULL);
  sec = tv.tv_sec + tv.tv_usec/1000000.0;
  return sec;
}

double *A, *B, *C, *D, *L, *U, *AB, *LC, *DU, *TOTAL;
int main(int argc,char*argv[]){
    if (argc < 2){
        printf("Faltan argumentos \n");
        return 0;
    }
    int i, j, k;
    int N = atol(argv[1]);
    unsigned long Total = N*N;
    int check = 1;
    double promedioL, promedioU, timetick;
    int elementosU = (N*N)-((N*(N-1))/2);
    A=(double*)malloc(sizeof(double)*N*N);
    B=(double*)malloc(sizeof(double)*N*N);
    C=(double*)malloc(sizeof(double)*N*N);
    D=(double*)malloc(sizeof(double)*N*N);
    L=(double*)malloc(sizeof(double)*N*N);
    U=(double*)malloc(sizeof(double)*elementosU);
    for(i=0;i<N;i++){       //Crea matrices
       for(j=0;j<N;j++){
           A[i*N+j]=1.0;
           B[i*N+j]=1.0;
           C[i*N+j]=1.0;
           D[i*N+j]=1.0;
           if(i==j){
               L[i*N+j]= 1.0;
               //U[i*N+j]= 1.0;
           } else if(i>j){
               //U[i*N+j]= 1.0;
               L[i*N+j]= 0.0;
           } else {
               //U[i*N+j]= 0.0;
               L[i*N+j]= 1.0;
           }
       }
    }
    for(i=0; i<elementosU; i++) U[i] = 1.0;
    promedioL = 0;
    promedioU = 0;

    timetick = dwalltime();
    for(i=0;i<N;i++){   //Calcula los promedios
       for(j=0;j<N;j++){
           promedioL+= L[i*N+j];
           //promedioU+= U[i*N+j];
       }
    }
    for(i=0; i<elementosU; i++) promedioU+= U[i];
    promedioL = promedioL / Total;
    promedioU = promedioU / Total;
    promedioL = promedioL * promedioU; //En promedioL queda el promedio de L por el de U.
    AB=(double*)malloc(sizeof(double)*N*N); //AB=A*B
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            AB[i*N+j]=0;
            for(k=0;k<N;k++){
	            AB[i*N+j]= AB[i*N+j] + A[i*N+k]*B[k+j*N];
            }
        }
    }
    LC=(double*)malloc(sizeof(double)*N*N); //LC= L*C
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            LC[i*N+j]=0;
            for(k=i;k<N;k++){
	            LC[i*N+j]= LC[i*N+j] + L[i*N+k]*C[k+j*N];
            }
        }
    }
    DU=(double*)malloc(sizeof(double)*N*N); //DU= D*U
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            DU[i*N+j]=0;
            for(k=0;k<=j;k++){
	            DU[i*N+j]= DU[i*N+j] + D[i*N+k]*U[k+j*(j+1)/2];
            }
        }
    }
    TOTAL=(double*)malloc(sizeof(double)*N*N); //TOTAL= (AB+LC+DU)*promedioL
    for(i=0;i<N;i++){
        for(j=0;j<N;j++)
        TOTAL[i*N+j]= (AB[i*N+j] + LC[i*N+j] + DU[i*N+j])*promedioL;
    }
    printf("Tiempo en segundos %f \n", dwalltime() - timetick);

    /* double resultado = TOTAL[0];
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
        check = check && (TOTAL[i*N+j]==resultado);
        }
    } */
	double resultado;
	for(i=0;i<N;i++){
		resultado = TOTAL[i*N]/promedioL;
        for(j=0;j<N;j++){
        check = check && (TOTAL[i*N+j]/promedioL==resultado+j);
        }
	}
	for(j=0;j<N;j++){
		resultado = TOTAL[j*N]/promedioL;
		for(i=0;i<N;i++){
		    check = check && (((TOTAL[i+j*N]/promedioL)-resultado-i)==0);
        }
	}
    if(check){
        printf("Multiplicacion de matriz correcta\n");
    }else{
        printf("Multiplicacion de matriz erroneo\n");
    }

    free(A);
    free(B);
    free(C);
    free(D);
    free(L);
    free(U);
    free(AB);
    free(LC);
    free(DU);
    free(TOTAL);
}
