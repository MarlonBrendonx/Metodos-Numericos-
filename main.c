#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

/**

  Trabalho2:Parte 2 -Realiza a solução do sistema Ax=b utilizando o método da eliminação de Gauss e Gauss-Jacobi

  A = |  4  1 −2 |   b = |3|
      |  2  3  0 |       |5|
      |  2 −2  6 |       |1|



  Autor       : Marlon Brendo Ramos
  Instituição : Universidade Federal de Uberlândia
  Curso       : Sistemas de Informação
  Disciplina  : Métodos Numéricos
  Tecnologias : Linguagem C e IDE codeblocks 17.12-Compilado em Debian Buster 10.5


  --------------------------------------------------------------------------------------------------------------------------
  Esse programa permite a resolução do sistema Ax=b utilizando o método da eliminação de Gauss e Gauss-Jacobi.
  O programa permite a opção de escolher a solução desejada e retornar as etapas dos processos iterativos da solução ou as
  matrizes correspondentes na saída padrão(STDOUT)



]




**/

float **alocamatriz(float vet[3][4]){


      float **matriz;  /* ponteiro para a matriz */
      int  i;    /* variavel auxiliar      */

      /* aloca as linhas da matriz */
      matriz = (float **) calloc (3, sizeof(float *));
      if (matriz == NULL) {
         printf ("** Erro: Memoria Insuficiente **");
         return (NULL);
      }

      /* aloca as colunas da matriz */
      for ( i = 0; i < 3; i++ ) {
          matriz[i] = (float*) calloc (4, sizeof(float));	/* m vetores de n floats */
          if (matriz[i] == NULL) {
             printf ("** Erro: Memoria Insuficiente **");
             return (NULL);
          }

       }

        for ( i = 0; i < 3; i++ ) {
            for (int  j = 0; j < 4; j++ ) {

                    matriz[i][j]=vet[i][j];
            }

        }



      return (matriz); /* retorna o ponteiro para a matriz */




}




void mostrarMatriz(float sistema[3][4]){

        for(int i=0;i<3;i++){

            for(int j=0;j<4;j++){

                    printf(" %.1f  ",sistema[i][j]);

            }
            printf("\n\n");
        }

}



int metodoGauss(float **sistema,int *retorno){
            int i,j,k;
            float m,y,aux;

            for(k=0;k<3;k++){
                    for(i=k+1;i<3;i++){

                       m=sistema[i][k]/sistema[k][k];
                       sistema[i][k]=0;

                       for(j=k+1 ;j<4;j++){

                            sistema[i][j]=sistema[i][j] - (m*sistema[k][j]);

                       }


                    }

            }
            printf("\nMatriz triangular superior:\n\n");
            for(int k=0;k<3;k++){

                    for(int j=0;j<4;j++){

                            printf(" %.2f ",sistema[k][j]);
                    }
                    printf("\n");
            }

            float s=0;
            j=0;
            float x[4];
            //Calcula o x3
            x[2]=sistema[2][3]/sistema[2][2];
            aux=x[2];
            for(k=2;k>=0;k--){
                    s=0;
                    for(j=k+1;j<4;j++){

                       s=s+ ( sistema[k][j] * x[j] );

                    }
                      x[k]=(sistema[k][3]-s)/sistema[k][k];
                      aux=x[k];

            }

            printf("\nAs soluções são: \n");
            for(int i=0; i<3; i++)
            {
                printf("\nx%d=%.3f\t",i,x[i]); /* x1, x2, x3 são as soluções do sistema*/
            }



}

float encontraMaiorvalor(float *anterior,float *proximo){

     int i,j,minant,minprox;
     float aux;

     //Coloca os vetores em módulo
     for(int i=0;i<3;i++){

        anterior[i] = fabs( anterior[i] );
        proximo[i]  = fabs( proximo[i] );

     }

     //Ordera os vetores
     for(i=0; i<3 ;i++){
            minant=i;
            minprox=i;

        for(j=i+1;j<3;j++){

              if(anterior[j] > anterior[minant])
                minant=j;

              if(proximo[j] > proximo[minprox])
                minprox=j;
        }

        if( i!=minant ){

              aux=anterior[i];
              anterior[i]=anterior[minant];
              anterior[minant]=aux;

        }
        if( i!= minprox){

             aux=proximo[i];
             proximo[i]=proximo[minprox];
             proximo[minprox]=aux;

        }

     }

     if( anterior[0] > proximo[0] )
        return anterior[0];

     return proximo[0];

}


float criterioParadaGaussjacobi(float *anterior,float *proximo){

            float diferenca;
            float maiordiferenca=0,Drk=0;
            int i=0;

            //Encontra max|X^k+1 - X^k|
            for( i=0; i<3 ; i++ ){

                    diferenca=fabs( proximo[i] - anterior[i] );

                    if( i==0 ){

                        maiordiferenca=diferenca;

                    }else if(diferenca > maiordiferenca){

                        maiordiferenca=diferenca;

                    }

            }

            //Encontra max|X^k+1|

           Drk=maiordiferenca/encontraMaiorvalor(anterior,proximo);

          return Drk;

}

float * insereVetoranterior(float x1,float x2,float x3){

            float *anterior=(float *)calloc(3,sizeof(float) );
            anterior[0] = x1;
            anterior[1] = x2;
            anterior[2] = x3;

            return anterior;
}

float * insereVetorproximo(float x1,float x2,float x3){

            float *proximo=(float *)calloc(3,sizeof(float) );
            proximo[0] = x1;
            proximo[1] = x2;
            proximo[2] = x3;

            return proximo;
}


void metodoGaussjacobi(float sistema[3][4]){


        float x1,x2,x3,Drk;
        float vetorinicial[3][1]={ 1.0,0.0,0.0 },precisao=0.001;
        float *anterior,*proximo;
        int i=0;


        printf("Isolando x1,x2,x3\n\n");
        printf("x1 = 0.75 - 0.25x2 + 0.5x3\n");
        printf("x2 = 1.66 - 0.66x1\n");
        printf("x3 = 0.16 - 0.33x1 + 0.33x2\n");

        //Isolando as váriavies x1,x2,x3
        x1=0.75 - ( 0.25*vetorinicial[1][0] ) + ( 0.50*vetorinicial[2][0]);
        x2=1.66 - ( 0.66*vetorinicial[0][0]);
        x3=0.16 - ( 0.33*vetorinicial[0][0] ) + (0.33*vetorinicial[1][0] );


        do{

                    anterior=insereVetoranterior(x1,x2,x3);
                    x1=0.75 - ( 0.25 * anterior[1]) + ( 0.5 * anterior[2]);
                    x2=1.66 - ( 0.66 * anterior[0]);
                    x3=0.16 - ( 0.33 * anterior[0] ) + (0.33 * anterior[1] );
                    proximo=insereVetorproximo(x1,x2,x3);

                    printf("---------------------Iteracao k = %d --------------------\n",i);
                    printf("X^(%d) = [%.3f,%.3f,%.3f]\n",i,anterior[0],anterior[1],anterior[2]);
                    printf("X1 = %.3f\nX2 = %.3f\nX3 = %.3f\n",x1,x2,x3);
                    printf("X^(%d) = [ %.3f,%.3f,%.3f]\n",i+1,proximo[0],proximo[1],proximo[2]);
                    Drk=criterioParadaGaussjacobi(anterior,proximo);
                    printf("Drk = %.3f\n",Drk);

                    printf("----------------------------------------------------------\n");

        i++;
        }while(Drk > precisao);

        printf("%d iteracoes-\nRaizes aproximadas:\n",i);
        printf("x1 = %.3f\n",x1);
        printf("x2 = %.3f\n",x2);
        printf("x3 = %.3f\n",x3);

}



int menu(){

     int op;
     printf("\n\nResolução do sistema: Ax=b\n\n");
     printf("A = 4   1  -2      b= 3  \n");
     printf("    2   3   0         5  \n");
     printf("    2  -2   6         1  \n\n");

     printf("[1]-Método da eliminacao de Gauss\n\n");
     printf("[2]-Método de Gauss-Jacobi \n\n");
     printf("[3]-Sair\n");
     scanf("%d",&op);


     return op;


}

int * alocaVetorretorno(){

         int *retorno=(int *)calloc(2,sizeof(int));
         retorno[0]=0;
         retorno[1]=0;

         return retorno;
}

int main()
{
    int op;
    float vet[3][4]={4,1,-2,3,2,3,0,5,2,-2,6,1};
    int *retorno=alocaVetorretorno();
    float **sistema=alocamatriz(vet);


    do{

        op=menu();
        system("clear");

        switch(op){


                case 1:
                metodoGauss(sistema,retorno);
                break;
                case 2:
                metodoGaussjacobi(sistema);
                break;

                default:
                if(op!=3)
                printf("Opção incorreta!\n");
                break;
        }

    }while(op!=3);



    return 0;
}
