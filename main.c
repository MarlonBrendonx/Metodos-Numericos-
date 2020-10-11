#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

/**

 Trabalho 2-Realiza a resolução da equação não linear F(x) = 2cos(x) - 0.5e^x
 utilizando os métodos iterativos da Bissecção,Falsa Posição e Newton-Raphson.



 Autor       : Marlon Brendo Ramos
 Instituição : Universidade Federal de Uberlândia
 Curso       : Sistemas de Informação
 Disciplina  : Métodos Numéricos
 Tecnologias : C e IDE codeblocks 17.12-Compilado em Debian Buster 10.5

 -------------------------------------------------------------------------------
 Este programa permite a resolução da equação não linear F(x) = 2cos(x) - 0.5e^x
 por 3 métodos diferentes:Bissecção,Falsa Posição e Newton-Raphson.O programa re-
 cebe como parâmetro os intervalos [a,b] e verifica através do teorema de Bolzano
 F(a)*F(b) < 0 se há uma raíz dentro desse intervalo,caso haja,é realizada o algo-
 ritmo escolhido respectivamente pelo usuário,retornando o valor aproximado da raíz
 na saída padrão(STDOUT).


  Exemplo:

  [1]-Método da bissecção
  [2]-Método da Posição Falsa
  [3]-Método de Newton-Raphson
  Opção:1
  ...
  Informe o intervalo onde deve se buscar a raiz aproximada (ex[0,1] = digite 0 1)

  Intervalo[a,b]: 0 1

  Resultado da raíz aproximada pelo método da Bissecção

  X=...



**/


double equacao(double x){

            //Encontra F(x) = 2cos(x) - 0.5e^x para x
            return ( 2*cos(x) )-  0.5* ( pow(2.7182,x) );
}

double derivadaEquacao(double x){

            // F(x) = 2cos(x) - 0.5e^x
            // F(x)' = -2*sen(x) - 0.5e^x
            return  ( -2*sin(x) ) - 0.5* ( pow(2.7182,x) ) ;

}


double * valoresIntervalo(){

            double *inter=(double *)calloc(2,sizeof(double) );
            scanf("%lf %lf",&inter[0],&inter[1]);
            return inter;
}

double metodoNewton(double a,double b){

            double precisao=0.001; //precisão E=10^-1
            double Xn,Xnprox,Fx,Fxderivada;
            int flag=0;
            bool   bolzano=false;
            int *intervalos,i=0;

            /**Verifica  se satisfaz o teorema de Bolzano
                Se f(a)*f(b) < 0
            **/

            bolzano= equacao(a) * equacao(b) < 0 ? true : false;

            if( !bolzano )
                return 0;

             do{

                //|F'(a)| < |F'(b)|
                if(!flag){
                    if( derivadaEquacao(a) < derivadaEquacao(b)  ? true : false ){

                        Xn=b;

                    }else if( derivadaEquacao(b) < derivadaEquacao(a)  ? true : false  ){

                        Xn=a;
                    }

                    flag=1;
                }
                printf("-------Iteração k = %d--------\n\n",i);
                printf("Xn = %.3lf\n",Xn);

                //F(x)
                Fx=equacao(Xn);
                //F'(x)
                Fxderivada=derivadaEquacao(Xn);
                Xn=( Xn - (Fx / Fxderivada) );


                printf("F(Xn) = %.3lf\n",Fx);
                printf("F'(Xn) = %.3lf\n",Fxderivada);
                printf("Xn+1 = %.3lf\n",Xn);
                printf("\n\n---------------------------------\n\n");
                i++;

            }while( fabs(equacao(Xn)) > precisao  );


            printf("Raíz aproximada: %.6lf\n\n",Xn);






}

double metodoPosicaoFalsa(double a,double b){


            double precisao=0.001; //precisão E=10^-3
            double x;
            bool   bolzano=false;
            int *intervalos,i=0;

            /**Verifica  se satisfaz o teorema de Bolzano
                Se f(a)*f(b) < 0
            **/
            bolzano= equacao(a) * equacao(b) < 0 ? true : false;

            if( !bolzano )
                return 0;

             do{


                x= ( ( a*equacao(b) )  - ( b*equacao(a) ) )/(equacao(b)-equacao(a));

                if( equacao(a) * equacao(x) < 0 ? true : false ){

                    b=x;

                }else if( equacao(x) * equacao(b) < 0 ? true : false ){

                    a=x;
                }

                printf("-------Iteração k = %d--------\n\n",i);
                printf("Xk            = %.3lf\n",x);
                printf("Intervalo     = [%.3lf,%.3lf]\n",a,b);
                printf("| F(%.3lf) |  = %.6lf\n",x,fabs(equacao(x)));
                printf("bk - ak       = %.6lf\n",b-a);
                printf("\n\n---------------------------------\n\n");

                i++;

            }while( fabs(equacao(x)) > precisao  );


            printf("Intervalo:[%.3lf,%.3lf]\n",a,b);
            printf("Raíz aproximada: %.6lf\n\n",x);



}


double metodoBisseccao(double a,double b){


            double precisao=0.001; //precisão E=10^-3
            double x;
            bool   bolzano=false;
            int *intervalos;

            /**Verifica  se satisfaz o teorema de Bolzano
                Se f(a)*f(b) < 0
            **/
            bolzano= equacao(a) * equacao(b) < 0 ? true : false;

            if( !bolzano )
                return 0;

            if(fabs(b-a) < precisao)
                printf("%.lf\n",(a+b)/2);

            /** iterações **/
            int i=0;
            do{


                x=(a+b)/2;

                if( equacao(a) * equacao(x) < 0 ? true : false ){

                    b=x;

                }else if( equacao(x) * equacao(b) < 0 ? true : false ){

                    a=x;
                }

                printf("-------Iteração k = %d--------\n\n",i);
                printf("Xk            = %.3lf\n",x);
                printf("Intervalo     = [%.3lf,%.3lf]\n",a,b);
                printf("| F(%.6lf) |  = %.6lf\n",x,fabs(equacao(x)));
                printf("bk - ak       = %.3lf\n",b-a);
                printf("\n\n---------------------------------\n");

                i++;
            }while( fabs(equacao(x)) > precisao  );


            printf("Intervalo:[%.3lf,%.3lf]\n",a,b);
            printf("Raíz aproximada: %.6lf\n\n",x);


}

* verificaIntervalos(double *intervalos){

       do{

             system("clear");
             printf("\033[1;31mNao ha raizes dentro desse intervalo!\n\033[0m");
             printf("Informe um novo intervalo para buscar a raiz aproximada \033[1;34m(ex[0,1] = digite 0 1)\n\n\033[0m");
             intervalos=valoresIntervalo();

       }while( !(equacao(intervalos[0]) * equacao(intervalos[1]) < 0 ? true : false) );

       return intervalos;


}


int menu(){

     int op;
     printf("Resolução da equaçâo não-linear: F(x) = 2cos(x) - 0.5e^x\n\n\n");
     printf("[1]-Método da bissecção\n\n");
     printf("[2]-Método da Posição Falsa\n\n");
     printf("[3]-Método de Newton-Raphson\n\n");
     printf("[4]-Sair\n");
     scanf("%d",&op);


     return op;


}

int main()
{

    int op;
    double a,b;
    double *intervalos;
    bool bolzano=false;
    do{

        op=menu();
        system("clear");

        switch(op){


                case 1:
                system("clear");
                printf("Informe o intervalo onde deve se buscar a raiz aproximada \033[1;34m(ex[0,1] = digite 0 1)\n\n\033[0m");
                intervalos=valoresIntervalo();
                bolzano=metodoBisseccao(intervalos[0],intervalos[1]);
                //Se não houver raizes dentre esse intervalo
                if( !bolzano ){

                    intervalos=verificaIntervalos(intervalos);
                    metodoBisseccao(intervalos[0],intervalos[1]);
                }

                break;
                case 2:
                system("clear");
                printf("Informe o intervalo onde deve se buscar a raiz aproximada \033[1;34m(ex[0,1] = digite 0 1)\n\n\033[0m");
                intervalos=valoresIntervalo();
                bolzano=metodoPosicaoFalsa(intervalos[0],intervalos[1]);
                //Se não houver raizes dentre esse intervalo
                if( !bolzano ){

                    intervalos=verificaIntervalos(intervalos);
                    metodoPosicaoFalsa(intervalos[0],intervalos[1]);
                }

                break;
                case 3:
                system("clear");
                printf("Informe o intervalo onde deve se buscar a raiz aproximada \033[1;34m(ex[0,1] = digite 0 1)\n\n\033[0m");
                intervalos=valoresIntervalo();
                bolzano=metodoNewton(intervalos[0],intervalos[1]);
                //Se não houver raizes dentre esse intervalo
                if( !bolzano ){

                    intervalos=verificaIntervalos(intervalos);
                    metodoNewton(intervalos[0],intervalos[1]);
                }

                break;


                default:
                if(op!=4)
                printf("Opção incorreta!\n");
                break;
        }

    }while(op!=4);

    return 0;
}
