#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#define mmax 100


int modello(double &a,double &b,double &R,double &S,double &T,double &R0,double &xe,int &c);
double f(double y,double z,double x,int k,int specie,double a,double b,int anno,double R,double S,double T,double R0,double xe,int c);
void scrivi_funzione_su_file(double pop[3],int n,double a,double b,int k,double R,double S,double T,double R0,double xe,int c);
void plotta(int n);


int main(){
	double a,b,pop[3],S,R,T,R0,xe;
	int scelta,n,k,c;
	printf("MODELLO PREDA-PREDATORE CON STRUTTURA D'ETA' SUL PREDATORE\n\n\n");
	printf("Inserire il numero di periodi su cui eseguire lo studio\n");
	scanf("%d",&n);
	printf("Inserire il numero di predatori maturi all'inizio dello studio (densita')\n");
	scanf("%lf",&pop[0]);
	if(pop[0]<0){
		pop[0]=0;
		printf("\nATTENZIONE. Valore negativo settato a 0.\n");
	}
	printf("Inserire il numero di predatori immaturi all'inizio dello studio (densita')\n");
	scanf("%lf",&pop[1]);
	if(pop[1]<0){
		pop[1]=0;
		printf("\nATTENZIONE. Valore negativo settato a 0.\n");
	}
	printf("Inserire il numero di prede all'inizio dello studio (densita')\n");
	scanf("%lf",&pop[2]);
	if(pop[2]<0){
		pop[2]=0;
		printf("\nATTENZIONE. Valore negativo settato a 0.\n");
	}
	k=modello(a,b,R,S,T,R0,xe,c);
	scrivi_funzione_su_file(pop,n,a,b,k,R,S,T,R0,xe,c);
	plotta(n);
}



int modello(double &a,double &b,double &R,double &S,double &T,double &R0,double &xe,int &c){
	int scelta;
	printf("\n"
	"(y=predatori maturi\tz=predatori immaturi\tx=prede)\n\n"
	"Modello di Maynard Smith-Slatkin:\n"
	"\tinverno:\n"
	"\ty(t+1)=y(t)*[1-exp(-a*x(t))]\n"
	"\tz(t+1)=z(t)*[1-exp(-b*x(t))]\n"
	"\tx(t+1)=x(t)-y(t+1)-z(t+1)\n\n"
	"\testate:\n"
	"\ty(t+1)=S*y(t)+z(t)\n"
	"\tx(t+1)=R*x(t)\n"
	"\tz(t+1)=T*y(t+1)\n\n");
		printf("\nInserire l'area di esplorazione dei predatori maturi:a=");
		scanf("%lf",&a);
		printf("\nInserire l'area di esplorazione dei predatori immaturi:b=");
		scanf("%lf",&b);
		printf("\nInserire l'indice di riproduzione dei predatori maturi:T=");
		scanf("%lf",&T);
		printf("\nInserire l'indice di sopravvivenza dei predatori maturi:S=");
		scanf("%lf",&S);
	printf("\nSelezionare R:\n"
	"\t1)R costante\n"
	"\t2)R=R0/(1+(R0-1)*pow(x/xe,c))\n");
	scanf("%d",&scelta); 
	while((scelta<1) || (scelta>2))
    {
      printf("\nIl valore inserito non e' nella lista.");
      printf("\nSelezionare R:\n");
      scanf("%d",&scelta);
    }
	switch(scelta){
	case 1:	printf("\nInserire l'indice di riproduzione delle prede:R=");
		scanf("%lf",&R);
		printf("\n");
		break;
	case 2:
			printf("\nInserire R0=");
			scanf("%lf",&R0);
			printf("\nInserire xe=");
			scanf("%lf",&xe);
			printf("\nInserire c=");
			scanf("%d",&c);
			printf("\n");
			break;
    default: return(0); break;	
	}
	return(scelta);
}




double f(double y,double z,double x,int k,int specie,double a,double b,int anno,double R,double S,double T,double R0,double xe,int c){
	if(anno%21==0){
		switch(specie){
			case 1:	return(S*y+z); break;
			case 2:	return(T*y); break;
			case 3:
				switch(k){
					case 1:	return(R*x); break;
					case 2: R=R0/(1+(R0-1)*pow(x/xe,c)); return(R*x);break;
					default: return 0; break;
				}
			default: return 0; break;
		}
	}else{
		switch(specie){
			case 1:	return(y*(1-exp(-a*x))); break;
			case 2:	return(z*(1-exp(-b*x))); break;
			case 3:	return(x-y-z); break;
			default: return 0; break;
		}
	}
}


void scrivi_funzione_su_file(double pop[3],int n,double a,double b,int k,double R,double S,double T,double R0,double xe,int c){
	FILE *file;
	file=fopen("dati_funzione.txt","w+");
	fprintf(file," %lf \t %lf \t %lf \t %lf \n",0,pop[0],pop[1],pop[2]/20);
		for(int i=1;i<=n;i++){
			pop[0]=f(pop[0],pop[1],pop[2],k,1,a,b,i,R,S,T,R0,xe,c);
			pop[1]=f(pop[0],pop[1],pop[2],k,2,a,b,i,R,S,T,R0,xe,c);
			pop[2]=f(pop[0],pop[1],pop[2],k,3,a,b,i,R,S,T,R0,xe,c);
			if(pop[2]<0){
				pop[2]=0;
			}
			if(n>=200){
				if(i%21==0){
					fprintf(file," %d \t %lf \t %lf \t %lf \n",i,pop[0],pop[1],pop[2]/20);
				}
			}else{
					fprintf(file," %d \t %lf \t %lf \t %lf \n",i,pop[0],pop[1],pop[2]/20);
				}
		}
	fclose(file);
}


void plotta(int n){
	FILE*file;
 	file=fopen("comandi.txt","w+");
	fprintf(file,"\n set samples %d;",n);
 	fprintf(file,"\n set xrange[%d:%d];",0,n);
// 	fprintf(file,"\n set yrange[%d:%d];",0,230);
  //	fprintf(file,"\n set title 'Preda-predatore con età';");
	fprintf(file,"\n plot 'dati_funzione.txt' u 1:2 with linespoints lw 2 lc rgb 'red' pt 7 ps 0.8 title 'Predatori maturi', 'dati_funzione.txt' u 1:3 with linespoints lw 2 lc rgb 'blue' pt 7 ps 0.8 title 'Predatori immaturi', 'dati_funzione.txt' u 1:4 with linespoints lw 2 lc rgb 'green' pt 7 ps 0.8 title 'Prede'");
	fclose(file);
	system("wgnuplot -persist comandi.txt");
}









