#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#define nmax 5000


double f(double s,double e,double m,double r,int comp,double beta,double alpha,double deltas,double omegas,double epsilon,double deltae,double omegae,double gamma,double deltam,double omegam,double deltar,double omegar);
void RK(double h,int n,double s[nmax],double e[nmax],double m[nmax],double r[nmax],double beta,double alpha,double deltas,double omegas,double epsilon,double deltae,double omegae,double gamma,double deltam,double omegam,double deltar,double omegar);
void scelta_dati(double &tin,double &T,double &s,double &e,double &m,double &r,int &n);
void scelta_coefficienti(double &beta,double &alpha,double &deltas,double &omegas,double &epsilon,double &deltae,double &omegae,double &gamma,double &deltam,double &omegam,double &deltar,double &omegar);




int main(){
	int n,ripeti=0;
	double tin,T,s[nmax],e[nmax],m[nmax],r[nmax],beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar;

	scelta_dati(tin,T,s[0],e[0],m[0],r[0],n);
	scelta_coefficienti(beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
	
	
	FILE*file1;
	file1=fopen("dati.txt","w+");
	RK((T-tin)/n,n,s,e,m,r,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
	for(int i=0;i<=n;i++){
		fprintf(file1," %lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",tin+i*(T-tin)/n,s[i],e[i],m[i],r[i],s[i]+e[i]+m[i]+r[i]);
	}
	fclose(file1);
	FILE*file2;
	file2=fopen("grafico.txt","w+");
	fprintf(file2,"\n set samples %d;",n);
 	fprintf(file2,"\n set xrange[%lf:%lf];",tin,T);
 	fprintf(file2,"\n set size 0.8,1;");
 	fprintf(file2,"\n set xlabel 'Tempo';");
 	fprintf(file2,"\n set ylabel 'Densità';");
	fprintf(file2,"\nplot 'dati.txt' u 1:2 with lines lw 3 lc rgb 'green'  title 'Suscettibili', 'dati.txt' u 1:3 with lines lw 3 lc rgb 'violet'  title 'Esposti', 'dati.txt' u 1:4 with lines lw 3 lc rgb 'red'  title 'Infetti', 'dati.txt' u 1:5 with lines lw 3 lc rgb 'blue'  title 'Recuperati', 'dati.txt' u 1:6 with lines lw 3 lc rgb 'black'  title 'Popolazione totale'");
	fclose(file2);
	system("wgnuplot --persist grafico.txt");
}





double f(double s,double e,double m,double r,int comp,double beta,double alpha,double deltas,double omegas,double epsilon,double deltae,double omegae,double gamma,double deltam,double omegam,double deltar,double omegar){
	switch(comp){
		case 1:return(-beta*s*m+alpha*r+deltas*(s+e+m+r)-omegas*s); break;
		case 2:return(beta*s*m-epsilon*e+deltae*(s+e+m+r)-omegae*e); break;
		case 3:return(epsilon*e-gamma*m+deltam*(s+e+m+r)-omegam*m); break;
		case 4:return(gamma*m-alpha*r+deltar*(s+e+m+r)-omegar*r); break;
	}
}



void RK(double h,int n,double s[nmax],double e[nmax],double m[nmax],double r[nmax],double beta,double alpha,double deltas,double omegas,double epsilon,double deltae,double omegae,double gamma,double deltam,double omegam,double deltar,double omegar){
	double k[16];
	for(int i=0;i<n;i++){
		k[0]=f(s[i],e[i],m[i],r[i],1,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		k[1]=f(s[i],e[i],m[i],r[i],2,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		k[2]=f(s[i],e[i],m[i],r[i],3,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		k[3]=f(s[i],e[i],m[i],r[i],4,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		
		
		k[4]=f(s[i]+h/2*k[0],e[i]+h/2*k[1],m[i]+h/2*k[2],r[i]+h/2*k[3],1,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		k[5]=f(s[i]+h/2*k[0],e[i]+h/2*k[1],m[i]+h/2*k[2],r[i]+h/2*k[3],2,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		k[6]=f(s[i]+h/2*k[0],e[i]+h/2*k[1],m[i]+h/2*k[2],r[i]+h/2*k[3],3,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		k[7]=f(s[i]+h/2*k[0],e[i]+h/2*k[1],m[i]+h/2*k[2],r[i]+h/2*k[3],4,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		
		
		k[8]=f(s[i]+h/2*k[4],e[i]+h/2*k[5],m[i]+h/2*k[6],r[i]+h/2*k[7],1,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		k[9]=f(s[i]+h/2*k[4],e[i]+h/2*k[5],m[i]+h/2*k[6],r[i]+h/2*k[7],2,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		k[10]=f(s[i]+h/2*k[4],e[i]+h/2*k[5],m[i]+h/2*k[6],r[i]+h/2*k[7],3,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		k[11]=f(s[i]+h/2*k[4],e[i]+h/2*k[5],m[i]+h/2*k[6],r[i]+h/2*k[7],4,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		
		
		k[12]=f(s[i]+h*k[8],e[i]+h*k[9],m[i]+h*k[10],r[i]+h*k[11],1,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		k[13]=f(s[i]+h*k[8],e[i]+h*k[9],m[i]+h*k[10],r[i]+h*k[11],2,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		k[14]=f(s[i]+h*k[8],e[i]+h*k[9],m[i]+h*k[10],r[i]+h*k[11],3,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		k[15]=f(s[i]+h*k[8],e[i]+h*k[9],m[i]+h*k[11],r[i]+h*k[11],4,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		
		
		s[i+1]=s[i]+h*(k[0]+2*k[4]+2*k[8]+k[12])/6.000000;
		if(s[i+1]<0){
			s[i+1]=0;
		}
		e[i+1]=e[i]+h*(k[1]+2*k[5]+2*k[9]+k[13])/6.000000;
		if(e[i+1]<0){
			e[i+1]=0;
		}
		m[i+1]=m[i]+h*(k[2]+2*k[6]+2*k[10]+k[14])/6.000000;
		if(m[i+1]<0){
			m[i+1]=0;
		}
		r[i+1]=r[i]+h*(k[3]+2*k[7]+2*k[11]+k[15])/6.000000;
		if(r[i+1]<0){
			r[i+1]=0;
		}
	}
	return;
}





void scelta_dati(double &tin,double &T,double &s,double &e,double &m,double &r,int &n){
	printf("\nScrivere gli estremi t0 e T dell' intervallo di tempo in cui calcolare la traiettoria\n");
	scanf("%lf\t%lf",&tin,&T);
	while(tin<0){
		printf("\nt0 deve essere non negativo.\n");
		printf("\nInserire nuovamente t0\n");
		scanf("%lf",&tin);
	}
	while(tin>=T){
		printf("\nIl secondo estremo deve essere maggiore del primo.");
		printf("\nInserire nuovamente T maggiore di t0\n");
		scanf("%lf",&T);
	}
	printf("\nInserire le componenti del dato iniziale s(t0)  e(t0)  m(t0)  r(t0)\n");
	scanf("%lf\t%lf\%lf\t%lf",&s,&e,&m,&r);	
	if(s<0){
		s=0;
		printf("\n!!!ATTENZIONE: la componente s deve essere non negativa, pertanto e' stata posta a 0!!!\n");
	}	
	if(e<0){
		e=0;
		printf("\n!!!ATTENZIONE: la componente e deve essere non negativa, pertanto e' stata posta a 0!!!\n");
	}	
	if(m<0){
		m=0;
		printf("\n!!!ATTENZIONE: la componente m deve essere non negativa, pertanto e' stata posta a 0!!!\n");
	}	
	if(r<0){
		r=0;
		printf("\n!!!ATTENZIONE: la componente r deve essere non negativa, pertanto e' stata posta a 0!!!\n");
	}
	printf("\nInserire il numero di passi\n");
	scanf("%d",&n);
	while(n<=0){
		printf("\nIl numero di passi deve essere un numero positivo.\n");
		printf("\nInserire il numero di passi\n");
		scanf("%d",&n);
	}
	printf("\nIl tempo che intercorre tra un passo e il successivo e' h=%lf\n",(T-tin)/n);
	if((T-tin)/n>=1){
		printf("\n\n!!!ATTENZIONE: h RISULTA SUPERIORE A 1, PERTANTO LA TRAIETTORIA CALCOLATA POTREBBE ESSERE POCO AFFIDABILE!!!\n\n");
	}
	return;
}


void scelta_coefficienti(double &beta,double &alpha,double &deltas,double &omegas,double &epsilon,double &deltae,double &omegae,double &gamma,double &deltam,double &omegam,double &deltar,double &omegar){
	printf("\nInserire il valore del coefficiente beta:\n");
	scanf("%lf",&beta);
	while(beta<0){
		printf("\nTutti i coefficienti devono essere non negativi.\n");
		printf("\nInserire nuovamente beta\n");	
		scanf("%lf",&beta);	
	}
	printf("\nInserire il valore del coefficiente alpha:\n");
	scanf("%lf",&alpha);
	while(alpha<0){
		printf("\nTutti i coefficienti devono essere non negativi.\n");
		printf("\nInserire nuovamente alpha\n");	
		scanf("%lf",&alpha);	
	}
	printf("\nInserire il valore del coefficiente deltas:\n");
	scanf("%lf",&deltas);
	while(deltas<0){
		printf("\nTutti i coefficienti devono essere non negativi.\n");
		printf("\nInserire nuovamente deltas\n");	
		scanf("%lf",&deltas);	
	}
	printf("\nInserire il valore del coefficiente omegas:\n");
	scanf("%lf",&omegas);
	while(omegas<0){
		printf("\nTutti i coefficienti devono essere non negativi.\n");
		printf("\nInserire nuovamente omegas\n");	
		scanf("%lf",&omegas);	
	}
	printf("\nInserire il valore del coefficiente epsilon:\n");
	scanf("%lf",&epsilon);
	while(epsilon<0){
		printf("\nTutti i coefficienti devono essere non negativi.\n");
		printf("\nInserire nuovamente epsilon\n");	
		scanf("%lf",&epsilon);	
	}
	printf("\nInserire il valore del coefficiente deltae:\n");
	scanf("%lf",&deltae);
	while(deltae<0){
		printf("\nTutti i coefficienti devono essere non negativi.\n");
		printf("\nInserire nuovamente deltae\n");	
		scanf("%lf",&deltae);	
	}
	printf("\nInserire il valore del coefficiente omegae:\n");
	scanf("%lf",&omegae);
	while(omegae<0){
		printf("\nTutti i coefficienti devono essere non negativi.\n");
		printf("\nInserire nuovamente omegae\n");	
		scanf("%lf",&omegae);	
	}
	printf("\nInserire il valore del coefficiente gamma:\n");
	scanf("%lf",&gamma);
	while(gamma<0){
		printf("\nTutti i coefficienti devono essere non negativi.\n");
		printf("\nInserire nuovamente gamma\n");	
		scanf("%lf",&gamma);	
	}
	printf("\nInserire il valore del coefficiente deltam:\n");
	scanf("%lf",&deltam);
	while(deltam<0){
		printf("\nTutti i coefficienti devono essere non negativi.\n");
		printf("\nInserire nuovamente deltam\n");	
		scanf("%lf",&deltam);	
	}
	printf("\nInserire il valore del coefficiente omegam:\n");
	scanf("%lf",&omegam);
	while(omegam<0){
		printf("\nTutti i coefficienti devono essere non negativi.\n");
		printf("\nInserire nuovamente omegam\n");	
		scanf("%lf",&omegam);	
	}
	printf("\nInserire il valore del coefficiente deltar:\n");
	scanf("%lf",&deltar);
	while(deltar<0){
		printf("\nTutti i coefficienti devono essere non negativi.\n");
		printf("\nInserire nuovamente deltar\n");	
		scanf("%lf",&deltar);	
	}
	printf("\nInserire il valore del coefficiente omegar:\n");
	scanf("%lf",&omegar);
	while(omegar<0){
		printf("\nTutti i coefficienti devono essere non negativi.\n");
		printf("\nInserire nuovamente omegar\n");	
		scanf("%lf",&omegar);	
	}
}
