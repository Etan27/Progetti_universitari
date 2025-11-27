#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#define nmax 5000


double f(double s,double e,double m,double r,int comp,double beta,double alpha,double deltas,double omegas,double epsilon,double deltae,double omegae,double gamma,double deltam,double omegam,double deltar,double omegar);
void RK(double h,int n,double s[nmax],double e[nmax],double m[nmax],double r[nmax],double beta,double alpha,double deltas,double omegas,double epsilon,double deltae,double omegae,double gamma,double deltam,double omegam,double deltar,double omegar);
void scelta_dati(double &tin,double &T,double &s,double &e,double &m,double &r,int &n);
void scelta_coefficienti(double &beta,double &alpha,double &deltas,double &omegas,double &epsilon,double &deltae,double &omegae,double &gamma,double &deltam,double &omegam,double &deltar,double &omegar);




// m sono gli infetti (per non usare la lettera i), sta per malati




// e[0] è diverso da 0, non per forza m[0]



int main(){
	int n,ripeti=0;
	double tin,T,s[nmax],e[nmax],m[nmax],r[nmax],beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar;
	do{
	//	scelta_dati(tin,T,s[0],e[0],m[0],r[0],n);
	//	scelta_coefficienti(beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
	
//	 tin=0;T=200;s[0]=9999;e[0]=1;m[0]=0;r[0]=0;n=1200;
//		
//		beta=2.8; alpha=0; deltas=0;omegas=0; epsilon=0.2; deltae=0; omegae=0; gamma=0.8; deltam=0; omegam=0; deltar=0; omegar=0;     	    s[0]=0.9999;e[0]=0.0001;m[0]=0;r[0]=0; T=90;tin=0;n=3251;													//1) Modello SEIR epidemico: https://www.math.uci.edu/~chenlong/CAMtips/Coronavirus/MathModelCoV19.html
//	
//		beta=2.8; alpha=0; deltas=0.03;omegas=0.03; epsilon=0.2; deltae=0; omegae=0.03; gamma=0.8; deltam=0; omegam=0.03; deltar=0; omegar=0.03;     	    s[0]=0.9999;e[0]=0.00001;m[0]=0;r[0]=0; T=200;tin=0;n=3251;									//2) è il 1) con l'introduzione di natalità in s e mortalità, ma con gli stessi coefficienti
//
		beta=2.8; alpha=0; deltas=0.02;omegas=0.015; epsilon=0.2; deltae=0; omegae=0.02; gamma=0.8; deltam=0; omegam=0.05; deltar=0; omegar=0.015;     	    s[0]=0.9999;e[0]=0.00001;m[0]=0;r[0]=0; T=200;tin=0;n=3251;									//3) è il 2) con mortalità differenti; popolazione totale in crescita
//
		beta=2.8; alpha=0; deltas=0;omegas=0; epsilon=0.2; deltae=0; omegae=0; gamma=0.8; deltam=0; omegam=0; deltar=0; omegar=0;     	    s[0]=0.7;e[0]=0.3;m[0]=0;r[0]=0; T=90;tin=0;n=3251;															//4) o 6) è il 1) con dati iniziali cambiati
//
//		beta=2.8; alpha=0; deltas=0.03;omegas=0.03; epsilon=0.2; deltae=0; omegae=0.03; gamma=0.8; deltam=0; omegam=0.03; deltar=0; omegar=0.03;     	    s[0]=0.7;e[0]=0.3;m[0]=0;r[0]=0; T=200;tin=0;n=3251; 										//5) o 7) è il 2) con dati iniziali diversi; !!!!!!!STESSO EQUILIBRIO E OSSERVAZIONI SU EQUILIBRI!!!!!!!
//	
		beta=2.8; alpha=0; deltas=0.02;omegas=0.015; epsilon=0.2; deltae=0; omegae=0.02; gamma=0.8; deltam=0; omegam=0.05; deltar=0; omegar=0.015;     	    s[0]=0.2;e[0]=0.3;m[0]=0.4;r[0]=0.1; T=200;tin=0;n=3251;										//6) o 8) è il 3) con dati iniziali diversi, forse poco interessante
		
		
//
//		beta=0.0028; alpha=0; deltas=0;omegas=0; epsilon=0.008; deltae=0; omegae=0; gamma=0.01; deltam=0; omegam=0; deltar=0; omegar=0;		tin=0;T=200;s[0]=9999;e[0]=1;m[0]=0;r[0]=0;n=3251;															//7) o 4) SEIR epid con densità superficiale
//	
//		beta=0.0004; alpha=0; deltas=0;omegas=0; epsilon=0.04; deltae=0; omegae=0; gamma=0; deltam=0; omegam=0; deltar=0; omegar=0;			tin=0;T=90;s[0]=9997;e[0]=3;m[0]=0;r[0]=0;n=3251;  															//8) o 5) sembra stupido, ma è un modo di utilizzare il modello per studiare il modello sir; può valere per l'hiv; calcolare il picco di e con formula per imax;
//	


  	beta=0.004; alpha=0.1; deltas=0.62;omegas=0.9; epsilon=0.7; deltae=0.1; omegae=0.9; gamma=0.2; deltam=0.02; omegam=1.4; deltar=0.2; omegar=0.75;		s[0]=997;e[0]=3;m[0]=0;r[0]=0; T=200;tin=0;n=3251;										//9) situazione totale e stabile; la malattia diventa endemica e si stabilizza tutto; la situazione sociale inizialmente spaventosa diventa la normalità
//	
		beta=0.0009; alpha=0.5259; deltas=0.5932;omegas=0.8899; epsilon=0.92; deltae=0.0034; omegae=1.7725; gamma=0.7885; deltam=0.0729; omegam=0.9776; deltar=0.2644; omegar=0.6941;		s[0]=997;e[0]=3;m[0]=0;r[0]=0; T=200;tin=0;n=3251;			//10) ; si può congetturare la presenza di un punto di equilibrio asintoticamente stabile
//	
//		beta=0.00019; alpha=0.04; deltas=0.045;omegas=0.13; epsilon=0.87; deltae=0.45; omegae=0.62; gamma=0.15; deltam=0.12; omegam=1.02; deltar=0.05; omegar=1.74;					s[0]=997;e[0]=3;m[0]=0;r[0]=0; T=200;tin=0;n=3251;					//11) ; interessante su 100 anni con la popolazione totale//			
//		
//		beta=0.004; alpha=0.1; deltas=0.62;omegas=0.9; epsilon=0.7; deltae=0.1; omegae=0.9; gamma=0.2; deltam=0.02; omegam=1.4; deltar=0.2; omegar=0.75;		s[0]=300;e[0]=300;m[0]=300;r[0]=1000; T=200;tin=0;n=3251;								//12) è il 9) con dati diversi
//
//		beta=0.004; alpha=0.1; deltas=0.62;omegas=0.9; epsilon=0.7; deltae=0.1; omegae=0.9; gamma=0.2; deltam=0.02; omegam=1.4; deltar=0.2; omegar=0.75;		s[0]=900;e[0]=1700;m[0]=1900;r[0]=2300; T=200;tin=0;n=3251;								//13) è il 9) con dati molto diversi
//
//		beta=0.004; alpha=0.1; deltas=0.62;omegas=0.9; epsilon=0.7; deltae=0.1; omegae=0.9; gamma=0.2; deltam=0.02; omegam=1.4; deltar=0.2; omegar=0.75;		s[0]=690;e[0]=1054;m[0]=500;r[0]=844; T=200;tin=0;n=3251;								//14) equilibrio del sistema 9)
//	
//	









//		beta=0.000769; alpha=0.732; deltas=0.731;omegas=0.027; epsilon=0.65; deltae=0.42; omegae=1.76; gamma=0.38; deltam=0.059; omegam=1.37; deltar=0.024; omegar=0.55;					s[0]=997;e[0]=3;m[0]=0;r[0]=0; T=100;tin=0;n=3251;			// ????
//	
//		beta=0.00053; alpha=0.37; deltas=0.46;omegas=1.06; epsilon=0.90; deltae=0.39; omegae=1.05; gamma=0.68; deltam=0.05; omegam=1.94; deltar=0.46; omegar=1.74;					s[0]=997;e[0]=3;m[0]=0;r[0]=0; T=200;tin=0;n=3251;				//è uguale a si può congetturare la presenza di un punto di equilibrio asintoticamente stabile
//	
	
	
	//	beta=2.8; alpha=0; deltas=0;omegas=0; epsilon=0.2; deltae=0; omegae=0; gamma=0.8; deltam=0; omegam=0; deltar=0; omegar=0;
	//	srand(time(NULL));
	//	beta=(float) rand()/RAND_MAX *7;beta=(float) rand()/RAND_MAX *7;beta=(float) rand()/RAND_MAX *7;beta=(float) rand()/RAND_MAX *7;beta=(float) rand()/RAND_MAX *0.00010; alpha=(float) rand()/RAND_MAX *1; deltas=(float) rand()/RAND_MAX *1;omegas=(float) rand()/RAND_MAX *1.00000001; epsilon=(float) rand()/RAND_MAX *0.7; deltae=(float) rand()/RAND_MAX *1; omegae=(float) rand()/RAND_MAX *1.0000001; gamma=(float) rand()/RAND_MAX *0.7; deltam=(float) rand()/RAND_MAX *1; omegam=(float) rand()/RAND_MAX *1.00000001; deltar=(float) rand()/RAND_MAX *1; omegar=(float) rand()/RAND_MAX *1.00000001;						s[0]=0.9999;e[0]=0.0001;m[0]=0;r[0]=0; T=100;tin=0;n=1200;
	//	printf("beta=%lf\talpha=%lf\tdeltas=%lf\tomegas=%lf\tepsilon=%lf\tdeltae=%lf\tomegae=%lf\tgamma=%lf\tdeltam=%lf\tomegam=%lf\tdeltar=%lf\tomegar=%lf\t",beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
	
	//	beta=0.005890; alpha=0.707877; deltas=0.120743;omegas=1.246132; epsilon=0.863247; deltae=0.3932; omegae=1.4969; gamma=0.4969; deltam=0.228; omegam=1.57; deltar=0.4422; omegar=0.71;		s[0]=997;e[0]=100;m[0]=100;r[0]=100; T=50;tin=0;n=1200;
	
	//	beta=; alpha=; deltas=;omegas=; epsilon=; deltae=; omegae=; gamma=; deltam=; omegam=; deltar=; omegar=;					s[0]=997;e[0]=3;m[0]=0;r[0]=0; T=200;tin=0;n=1200;
	
	//	beta=0.00022; alpha=0.83; deltas=0.073;omegas=0.49; epsilon=0.77; deltae=0.14; omegae=0.52; gamma=0.46; deltam=0.57; omegam=1.32; deltar=0.15; omegar=1.25;					s[0]=997;e[0]=100;m[0]=100;r[0]=100; T=200;tin=0;n=1200;	
	

	
	
		
//				beta=0.4; alpha=0.; deltas=0.005;omegas=0.; epsilon=0.2; deltae=0.; omegae=0.; gamma=0.; deltam=0.; omegam=0.005; deltar=0.; omegar=0.;		s[0]=995;e[0]=5;m[0]=0;r[0]=0; T=200;tin=0;n=2531;
		
		
		
		
	
		
		FILE*file1;
		file1=fopen("dati.txt","w+");
		RK((T-tin)/n,n,s,e,m,r,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		for(int i=0;i<=n;i++){
	//		fprintf(file1," %lf \t %lf \t %lf \t %lf \t %lf \n",tin+i*(T-tin)/n,s[i],e[i],m[i],r[i]);
			fprintf(file1," %lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",tin+i*(T-tin)/n,s[i],e[i],m[i],r[i],s[i]+e[i]+m[i]+r[i]);
		}
		fclose(file1);
		FILE*file2;
		file2=fopen("grafico.txt","w+");
		fprintf(file2,"\n set samples %d;",n);
	 	fprintf(file2,"\n set xrange[%lf:%lf];",tin,T);
	 	fprintf(file2,"\n set yrange[0:7900];");
	 	fprintf(file2,"\n set size 0.8,1;");
	 	fprintf(file2,"\n set xlabel 'Tempo';");
	 	fprintf(file2,"\n set ylabel 'Densità';");
		fprintf(file2,"\nplot 'dati.txt' u 1:2 with lines lw 3 lc rgb 'green'  title 'Suscettibili', 'dati.txt' u 1:3 with lines lw 3 lc rgb 'violet'  title 'Esposti', 'dati.txt' u 1:4 with lines lw 3 lc rgb 'red'  title 'Infetti', 'dati.txt' u 1:5 with lines lw 3 lc rgb 'blue'  title 'Recuperati'");
	//	fprintf(file2,"\nplot 'dati.txt' u 1:2 with lines lw 3 lc rgb 'green'  title 'Suscettibili', 'dati.txt' u 1:3 with lines lw 3 lc rgb 'violet'  title 'Esposti', 'dati.txt' u 1:4 with lines lw 3 lc rgb 'red'  title 'Infetti', 'dati.txt' u 1:5 with lines lw 3 lc rgb 'blue'  title 'Recuperati', 'dati.txt' u 1:6 with lines lw 3 lc rgb 'black'  title 'Popolazione totale'");
		fclose(file2);
		system("start C:\\Adesso\\gnuplot\\bin\\wgnuplot.exe --persist grafico.txt");
		printf("\n\nSi vuole far ripartire il programma?\n0)Si\t1)No\n");
		scanf("%d",&ripeti);
	}while(ripeti==0);
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
		
		//printf("\n%lf\t%lf\t%lf\t%lf",k[0],k[1],k[2],k[3]);
		
		k[4]=f(s[i]+h/2*k[0],e[i]+h/2*k[1],m[i]+h/2*k[2],r[i]+h/2*k[3],1,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		k[5]=f(s[i]+h/2*k[0],e[i]+h/2*k[1],m[i]+h/2*k[2],r[i]+h/2*k[3],2,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		k[6]=f(s[i]+h/2*k[0],e[i]+h/2*k[1],m[i]+h/2*k[2],r[i]+h/2*k[3],3,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		k[7]=f(s[i]+h/2*k[0],e[i]+h/2*k[1],m[i]+h/2*k[2],r[i]+h/2*k[3],4,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		
		//printf("\n%lf\t%lf\t%lf\t%lf",k[4],k[5],k[6],k[7]);
		
		k[8]=f(s[i]+h/2*k[4],e[i]+h/2*k[5],m[i]+h/2*k[6],r[i]+h/2*k[7],1,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		k[9]=f(s[i]+h/2*k[4],e[i]+h/2*k[5],m[i]+h/2*k[6],r[i]+h/2*k[7],2,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		k[10]=f(s[i]+h/2*k[4],e[i]+h/2*k[5],m[i]+h/2*k[6],r[i]+h/2*k[7],3,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		k[11]=f(s[i]+h/2*k[4],e[i]+h/2*k[5],m[i]+h/2*k[6],r[i]+h/2*k[7],4,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		
		//printf("\n%lf\t%lf\t%lf\t%lf",k[8],k[9],k[10],k[11]);
		
		k[12]=f(s[i]+h*k[8],e[i]+h*k[9],m[i]+h*k[10],r[i]+h*k[11],1,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		k[13]=f(s[i]+h*k[8],e[i]+h*k[9],m[i]+h*k[10],r[i]+h*k[11],2,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		k[14]=f(s[i]+h*k[8],e[i]+h*k[9],m[i]+h*k[10],r[i]+h*k[11],3,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		k[15]=f(s[i]+h*k[8],e[i]+h*k[9],m[i]+h*k[11],r[i]+h*k[11],4,beta,alpha,deltas,omegas,epsilon,deltae,omegae,gamma,deltam,omegam,deltar,omegar);
		
		//printf("\n%lf\t%lf\t%lf\t%lf",k[12],k[13],k[14],k[15]);
		
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
