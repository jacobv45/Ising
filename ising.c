#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main(void){
	int Lx=10,Ly=10;
	int N=Lx*Ly;
	int s[Lx][Ly];
	int i,j,iNext,jNext,step,step2,thermstep=10000;
	double h=0,J=1;
	double kT=1.5;
	double E=0,Etemp=0,C=0,X=0;
	double S=0;
	double Esumtemp=0,Msumtemp=0,Esum=0,Msum=0,M2sum=0,E2sum=0;
	double Eavg=0,Mavg=0,E2avg=0,M2avg=0,Cavg=0,Xavg=0;
	int r1,r2;
	double test;		
	double z;
	double ztemp;

srand(time(NULL));


//initialize grid	
for(i=0;i<Lx;i++){
	for(j=0;j<Ly;j++){
		test=(double)rand() / (double)RAND_MAX;
		if(test>0.5){
			s[i][j]=-1;
		}
		else{
			s[i][j]=1;
		}
	}
}

for(h=0;h<1.1;h+=0.1){
for(kT=1.5;kT<7.1;kT+=0.25){
	Esum=0;
	Msum=0;
	E2sum=0;
	M2sum=0;	
//begin metropolis steps
for(step=0;step<(10000+thermstep);step++){
for(step2=0;step2<N;){	
	E=0;
	Etemp=0;	
	r1=rand()%Lx;
	r2=rand()%Ly;
	test=(double)rand() / (double)RAND_MAX;	
	
	for(i=0;i<Lx;i++){
		if(i==Lx-1){iNext=0;}
		else{iNext=i+1;}
		
		for(j=0;j<Ly;j++){
			if(j==Ly-1){jNext =  0;}
			else{jNext = j+1;}
			
			E+=-J*s[i][j]*(s[iNext][j]+s[i][jNext])-h*s[i][j];
//			printf("%d",s[i][j]);
		}
//		printf("\n");
	}
//	printf("\n");
	s[r1][r2]=-1*s[r1][r2];
	
	for(i=0;i<Lx;i++){
		if(i==Lx-1){iNext=0;}
		else{iNext=i+1;}
		for(j=0;j<Ly;j++){
			if(j==Ly-1){jNext=0;}
			else{jNext=j+1;}		
	
			Etemp+=-1*s[i][j]*(s[iNext][j]+s[i][jNext])-h*s[i][j];
//			printf("%d",s[i][j]);
		}
//		printf("\n");
	}
	z=exp(-1*(Etemp-E)/kT);
//	printf("%lf %lf %lf %lf\n",E,Etemp,z,test);
	
	if(Etemp<E){
		E=Etemp;
		step2+=1;
	}

	else if(z>test){
		E=Etemp;
		step2+=1;
	}
	else if(z<test){
		s[r1][r2]=-1*s[r1][r2];
		step2+=1;
	}
	else{
		printf("Error\n");
	}
}
//with new uncorrelated state we can calculate observables and averages 
	S=0;
	for(i=0;i<Lx;i++){
		for(j=0;j<Ly;j++){
			S=S+s[i][j];
		}
	}
	
	//this next step gives 5000 steps to allow for thermalization 	
	if(step>thermstep){	
//	printf("%d %lf %lf\n",step,E/N,S/N);	
		Msum+=S/N;
		Esum+=E/N;
		M2sum+=(S/N)*(S/N);
		E2sum+=(E/N)*(E/N);	
	}
}
Eavg=Esum/(step-thermstep);
Mavg=fabs(Msum/(step-thermstep));
E2avg=E2sum/(step-thermstep);
M2avg=M2sum/(step-thermstep);

Cavg=(E2avg-Eavg*Eavg)/kT/kT;
Xavg=(M2avg-Mavg*Mavg)/kT;

if(Mavg >0.1){
	Mavg=1;
}
else if(Mavg <0.1){
	Mavg=0;
}
//printf("%lf %lf %lf %lf %lf %lf %lf\n",kT,Eavg,Mavg,E2avg,M2avg,Cavg,Xavg);
printf("%lf %lf %lf\n",kT,h,Mavg);
}
}
return 0;
}


