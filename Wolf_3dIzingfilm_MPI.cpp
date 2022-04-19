#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "mersenne.cpp"
#include "mpi.h"
long N_flip,N_prog,count=0,maxi,istart,isum,E,E_lay[20],M,M_lay[20];
long L,l; /* Ising model parameters */
double T,T_start=3.1,T_end=3.3,T_step=0.001,Phi,J;
double ecum=0,e2cum=0,mcum=0,m2cum=0,m4cum=0,Phicum=0;
double elaycum[20],elay2cum[20],mlaycum[20],mlay2cum[20],mlay4cum[20];
double padd; /* 1 - exp(-2*J/T) */
double t1,t2;
long spin[128][128][20]; /* Lattice of spins */
int NumOfThread,ProcRank;
#define N (L*L*l) /* Total number of spins */
const double PI=3.14159265;
FILE *fp,*fp1,*fp2[20],*config;
/****************************************************/
void period(long i,long j,long k)
{  long left,right,down,up,zdown,zup;
   if(i==0) left=spin[L-1][j][k];
   else left=spin[i-1][j][k];
   if(i==L-1) right=spin[0][j][k];
   else right=spin[i+1][j][k];
   if(j==0) down=spin[i][L-1][k];
   else down=spin[i][j-1][k];
   if(j==L-1) up=spin[i][0][k];
   else up=spin[i][j+1][k];
   if(k==0) zdown=0;
   else zdown=spin[i][j][k-1];
   if(k==L-1) zup=0;
   else zup=spin[i][j][k+1];
   isum=left+right+down+up+zdown+zup;
}
/****************************************************/
void FM_start(long m,long p)
{ long i,j,k;
  for(i=0;i<m;++i)
    for(j=0;j<m;++j)
      for(k=0;k<p;++k) spin[i][j][k]=1;
}
/*****************************************************/
double multi(double a,double b)
{
          return a*b;
}
/****************************************************/
void measures()
{
    long p,m,v;
    double x_cos=0,x_sin=0,y_cos=0,y_sin=0,z_cos=0,z_sin=0;
    E=0;M=0;

    for(v=0;v<l;++v)
    {  E_lay[v]=0;
       M_lay[v]=0;
    }


    for(p=0;p<L;++p)
      for(m=0;m<L;++m)
        for(v=0;v<l;++v)
        {  period(p,m,v);
           E+=-0.5*spin[p][m][v]*(J*isum);
           E_lay[v]+=-0.5*spin[p][m][v]*(J*isum);

           M+=spin[p][m][v];
           M_lay[v]+=spin[p][m][v];

           x_cos+=spin[p][m][v]*cos(2*PI*(p+1)/L);
           x_sin+=spin[p][m][v]*sin(2*PI*(p+1)/L);
           y_cos+=spin[p][m][v]*cos(2*PI*(m+1)/L);
           y_sin+=spin[p][m][v]*sin(2*PI*(m+1)/L);
           z_cos+=spin[p][m][v]*cos(2*PI*(v+1)/l);
           z_sin+=spin[p][m][v]*sin(2*PI*(v+1)/l);
        }
    Phi=(x_cos*x_cos+x_sin*x_sin+y_cos*y_cos+
              y_sin*y_sin+z_cos*z_cos+z_sin*z_sin)/3;
}
/****************************************************/
void data()              /*накапливание данных после каждого шага Монте-Карло на спин*/
{  long i;

   ecum+=E;
   e2cum+=pow(E,2);

   mcum+=fabs(M);
   m2cum+=pow(M,2);
   m4cum+=pow(M,4);

   for(i=0;i<l;++i)
   {
      elaycum[i]+=E_lay[i];
      elay2cum[i]+=pow(E_lay[i],2);

      mlaycum[i]+=fabs(M_lay[i]);
      mlay2cum[i]+=pow(M_lay[i],2);
      mlay4cum[i]+=pow(M_lay[i],4);
   }

   Phicum+=Phi;
}
/****************************************************/
void output()
{ double znorm,eav,e2av,mav,m2av,m4av,Phiav,ksi2_xy,ksi2_z,U;
  double elayav[20],elay2av[20],mlayav[20],mlay2av[20],mlay4av[20];
  znorm=1/multi(N,count);
  eav=ecum*znorm;
  e2av=e2cum*znorm;
  mav=mcum*znorm;
  m2av=m2cum*znorm;
  m4av=m4cum*znorm;
  Phiav=Phicum*znorm;
  U=0.5*(3-m4av/(N*m2av*m2av));
  ksi2_xy=sqrt((0.25/pow(sin(PI/L),2))*(m2av/Phiav-1));
  ksi2_z=sqrt((0.25/pow(sin(PI/l),2))*(m2av/Phiav-1));

  fprintf(fp1,"%f	%f	%f	%f	%f	%f	%f",T,eav,e2av,mav,m2av,m4av,U);
  fprintf(fp1,"	%f	%f	%f",Phiav,ksi2_xy,ksi2_z);
  fprintf(fp1,"\n");
  for(long f=0;f<l;++f)
  {  elayav[f]=elaycum[f]*znorm*l;
     elay2av[f]=elay2cum[f]*znorm*l;
     mlayav[f]=mlaycum[f]*znorm*l;
     mlay2av[f]=mlay2cum[f]*znorm*l;
     mlay4av[f]=mlay4cum[f]*znorm*l;
     U=0.5*(3-mlay4av[f]/(N*mlay2av[f]*mlay2av[f]/l));

     fprintf(fp2[f],"	%f	%f",elayav[f],elay2av[f]);
     fprintf(fp2[f],"	%f	%f	%f",mlayav[f],mlay2av[f],mlay4av[f]);
     fprintf(fp2[f],"	%f",U);
     fprintf(fp2[f],"\n");
  }
}
/****************************************************/
void output_info()
{
    fprintf(fp,"\nL=%ld  l=%ld",L,l);
    fprintf(fp,"\nN_flip=%ld",N_flip);
    fprintf(fp,"\nN_prog=%ld",N_prog);
    fprintf(fp,"\nчисло шагов Монте-Карло на спин=%ld",maxi);
    fprintf(fp,"\nчисло начальных отбрасываемых конфигураций=%ld",istart);
    fprintf(fp,"\nT_start=%f	T_end=%f	T_step=%f",T_start,T_end,T_step);
    fprintf(fp,"\nt=%f",(t2-t1));
}
/*********************************************************/
void Reading()
{
    config=fopen("config.cfg","r");
    fscanf(config,"%ld %ld %lf %lf %lf %ld %ld %lf	%ld %ld"
         ,&L,&l,&T_start,&T_end,&T_step,&maxi,&istart,&J,&N_flip,&N_prog);
    fclose(config);
}
/**********************BEGIN****************************/

int main(int argc,char *argv[])
{ long m,f,p;
  char s[50];
  long sp; /* Stack pointer */
  long i,j,k,oldspin,newspin,current,nn;
  t1=clock();
  Reading();

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &NumOfThread);
  MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);


  for(p=0;p<N_prog;++p)
  {  sprintf(s,"Ln4=%ld ProcRank=%d N_prog=%ld.dat",L,ProcRank,p+1);
     fp=fopen(s,"a");
     sprintf(s,"stat_Ln4=%ld ProcRank=%d N_prog=%ld.dat",L,ProcRank,p+1);
     fp1=fopen(s,"a");
     for(f=0;f<l;++f)
     {  sprintf(s,"laystat_Ln4=%ld l=%ld ProcRank=%d N_prog=%ld.dat",L,f+1,ProcRank,p+1);
        fp2[f]=fopen(s,"a");
     }
     FM_start(L,l);
     CRandomMersenne rg(time(NULL)*rand()+ProcRank);

     for(T=T_start;T<=T_end;T+=T_step)
     { mcum=m2cum=m4cum=0;
       ecum=e2cum=0;
       Phicum=0;
       count=0;
       padd=1-exp(-2*J/T);
       for(f=0;f<l;++f)
       {  elaycum[f]=elay2cum[f]=0;
          mlaycum[f]=mlay2cum[f]=mlay4cum[f]=0;
       }
       for(m=0;m<maxi;++m)
       {  for(f=0;f<N_flip;++f)
          {  long stack[N]; /* LIFO stack */

              /* Put a random seed spin site onto a stack */
             i = rg.IRandom(0,L-1); j = rg.IRandom(0,L-1);k = rg.IRandom(0,l-1);
             stack[0] = i*L + j+ k*L*(L+l);
             sp =1;
            /* Flip the seed and remember the old & new spins */
            oldspin = spin[i][j][k]; newspin = -spin[i][j][k];
            spin[i][j][k] = newspin;

            while (sp) {
            /* Pop a site off the stack */
            current = stack[--sp];
            i = (current/L)%(L+l);
            j = current%L;
            k = current/(L*(L+l));
            /* Check the neighbors */
            if ((nn=i+1) >= L) nn -= L; /* South neighbor */
            if (spin[nn][j][k] == oldspin)
            if (rg.Random() < padd) {
            stack[sp++] = nn*L + j+k*L*(L+l);
            spin[nn][j][k] = newspin;
            }

            if ((nn=i-1) < 0) nn += L; /* North neighbor */
            if (spin[nn][j][k] == oldspin)
            if (rg.Random() < padd) {
            stack[sp++] = nn*L + j+k*L*(L+l);
            spin[nn][j][k] = newspin;
            }

            if ((nn=j+1) >= L) nn -= L; /* East neighbor */
            if (spin[i][nn][k] == oldspin)
            if (rg.Random() < padd) {
            stack[sp++] = i*L + nn+k*L*(L+l);
            spin[i][nn][k] = newspin;
            }

            if ((nn=j-1) < 0) nn += L; /* West neighbor */
            if (spin[i][nn][k] == oldspin)
            if (rg.Random() < padd) {
            stack[sp++] = i*L + nn+k*L*(L+l);
            spin[i][nn][k] = newspin;
            }

            if ((nn=k+1) >= l) ; /* zDown neighbor */
            else
            {
                if (spin[i][j][nn] == oldspin)
                if (rg.Random() < padd) {
                stack[sp++] = i*L + j+nn*L*(L+l);
                spin[i][j][nn] = newspin;
                }
            }

            if ((nn=k-1) < 0) ; /* zUp neighbor */
            else
            {
                if (spin[i][j][nn] == oldspin)
                if (rg.Random() < padd) {
                stack[sp++] = i*L + j+nn*L*(L+l);
                spin[i][j][nn] = newspin;
                }
            }

            } /* End while stack is not empty */

         }



          if(m>=istart)
          {  measures();
             data();
             ++count;
          }

       }
       output();
     }
     t2=clock();
     output_info();
     fclose(fp);
     fclose(fp1);
     for(f=0;f<l;++f) fclose(fp2[f]);
  }
  MPI_Finalize();
}
