#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double gam = 1.4;

typedef struct Vars{
  double mass;
  double velocity;
  double energy;
  double press;
} Vars;

typedef struct Flux{
  double rhov;
  double mom;
  double energy;
} Flux;

double Et(double rho, double vel, double P){
  return(.5*rho*pow(vel/rho,2) + P/(gam-1.));
}

void init_sys(int N, Vars * U){
  int i;
  for(i=0;i<N;++i){
    if(i < N/2){
      U[i].mass = 1.0;
      U[i].velocity = 0.0;
      U[i].press = 1.0;
      //printf("yes %d\n",i);
    }
    else{
      U[i].mass = 0.1;
      U[i].velocity = 0.0;
      U[i].press = 0.125;
      U[i].energy = Et(U[i].mass, U[i].velocity, U[i].press);
    }
  }
    for(i=0;i<N;++i) U[i].energy = Et(U[i].mass, U[i].velocity, U[i].press);
}

double lambda(Vars U, double pm){
  double c = sqrt(gam*U.press/U.mass);
  //if(pm<0) printf("%f %f %f\n", c, U.press, U.mass);
  return(U.velocity/U.mass + pm*c);
}
double MAX(double a, double b, double c){
  double max = a;
  if(max < b) max = b;
  if(max < c) max = c;
  return(max);
}
double alpha(Vars UL, Vars UR, double pm){
  double alph = MAX(0., pm*lambda(UR, pm), pm*lambda(UL,pm));
  //printf("%f %f %f %f\n", UL.press, UL.mass, UR.press, UR.mass);
  return(alph);
}

void make_F(int N, Vars * U, Flux * F){
  int i;
  for(i=0;i<N;++i){
    F[i].rhov = U[i].velocity;
    F[i].mom = pow(U[i].velocity,2)/U[i].mass + U[i].press;
    F[i].energy = (U[i].energy + U[i].press)*U[i].velocity/U[i].mass;
    //printf("Pressure %d  %f\n", i, U[i].press);
  }
}

double HLL(int N, Vars * U, Flux * F, Flux * F_HLL){
  int i; double maxalph = 0.;
  for(i=0;i<N-1;++i){
    double alphap = alpha(U[i], U[i+1], 1.);
    double alpham = alpha(U[i], U[i+1], -1.);
    F_HLL[i].rhov = (alphap*F[i].rhov + alpham*F[i+1].rhov - alphap*alpham*(U[i+1].mass - U[i].mass))/(alphap + alpham);
    F_HLL[i].mom = (alphap*F[i].mom + alpham*F[i+1].mom - alphap*alpham*(U[i+1].velocity - U[i].velocity))/(alphap + alpham);
    F_HLL[i].energy = (alphap*F[i].energy + alpham*F[i+1].energy - alphap*alpham*(U[i+1].energy - U[i].energy))/(alphap + alpham);
    maxalph = MAX(maxalph, alphap, alpham);
    //printf("%d %f  %f\n", i, alphap, alpham);
  }
  return(maxalph);
}

double advance_system(int N, Vars * U, Flux * F, Flux * F_HLL, double dx){
  make_F(N, U, F);//get the fluxes at all points
  double maxalpha = HLL(N, U, F, F_HLL);//calculate the flux at the interfaces
  double dt = 0.5*dx/maxalpha;
  Vars U_n[N]; int i;
  double coef = dt/dx;
  for(i=1;i<N-1;++i){
    U_n[i].mass = U[i].mass - coef*(F_HLL[i].rhov - F_HLL[i-1].rhov);
    U_n[i].velocity = U[i].velocity - coef*(F_HLL[i].mom - F_HLL[i-1].mom);
    U_n[i].energy = U[i].energy - coef*(F_HLL[i].energy - F_HLL[i-1].energy);
    U_n[i].press = (gam-1.)*(U_n[i].energy - .5*pow(U_n[i].velocity,2)/U_n[i].mass);
  }
  for(i=1;i<N-1;++i){
    U[i] = U_n[i];/*
    U[i].mass = U_n[i].mass;
    U[i].velocity = U_n[i].velocity;
    U[i].energy = U_n[i].energy;
    U[i].press = U_n[i].press;*/
  }
  return(dt);
}

void Print_Vars(int N, Vars * U){
  int i; for(i=0; i<N; ++i){
    printf("Cons %f  %f  %f  %f\n", U[i].mass, U[i].velocity, U[i].energy, U[i].press);
  }
}
void Print_Flux(int N, Flux * F){
  int i;for(i=0;i<N;++i) printf("Flux %f  %f  %f\n", F[i].rhov, F[i].mom, F[i].energy);
}

void Write_Cons(int N, Vars * U, double dx, FILE * fid){
  int i; for(i=0;i<N;++i){
    fprintf(fid, "%e %e  %e  %e  %e\n", -0.5 + i*dx, U[i].mass, U[i].velocity/U[i].mass, U[i].energy, U[i].press);
  }
}
int main(void){
  FILE * fid, * finit;
  fid = fopen("data.dat", "w");
  finit = fopen("init_data.dat", "w");
  int N = 200; Vars U[N]; Flux F[N]; Flux F_HLL[N-1];
  double T, t, dt, dx; T = .14; t = 0; dx = 1./(double)N;
  init_sys(N, U); Write_Cons(N, U, dx, finit);
  while(t<T){
    dt = advance_system(N, U, F, F_HLL, dx);
    t += dt;
  }
  Write_Cons(N, U, dx, fid);
}
