#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double gam = 1.4;
double thet = 1.5;

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

double sgn(double x){
  if(x < 0) return(-1.);
  else return(1.);
}

double min(double x, double y, double z){
  double num = x;
  if(num > y) num = y;
  if(num > z) num = z;
  return(num);
}
double minmod(double x, double y, double z){
  return(.25 * fabs(sgn(x) + sgn(y))*fabs(sgn(x) + sgn(z))*min(fabs(x), fabs(y), fabs(z)));
}

Vars U_L(Vars Ui, Vars Uimo, Vars Uipo){
  Vars UL;
  UL.mass = Ui.mass + 0.5*minmod(thet*(Ui.mass - Uimo.mass), 0.5*(Uipo.mass - Uimo.mass),thet*(Uipo.mass-Ui.mass));
  UL.velocity = Ui.velocity + 0.5*minmod(thet*(Ui.velocity - Uimo.velocity), 0.5*(Uipo.velocity - Uimo.velocity),thet*(Uipo.velocity-Ui.velocity));
  UL.energy = Ui.energy + 0.5*minmod(thet*(Ui.energy - Uimo.energy), 0.5*(Uipo.energy - Uimo.energy),thet*(Uipo.energy-Ui.energy));
  UL.press = Ui.press + 0.5*minmod(thet*(Ui.press - Uimo.press), 0.5*(Uipo.press - Uimo.press),thet*(Uipo.press-Ui.press));
  return(UL);
}

Vars U_R(Vars Ui, Vars Uipo, Vars Uipt){
  Vars UR;
  UR.mass = Ui.mass - 0.5*minmod(thet*(Uipo.mass-Ui.mass),.5*(Uipt.mass-Ui.mass),thet*(Uipt.mass-Uipo.mass));
  UR.velocity = Ui.velocity - 0.5*minmod(thet*(Uipo.velocity-Ui.velocity),.5*(Uipt.velocity-Ui.velocity),thet*(Uipt.velocity-Uipo.velocity));
  UR.energy = Ui.energy - 0.5*minmod(thet*(Uipo.energy-Ui.energy),.5*(Uipt.energy-Ui.energy),thet*(Uipt.energy-Uipo.energy));
  UR.press = Ui.press - 0.5*minmod(thet*(Uipo.press-Ui.press),.5*(Uipt.press-Ui.press),thet*(Uipt.press-Uipo.press));
  return(UR);
}

void init_sys(int N, Vars * U){
  int i;
  for(i=0;i<N+4;++i){
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

Flux get_F(Vars U){
  Flux F;  
  F.rhov = U.velocity;
  F.mom = pow(U.velocity,2)/U.mass + U.press;
  F.energy = (U.energy + U.press)*U.velocity/U.mass;
  return(F);
}
Vars get_UL(int i, Vars * U){
  Vars UL = U_L(U[i], U[i-1], U[i+1]);
  return(UL);
}
Vars make_UR(int i, Vars * U){
  Vars UR = U_R(U[i], U[i+1], U[i+2]);
  return(UR);
}

double HLL(int N, Vars * U, Flux * F_HLL){
  int i; double maxalph = 0.;
  for(i=2;i<N+3;++i){
    Vars UL = U_L(U[i], U[i-1], U[i+1]); Vars UR = U_R(U[i], U[i+1], U[i+2]); //Calculates the conserved variables at the interfaces
    double alphap = alpha(UL, UR, 1.);
    double alpham = alpha(UL, UR, -1.);
    Flux FL = get_F(UL); Flux FR = get_F(UR); //Calculates the Fluxes from the left and right at the interface
    F_HLL[i-2].rhov = (alphap*FL.rhov + alpham*FR.rhov - alphap*alpham*(UR.mass - UL.mass))/(alphap + alpham);
    F_HLL[i-2].mom = (alphap*FL.mom + alpham*FR.mom - alphap*alpham*(UR.velocity - UL.velocity))/(alphap + alpham);
    F_HLL[i-2].energy = (alphap*FL.energy + alpham*FR.energy - alphap*alpham*(UR.energy - UL.energy))/(alphap + alpham);
    maxalph = MAX(maxalph, alphap, alpham);
    printf("alpha %d %f  %f %f\n", i, alphap, alpham, maxalph);
  }
  return(maxalph);
}

double advance_system(int N, Vars * U, Flux * F_HLL, double dx){
  double maxalpha = HLL(N, U, F_HLL);//calculate the flux at the interfaces
  double dt = 0.5*dx/maxalpha;//calculate the time step
  Vars U_n[N]; int i;
  double coef = dt/dx;
  for(i=2;i<N+2;++i){
    U_n[i].mass = U[i].mass - coef*(F_HLL[i-1].rhov - F_HLL[i-2].rhov);
    U_n[i].velocity = U[i].velocity - coef*(F_HLL[i-1].mom - F_HLL[i-2].mom);
    U_n[i].energy = U[i].energy - coef*(F_HLL[i-1].energy - F_HLL[i-2].energy);
    U_n[i].press = (gam-1.)*(U_n[i].energy - .5*pow(U_n[i].velocity,2)/U_n[i].mass);
  }
  for(i=2;i<N+2;++i){
    U[i] = U_n[i];
    /*U[i].mass = U_n[i].mass;
    U[i].velocity = U_n[i].velocity;
    U[i].energy = U_n[i].energy;
    U[i].press = U_n[i].press;*/
  }
  //printf("max %f\n", maxalpha);
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
  int i; for(i=0;i<N+4;++i){
    fprintf(fid, "%e %e  %e  %e  %e\n", -0.5 + i*dx, U[i].mass, U[i].velocity/U[i].mass, U[i].energy, U[i].press);
  }
}
int main(void){
  FILE * fid, * finit;
  fid = fopen("data.dat", "w");
  finit = fopen("init_data.dat", "w");
  int N = 50; Vars U[N+4]; Flux F_HLL[N+1];
  double T, t, dt, dx; T = .14; t = 0; dx = 1./(double)N;
  init_sys(N, U); Write_Cons(N, U, dx, finit);
  //int j;for(j=0;j<2;++j){
  while(t<T){
    dt = advance_system(N, U, F_HLL, dx);
    t += dt; printf("dt %f\n", dt);
    break;
  }
  Print_Flux(N+1, F_HLL);
  Write_Cons(N, U, dx, fid);
}
