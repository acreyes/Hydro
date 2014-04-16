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
  else if(x >0) return(1.);
  else return(0.);
}

double min(double x, double y, double z){
  double num = x;
  if(num > y) num = y;
  if(num > z) num = z;
  return(num);
}
double minmod(double x, double y, double z){
  double minm = .25*fabs(sgn(x) + sgn(y))*(sgn(x) + sgn(z))*min(fabs(x), fabs(y), fabs(z));
  //if(minm != 0. && minm!= -0.) printf("x y z = %f %f %f minmod = %e\n", x, y, z, minm);
  return(minm);
}

Vars U_L(Vars Ui, Vars Uimo, Vars Uipo, double dx){
  Vars UL;
  UL.mass = Ui.mass + 0.5*dx*minmod(thet*(Ui.mass - Uimo.mass), 0.5*(Uipo.mass - Uimo.mass),thet*(Uipo.mass-Ui.mass));
  UL.velocity = Ui.velocity + 0.5*dx*minmod(thet*(Ui.velocity - Uimo.velocity), 0.5*(Uipo.velocity - Uimo.velocity),thet*(Uipo.velocity-Ui.velocity));
  UL.press = Ui.press + 0.5*dx*minmod(thet*(Ui.press - Uimo.press), 0.5*(Uipo.press - Uimo.press),thet*(Uipo.press-Ui.press));
  //  printf("%f\n", (UL.mass - Ui.mass)/Ui.mass);
  UL.energy = Et(UL.mass, UL.velocity, UL.press);
  return(UL);
}

Vars U_R(Vars Ui, Vars Uipo, Vars Uipt, double dx){
  Vars UR;
  UR.mass = Uipo.mass - 0.5*dx*minmod(thet*(Uipo.mass-Ui.mass),0.5*(Uipt.mass-Ui.mass),thet*(Uipt.mass-Uipo.mass));
  UR.velocity = Uipo.velocity - 0.5*dx*minmod(thet*(Uipo.velocity-Ui.velocity),.5*(Uipt.velocity-Ui.velocity),thet*(Uipt.velocity-Uipo.velocity));
  UR.press = Uipo.press - 0.5*dx*minmod(thet*(Uipo.press-Ui.press), 0.5*(Uipt.press-Ui.press), thet*(Uipt.press-Uipo.press));
  UR.energy = Et(UR.mass, UR.velocity, UR.press);
  return(UR);
}

void init_sys(int N, Vars * U, double dx, int type){
  if(type == 1){
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
  else if(type == 2){
    int i;
    for(i=0;i<N+4;++i){
      double x = (double)i*dx;
      if(x < 1./4.){
	U[i].mass = 3.856;
	U[i].velocity = 5.629;
	U[i].press = 10.33;
	//printf("yes %d\n",i);
      }
      else{
	U[i].mass = 1. - 0.2*sin(32*x);
	U[i].velocity = 0.0;
	U[i].press = 1.;
	U[i].energy = Et(U[i].mass, U[i].velocity, U[i].press);
      }
    }
    for(i=0;i<N;++i) U[i].energy = Et(U[i].mass, U[i].velocity, U[i].press);
  }
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

double HLL(int N, Vars * U, Flux * F_HLL, double dx){
  int i; double maxalph = 0.;
  for(i=2;i<N+3;++i){
    Vars UL = U_L(U[i-1], U[i-2], U[i], 1.); Vars UR = U_R(U[i-1], U[i], U[i+1], 1.); //Calculates the conserved variables at the interfaces
    //Vars UL, UR; UL = U[i-1]; UR = U[i];
    double alphap = alpha(UL, UR, 1.);
    double alpham = alpha(UL, UR, -1.);
    Flux FL = get_F(UL); Flux FR = get_F(UR); //Calculates the Fluxes from the left and right at the interface
    //Flux FL = get_F(U[i-1]); Flux FR = get_F(U[i]);
    F_HLL[i-2].rhov = (alphap*FL.rhov + alpham*FR.rhov - alphap*alpham*(UR.mass - UL.mass))/(alphap + alpham);
    F_HLL[i-2].mom = (alphap*FL.mom + alpham*FR.mom - alphap*alpham*(UR.velocity - UL.velocity))/(alphap + alpham);
    F_HLL[i-2].energy = (alphap*FL.energy + alpham*FR.energy - alphap*alpham*(UR.energy - UL.energy))/(alphap + alpham);
    maxalph = MAX(maxalph, alphap, alpham);
    //printf("alpha %d %f  %f %f\n", i, alphap, alpham, maxalph);
  }
  return(maxalph);
}

Flux hll( Vars * U, int i, double dx){
  Flux F_HLL;
  Vars UL = U_L(U[i-1], U[i-2], U[i], dx); Vars UR = U_R(U[i-1], U[i], U[i+1], dx); //Calculates the conserved variables at the interfaces
  //Vars UL, UR; UL = U[i-1]; UR=U[i];
  double alphap = alpha(UL, UR, 1.);
  double alpham = alpha(UL, UR, -1.);
  Flux FL = get_F(UL); Flux FR = get_F(UR); //Calculates the Fluxes from the left and right at the interface
  F_HLL.rhov = (alphap*FL.rhov + alpham*FR.rhov - alphap*alpham*(UR.mass - UL.mass))/(alphap + alpham);
  F_HLL.mom = (alphap*FL.mom + alpham*FR.mom - alphap*alpham*(UR.velocity - UL.velocity))/(alphap + alpham);
  F_HLL.energy = (alphap*FL.energy + alpham*FR.energy - alphap*alpham*(UR.energy - UL.energy))/(alphap + alpham);
  return(F_HLL);
}

Flux L(Flux FL, Flux FR, double dx){
  Flux LU;
  LU.rhov = -1.*(FR.rhov - FL.rhov)/dx;
  LU.mom = -1.*(FR.mom - FL.mom)/dx;
  LU.energy = -1.*(FR.energy - FL.energy)/dx;
  return(LU);
}

Vars make_U1(double dt, Vars U, Flux LU){
  Vars U1;
  U1.mass = U.mass + dt*LU.rhov;
  U1.velocity = U.velocity + dt*LU.mom;
  U1.energy = U.energy + dt*LU.energy;
  U1.press = (gam-1.)*(U1.energy - .5*pow(U1.velocity,2)/U1.mass);
  return(U1);
}

Vars make_U2(double dt, Vars U, Vars U1, Flux LU1){
  Vars U2;
  U2.mass = 0.75*U.mass + 0.25*U1.mass + 0.25*dt*LU1.rhov;
  U2.velocity = 0.75*U.velocity + 0.25*U1.velocity + 0.25*dt*LU1.mom;
  U2.energy = 0.75*U.energy + 0.25*U1.energy + 0.25*dt*LU1.energy;
  U2.press = (gam-1.)*(U2.energy - .5*pow(U2.velocity,2)/U2.mass);
  return(U2);
}

Vars make_UN(double dt, Vars U, Vars U2, Flux LU2){
  Vars UN;
  UN.mass = U.mass/3. + 2.*U2.mass/3. + 2.*dt*LU2.rhov/3.;
  UN.velocity = U.velocity/3. + 2.*U2.velocity/3. + 2.*dt*LU2.mom/3.;
  UN.energy = U.energy/3. + 2.*U2.energy/3. + 2.*dt*LU2.energy/3.;
  UN.press = (gam-1.)*(UN.energy - .5*pow(UN.velocity,2)/UN.mass);
  return(UN);
}

double Pressure(Vars U){
  double press;
  press = (gam-1.)*(U.energy - .5*pow(U.velocity,2)/U.mass);
  return(press);
}



double advance_system(int N, Vars * U, Flux * F_HLL, double dx, int type){
  double maxalpha = HLL(N, U, F_HLL, dx);//calculate the flux at the interfaces
  double dt = 0.5*dx/maxalpha;//calculate the time step
  Vars U_n[N+4], Un[N+4], U1, U2; int i;
  Flux FL, FR, LU;
  for(i=0;i<N+4;++i){
    U_n[i] = U[i];
  }
  for(i=2;i<N+2;++i){//make U1
    //FL = hll(U, i, dx); FR = hll(U, i+1, dx);
    FL = F_HLL[i-2]; FR = F_HLL[i-1];
    LU = L(FL, FR, dx);
    U_n[i] = make_U1(dt, U[i], LU);
  }
  for(i=0;i<N+4;++i){
    Un[i] = U_n[i];
  }
  for(i=2;i<N+2;++i){//make U2 U_n = U1
    FL = hll(U_n, i, 1.); FR = hll(U_n, i+1, 1.);
    //LU = L(FL, FR, dx);
    Un[i] = make_U2(dt, U[i], U_n[i], L(FL, FR, dx));
  }
  for(i=2;i<N+2;++i){//make U_n // Un = U2
    FL = hll(Un, i, 1.); FR = hll(Un, i+1, 1.);
    LU = L(FL, FR, dx);
    U_n[i] = make_UN(dt, U[i], Un[i], LU);
    }
  for(i=2;i<N+2;++i){
    U[i] = U_n[i];
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
  int i; for(i=0;i<N+4;++i){
    fprintf(fid, "%e %e  %e  %e  %e\n", -0.5 + i*dx, U[i].mass, U[i].velocity/U[i].mass, U[i].energy, U[i].press);
  }
}
int main(void){
  FILE * fid, * finit;
  fid = fopen("data_2nd.dat", "w");
  finit = fopen("init_data.dat", "w");
  int N = 300; Vars U[N+4]; Flux F_HLL[N+1];
  double T, t, dt, dx; T = .178; t = 0; dx = 1./(double)N;
  int type = 2;
  init_sys(N, U, dx, type); Write_Cons(N, U, dx, finit);
  //int j;for(j=0;j<2;++j){
  while(t<T){
    dt = advance_system(N, U, F_HLL, dx, type);
    t += dt; //printf("****t %f*****\n", t);
    //break;
  }
  //Print_Flux(N+1, F_HLL);
  Write_Cons(N, U, dx, fid);
}
