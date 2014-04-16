#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double gam = 1.4;//adiabatic index
double thet = 1.5;//slope limiter

typedef struct Vars{//conserved variables
  double mass;
  double xvelocity;
  double yvelocity;
  double energy;
  double press;
} Vars;

typedef struct Flux{
  double rhov;
  double momx;
  double momy;
  double energy;
} Flux;


void Print(int Nx, int Ny, Vars * U, double dx, double dy);

double Et(double rho, double xvel, double yvel, double P){//calculates the energy from pressure and other conserved variables
  return(.5*rho*(pow(xvel/rho,2)+pow(yvel/rho,2)) + P/(gam-1.));
}

double sgn(double x){//returns sign of the number
  if(x < 0) return(-1.);
  else if(x >0) return(1.);
  else return(0.);
}

double min(double x, double y, double z){//returns minium of set
  double num = x;
  if(num > y) num = y;
  if(num > z) num = z;
  return(num);
}
double minmod(double x, double y, double z){//slope limiter
  double minm = .25*fabs(sgn(x) + sgn(y))*(sgn(x) + sgn(z))*min(fabs(x), fabs(y), fabs(z));
  return(minm);
}

Vars U_L(Vars Ui, Vars Uimo, Vars Uipo){//interpolates the conserved variables at interface for the left state
  Vars UL;
  UL.mass = Ui.mass + 0.5*minmod(thet*(Ui.mass - Uimo.mass), 0.5*(Uipo.mass - Uimo.mass),thet*(Uipo.mass-Ui.mass));
  UL.xvelocity = Ui.xvelocity + 0.5*minmod(thet*(Ui.xvelocity - Uimo.xvelocity), 0.5*(Uipo.xvelocity - Uimo.xvelocity),thet*(Uipo.xvelocity-Ui.xvelocity));
  UL.yvelocity = Ui.yvelocity + 0.5*minmod(thet*(Ui.yvelocity - Uimo.yvelocity), 0.5*(Uipo.yvelocity - Uimo.yvelocity),thet*(Uipo.yvelocity-Ui.yvelocity));
  UL.press = Ui.press + 0.5*minmod(thet*(Ui.press - Uimo.press), 0.5*(Uipo.press - Uimo.press),thet*(Uipo.press-Ui.press));
  UL.energy = Et(UL.mass, UL.xvelocity, UL.yvelocity, UL.press);
  return(UL);
}

Vars U_R(Vars Ui, Vars Uipo, Vars Uipt){//does U_L for the right state
  Vars UR;
  UR.mass = Uipo.mass - 0.5*minmod(thet*(Uipo.mass-Ui.mass),0.5*(Uipt.mass-Ui.mass),thet*(Uipt.mass-Uipo.mass));
  UR.xvelocity = Uipo.xvelocity - 0.5*minmod(thet*(Uipo.xvelocity-Ui.xvelocity),.5*(Uipt.xvelocity-Ui.xvelocity),thet*(Uipt.xvelocity-Uipo.xvelocity));
  UR.yvelocity = Uipo.yvelocity - 0.5*minmod(thet*(Uipo.yvelocity-Ui.yvelocity),.5*(Uipt.yvelocity-Ui.yvelocity),thet*(Uipt.yvelocity-Uipo.yvelocity));
  UR.press = Uipo.press - 0.5*minmod(thet*(Uipo.press-Ui.press), 0.5*(Uipt.press-Ui.press), thet*(Uipt.press-Uipo.press));
  UR.energy = Et(UR.mass, UR.xvelocity, UR.yvelocity, UR.press);
  return(UR);
}

void init_sys(int Nx, int Ny, Vars * U, double dx, double dy){//various initial conditions
  int i, j;
  for(i=0;i<Nx;++i){
    for(j=0;j<Ny;++j){
      if((double)i*dx <= .5){
	U[i*Ny+j].mass = 1.0;
	U[i*Ny+j].xvelocity = 0.;
	U[i*Ny+j].yvelocity = 0.;
	U[i*Ny+j].press = 2.5;
      }
      else{
	U[i*Ny+j].mass = .1;
	U[i*Ny+j].xvelocity = -0.;
	U[i*Ny+j].yvelocity = 0.;
	U[i*Ny+j].press = .125;
      }
      //if(i*dx == .5) printf("here %f\n", U[i*Ny+j].mass);
    }
    //if(i*dx == .5) Print(Nx, Ny, U, dx, dy);
  }
  for(i=0;i<Nx;++i){
    for(j=0;j<Ny;++j){
      U[i*Ny+j].energy = Et(U[i*Ny+j].mass, U[i*Ny+j].xvelocity, U[i*Ny+j].yvelocity, U[i*Ny+j].press);
      //if(i*dx == .5) printf("%f\n", U[i*Ny+j].mass);
    }
    }
  //Print(Nx, Ny, U, dx, dy);
}

double lambda(double vel, Vars U, double pm){
  double c = sqrt(gam*U.press/U.mass);
  return(vel + pm*c);
}
double MAX(double a, double b, double c){
  double max = a;
  if(max < b) max = b;
  if(max < c) max = c;
  return(max);
}

double alpha(Vars UL, Vars UR, double pm, int xy){
  double alph;
  if(xy == 1){// x integration
    alph = MAX(0., pm*lambda(UR.xvelocity, UR, pm), pm*lambda(UL.xvelocity, UL,pm));
  }
  else if (xy == -1){// y integration
    alph = MAX(0., pm*lambda(UR.yvelocity, UR, pm), pm*lambda(UL.yvelocity, UL,pm));
  }
  return(alph);
  }

Flux get_F(Vars U){//get the flux from the conserved variables
  Flux F;  
  F.rhov = U.xvelocity;
  F.momx = pow(U.xvelocity,2)/U.mass + U.press;
  F.momy = U.xvelocity*U.yvelocity/U.mass;
  F.energy = (U.energy + U.press)*U.xvelocity/U.mass;
  return(F);
}

Flux get_G(Vars U){
  Flux G;
  G.rhov = U.xvelocity;
  G.momx = U.yvelocity*U.xvelocity/U.mass;
  G.momy = pow(U.xvelocity,2)/U.mass + U.press;
  G.energy = (U.energy + U.press)*U.yvelocity/U.mass;

}
/*
double Maxalpha(int N, Vars * U){//calculates maxalpha for stability condition
  int i; double maxalph = 0.;
  for(i=2;i<N+3;++i){
    Vars UL = U_L(U[i-1], U[i-2], U[i]); Vars UR = U_R(U[i-1], U[i], U[i+1]); //Calculates the conserved variables at the interfaces
    double alphap = alpha(UL, UR, 1.);
    double alpham = alpha(UL, UR, -1.);
    maxalph = MAX(maxalph, alphap, alpham);
  }
  return(maxalph);
}
*/
Flux Fhll( Vars Uimo, Vars Ui, Vars Uipo, Vars Uipt){//calculates the HLL flux for a given interface in x
  Flux F_HLL;
  Vars UL = U_L(Ui, Uimo, Uipo); Vars UR = U_R(Ui, Uipo, Uipt); //Calculates the conserved variables at the interfaces
  double alphap = alpha(UL, UR, 1., 1);
  double alpham = alpha(UL, UR, -1., 1);
  Flux FL = get_F(UL); Flux FR = get_F(UR); //Calculates the Fluxes from the left and right at the interface
  F_HLL.rhov = (alphap*FL.rhov + alpham*FR.rhov - alphap*alpham*(UR.mass - UL.mass))/(alphap + alpham);
  F_HLL.momx = (alphap*FL.momx + alpham*FR.momx - alphap*alpham*(UR.xvelocity - UL.xvelocity))/(alphap + alpham);
  F_HLL.momy = (alphap*FL.momy + alpham*FR.momy - alphap*alpham*(UR.yvelocity - UL.yvelocity))/(alphap + alpham);
  F_HLL.energy = (alphap*FL.energy + alpham*FR.energy - alphap*alpham*(UR.energy - UL.energy))/(alphap + alpham);
  return(F_HLL);
}

Flux Ghll( Vars Uimo, Vars Ui, Vars Uipo, Vars Uipt){//calculates the HLL flux for a given interface in y
  Flux G_HLL;
  Vars UL = U_L(Ui, Uimo, Uipo); Vars UR = U_R(Ui, Uipo, Uipt); //Calculates the conserved variables at the interfaces
  double alphap = alpha(UL, UR, 1., -1);
  double alpham = alpha(UL, UR, -1., -1);
  Flux FL = get_G(UL); Flux FR = get_G(UR); //Calculates the Fluxes from the left and right at the interface
  G_HLL.rhov = (alphap*FL.rhov + alpham*FR.rhov - alphap*alpham*(UR.mass - UL.mass))/(alphap + alpham);
  G_HLL.momx = (alphap*FL.momx + alpham*FR.momx - alphap*alpham*(UR.xvelocity - UL.xvelocity))/(alphap + alpham);
  G_HLL.momy = (alphap*FL.momy + alpham*FR.momy - alphap*alpham*(UR.yvelocity - UL.yvelocity))/(alphap + alpham);
  G_HLL.energy = (alphap*FL.energy + alpham*FR.energy - alphap*alpham*(UR.energy - UL.energy))/(alphap + alpham);
  return(G_HLL);
}

Flux L(Flux FL, Flux FR,  Flux GL, Flux GR, double dy, double dx){
  Flux LU;
  LU.rhov = -1.*(FR.rhov - FL.rhov)/dx - 1.*(GR.rhov - GL.rhov)/dy;
  LU.momx = -1.*(FR.momx - FL.momx)/dx - 1.*(GR.momx - GL.momx)/dy;
  LU.momy = -1.*(FR.momy - FL.momy)/dx - 1.*(GR.momy - GL.momy)/dy;
  LU.energy = -1.*(FR.energy - FL.energy)/dx - 1.*(GR.energy - GL.energy)/dy;
  return(LU);
}

Vars make_U1(double dt, Vars U, Flux LU){
  Vars U1;
  U1.mass = U.mass + dt*LU.rhov;
  U1.xvelocity = U.xvelocity + dt*LU.momx;
  U1.yvelocity = U.yvelocity + dt*LU.momy;
  U1.energy = U.energy + dt*LU.energy;
  U1.press = (gam-1.)*(U1.energy - .5*(pow(U1.yvelocity,2) + pow(U1.xvelocity,2))/U1.mass);
  return(U1);
}
/*
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
*/

void Set_B2A(int Nx, int Ny, Vars * A, Vars * B){
  int i, j;
  for(i=0;i<Nx;++i){
    for(j=0;j<Ny;++j){
      B[i*Ny+j] = A[i*Ny+j];
    }
  }
}

double advance_system(int Nx, int Ny, Vars * U, double dx, double dy){
  //double maxalpha = Maxalpha(Nx, Ny, U);//calculate the maxalpha
  double dt = 0.001;//0.5*dx/maxalpha;//calculate the time step
  Vars U_n[Nx*Ny], Un[Nx*Ny]; int i, j;
  Flux FL, FR, GL, GR;
  Set_B2A(Nx, Ny, U, U_n);
  for(i=2;i<Nx-2;++i){
    for(j=2;j<Ny-2;++j){
      FL = Fhll(U[(i-2)*Ny+j],U[(i-1)*Ny+j],U[i*Ny+j],U[(i+1)*Ny+j]); 
      FR = Fhll(U[(i-1)*Ny+j],U[i*Ny+j],U[(i+1)*Ny+j],U[(i+2)*Ny+j]);
      GL = Ghll(U[i*Ny+j-2],U[i*Ny+j-1],U[i*Ny+j],U[i*Ny+j+1]);
      GR = Ghll(U[i*Ny+j-1],U[i*Ny+j],U[i*Ny+j+1],U[i*Ny+j+2]);
      U_n[i*Ny+j] = make_U1(dt, U[i*Ny+j], L(FL, FR, GL, GR, dy, dx));
      
    }
  }
  for(i=2;i<Nx;++i){
    for(j=2;j<Ny;++j){
      U[i*Ny+j] = U_n[i*Ny+j];
    }
  }
  for(i=0;i<Nx;++i){
      U[i*Ny+1] = U[i*Ny+0];
      U[i*Ny+Nx-1] = U[i*Ny+Nx-2];
      U[i*Ny+0] = U[i*Ny+Nx-1];
      U[i*Ny+Nx-2] = U[i*Ny+Nx-3];
    }
  /*
  for(i=0;i<N+4;++i){
    U_n[i] = U[i];
  }
  for(i=2;i<N+2;++i){//make U1
    FL = hll(U, i); FR = hll(U, i+1);
    U_n[i] = make_U1(dt, U[i], L(FL, FR, dx));
  }
  for(i=0;i<N+4;++i){
    Un[i] = U_n[i];
  }
  for(i=2;i<N+2;++i){//make U2 U_n = U1
    FL = hll(U_n, i); FR = hll(U_n, i+1);
    Un[i] = make_U2(dt, U[i], U_n[i], L(FL, FR, dx));
  }
  for(i=2;i<N+2;++i){//make U_n // Un = U2
    FL = hll(Un, i); FR = hll(Un, i+1);
    U_n[i] = make_UN(dt, U[i], Un[i], L(FL, FR, dx));
    }
  for(i=2;i<N+2;++i){
    U[i] = U_n[i];
  }
  */
  return(dt);
}

void Write_Cons(int Nx, int Ny, Vars * U, double dx, double dy, FILE * fid){
  int i, j; for(i=0;i<Nx;++i){
    for(j=2;j<Ny-2;++j){
      // printf("%d %d\n", i, j);
      fprintf(fid, "%e %e  %e  %e  %e  %e  %e\n", i*dx, j*dy, U[i*Ny+j].mass, U[i*Ny+j].xvelocity/U[i*Ny+j].mass, U[i*Ny+j].yvelocity/U[i*Ny+j].mass, U[i*Ny+j].energy, U[i*Ny+j].press);
    }
    fprintf(fid,"\n");
  }
}
void Print(int Nx, int Ny, Vars * U, double dx, double dy){
  int i, j; for(i=0;i<Nx;++i){
    for(j=0;j<Ny;++j){
      //if(i*dx == .5) printf("%f\n", U[i*Ny+j].mass);
      printf( "%e %e  %e  %e  %e  %e  %e\n", i*dx, j*dy, U[i*Ny+j].mass, U[i*Ny+j].xvelocity/U[i*Ny+j].mass, U[i*Ny+j].yvelocity/U[i*Ny+j].mass, U[i*Ny+j].energy, U[i*Ny+j].press);
    }
    printf("\n");
  }
}

int main(void){
  FILE * fid, * finit;
  fid = fopen("data_2D.dat", "w");
  finit = fopen("init_data_2D.dat", "w");
  int Nx = 100; int Ny = 100; 
  Vars U[(Nx)*(Ny)]; double delt = 0.001;
  double T, t, dt, dx, dy; T = .15; t = 0; dx = 1./(double)(Nx); dy = 1./(double)(Ny);
  int type = 1;
  init_sys(Nx, Ny, U, dx, dy);  Write_Cons(Nx, Ny, U, dx, dy, finit); fclose(finit);
  int count = 0; char str[80];
  while(t<T){
    dt = advance_system(Nx, Ny, U, dx, dy);
    t += dt;
    /*if(count % 2 == 0){
      sprintf(str, "T_%d.dat", count);
      fid = fopen(str, "w");
      Write_Cons(N, U, dx, fid);
      fclose(fid);
    }
    count += 1;*/
    //break;
  }
  //printf("here\n");
  //Print(Nx, Ny, U, dx, dy);
  Write_Cons(Nx, Ny, U, dx, dy, fid);
  fclose(fid);
}
