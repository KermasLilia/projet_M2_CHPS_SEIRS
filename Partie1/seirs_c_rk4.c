#include <stdio.h>
#define RHO (1.0/365.0)
#define BETA 0.5
#define SIGMA (1.0/3.0)
#define GAMMA (1.0/7.0)
typedef struct { double S,E,I,R; } St;
St d(St y) { 
    St r; 
    r.S=RHO*y.R - BETA*y.I*y.S; 
    r.E=BETA*y.I*y.S - SIGMA*y.E; 
    r.I=SIGMA*y.E - GAMMA*y.I; 
    r.R=GAMMA*y.I - RHO*y.R; 
    return r; 
}
int main() {
    St y={0.999, 0.0, 0.001, 0.0}; double t=0.0, dt=0.01;
    FILE *f = fopen("results_c_rk4.csv", "w");
    fprintf(f, "jour,S,E,I,R\n");
    for(int i=0; i<=73000; i++) {
        fprintf(f, "%f,%f,%f,%f,%f\n", t, y.S, y.E, y.I, y.R);
        St k1=d(y);
        St y2={y.S+k1.S*dt/2, y.E+k1.E*dt/2, y.I+k1.I*dt/2, y.R+k1.R*dt/2}; St k2=d(y2);
        St y3={y.S+k2.S*dt/2, y.E+k2.E*dt/2, y.I+k2.I*dt/2, y.R+k2.R*dt/2}; St k3=d(y3);
        St y4={y.S+k3.S*dt, y.E+k3.E*dt, y.I+k3.I*dt, y.R+k3.R*dt}; St k4=d(y4);
        y.S+=(dt/6)*(k1.S+2*k2.S+2*k3.S+k4.S);
        y.E+=(dt/6)*(k1.E+2*k2.E+2*k3.E+k4.E);
        y.I+=(dt/6)*(k1.I+2*k2.I+2*k3.I+k4.I);
        y.R+=(dt/6)*(k1.R+2*k2.R+2*k3.R+k4.R);
        t+=dt;
    }
    fclose(f); return 0;
}
