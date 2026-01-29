#include <stdio.h>
#define RHO (1.0/365.0)
#define BETA 0.5
#define SIGMA (1.0/3.0)
#define GAMMA (1.0/7.0)
int main() {
    double S=0.999, E=0.0, I=0.001, R=0.0, t=0.0, dt=0.01;
    FILE *f = fopen("results_c_euler.csv", "w");
    fprintf(f, "jour,S,E,I,R\n");
    for(int i=0; i<=73000; i++) {
        fprintf(f, "%f,%f,%f,%f,%f\n", t, S, E, I, R);
        double dS = RHO*R - BETA*I*S;
        double dE = BETA*I*S - SIGMA*E;
        double dI = SIGMA*E - GAMMA*I;
        double dR = GAMMA*I - RHO*R;
        S+=dS*dt; E+=dE*dt; I+=dI*dt; R+=dR*dt; t+=dt;
    }
    fclose(f); return 0;
}
