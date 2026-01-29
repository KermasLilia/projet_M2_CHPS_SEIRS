#!/bin/bash

# ==============================================================================
# SCRIPT D'INSTALLATION COMPLET - PROJET SEIRS M2 CHPS
# ==============================================================================

echo ">>> Création de l'arborescence..."
mkdir -p Partie1
mkdir -p Partie2

# ==============================================================================
# PARTIE 1 : ÉQUATIONS DIFFÉRENTIELLES (ODE)
# ==============================================================================
echo ">>> Génération des fichiers Partie 1..."

# --- Python Euler ---
cat << 'EOF' > Partie1/seirs_python_euler.py
import csv
RHO, BETA, SIGMA, GAMMA = 1.0/365.0, 0.5, 1.0/3.0, 1.0/7.0
T_MAX, DT = 730.0, 0.01
S, E, I, R = 0.999, 0.0, 0.001, 0.0
t = 0.0

with open("results_py_euler.csv", "w", newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["jour","S","E","I","R"])
    steps = int(T_MAX/DT)
    for _ in range(steps+1):
        writer.writerow([t,S,E,I,R])
        dS = RHO*R - BETA*I*S
        dE = BETA*I*S - SIGMA*E
        dI = SIGMA*E - GAMMA*I
        dR = GAMMA*I - RHO*R
        S += dS*DT; E += dE*DT; I += dI*DT; R += dR*DT
        t += DT
print("Python Euler terminé.")
EOF

# --- Python RK4 ---
cat << 'EOF' > Partie1/seirs_python_rk4.py
import csv
RHO, BETA, SIGMA, GAMMA = 1.0/365.0, 0.5, 1.0/3.0, 1.0/7.0
T_MAX, DT = 730.0, 0.01
state = [0.999, 0.0, 0.001, 0.0] # S, E, I, R

def derivs(y):
    s, e, i, r = y
    return [RHO*r - BETA*i*s, BETA*i*s - SIGMA*e, SIGMA*e - GAMMA*i, GAMMA*i - RHO*r]

with open("results_py_rk4.csv", "w", newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["jour","S","E","I","R"])
    t = 0.0
    steps = int(T_MAX/DT)
    for _ in range(steps+1):
        writer.writerow([t]+state)
        k1 = derivs(state)
        y2 = [s + k*DT/2 for s,k in zip(state, k1)]
        k2 = derivs(y2)
        y3 = [s + k*DT/2 for s,k in zip(state, k2)]
        k3 = derivs(y3)
        y4 = [s + k*DT for s,k in zip(state, k3)]
        k4 = derivs(y4)
        state = [s + (DT/6)*(k1[i]+2*k2[i]+2*k3[i]+k4[i]) for i,s in enumerate(state)]
        t += DT
print("Python RK4 terminé.")
EOF

# --- C Euler ---
cat << 'EOF' > Partie1/seirs_c_euler.c
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
EOF

# --- C RK4 ---
cat << 'EOF' > Partie1/seirs_c_rk4.c
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
EOF

# --- Notebook Partie 1 ---
cat << 'EOF' > Partie1/gen_nb1.py
import json
nb = {
 "cells": [
  {"cell_type": "markdown", "source": ["# Analyse Partie 1 : Comparaison ODE"], "metadata": {}},
  {"cell_type": "code", "source": [
    "import pandas as pd\nimport matplotlib.pyplot as plt\nimport seaborn as sns\nplt.style.use('seaborn-v0_8-whitegrid')\nplt.rcParams['figure.figsize']=(10,6)\n\n# Chargement\ndf_py_rk4 = pd.read_csv('results_py_rk4.csv')\ndf_py_euler = pd.read_csv('results_py_euler.csv')\ndf_c_rk4 = pd.read_csv('results_c_rk4.csv')\n\n# Graphique 1 : Dynamique Globale (Couleurs significatives)\nplt.figure()\nplt.plot(df_py_rk4['jour'], df_py_rk4['S'], label='S (Susceptible)', color='#1f77b4', linewidth=2)\nplt.plot(df_py_rk4['jour'], df_py_rk4['E'], label='E (Exposed)', color='#ff7f0e', linewidth=2)\nplt.plot(df_py_rk4['jour'], df_py_rk4['I'], label='I (Infected)', color='#d62728', linewidth=2)\nplt.plot(df_py_rk4['jour'], df_py_rk4['R'], label='R (Recovered)', color='#2ca02c', linewidth=2)\nplt.title('Dynamique SEIRS (Méthode RK4)', fontsize=14)\nplt.xlabel('Jours')\nplt.ylabel('Proportion')\nplt.legend(frameon=True)\nplt.show()\n\n# Graphique 2 : Comparaison Euler vs RK4 (Zoom pic)\nplt.figure()\nplt.plot(df_py_euler['jour'], df_py_euler['I'], label='Euler (I)', color='gray', linestyle='--')\nplt.plot(df_py_rk4['jour'], df_py_rk4['I'], label='RK4 (I)', color='#d62728', linestyle='-')\nplt.xlim(0, 150)\nplt.title('Comparaison Précision : Euler vs RK4', fontsize=14)\nplt.xlabel('Jours')\nplt.ylabel('Infectés')\nplt.legend()\nplt.show()"
  ], "metadata": {}, "outputs": [], "execution_count": None}
 ],
 "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}}, "nbformat": 4, "nbformat_minor": 4
}
with open('Partie1/partie1_analyse.ipynb', 'w') as f: json.dump(nb, f, indent=1)
EOF
python3 Partie1/gen_nb1.py
rm Partie1/gen_nb1.py

# ==============================================================================
# PARTIE 2 : SYSTÈMES MULTI-AGENTS (SMA)
# ==============================================================================
echo ">>> Génération des fichiers Partie 2..."

# --- C Code ---
cat << 'EOF' > Partie2/seirs_abm.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define N 20000
#define GRID 300
#define TMAX 730
#define REPS 30
#define BETA 0.5
typedef struct { int state; int t; int x, y; int dE, dI, dR; } Ag;
int grid[GRID][GRID];
unsigned long mt[624]; int mti=625;
void init_mt(unsigned long s) { mt[0]=s; for(mti=1;mti<624;mti++) mt[mti]=(1812433253UL*(mt[mti-1]^(mt[mti-1]>>30))+mti); mti=0; }
unsigned long rand_mt() { 
    unsigned long y; if(mti>=624) init_mt(5489); 
    y=mt[mti++]; y^=(y>>11); y^=(y<<7)&0x9d2c5680UL; y^=(y<<15)&0xefc60000UL; y^=(y>>18); return y; 
}
double rand_real() { return rand_mt()*(1.0/4294967296.0); }
int nExp(double m) { return (int)(-m*log(1.0-rand_real())); }

void run(int id) {
    Ag *ags = malloc(N*sizeof(Ag));
    char fn[64]; sprintf(fn, "res_c_%02d.csv", id); FILE *f=fopen(fn, "w");
    fprintf(f, "jour,S,E,I,R\n");
    init_mt(id*123+1);
    for(int i=0;i<N;i++) {
        ags[i].x=rand_mt()%GRID; ags[i].y=rand_mt()%GRID; ags[i].t=0;
        ags[i].dE=nExp(3.0); ags[i].dI=nExp(7.0); ags[i].dR=nExp(365.0);
        ags[i].state = (i<20) ? 2 : 0;
    }
    for(int t=0; t<=TMAX; t++) {
        int S=0,E=0,I=0,R=0;
        for(int x=0;x<GRID;x++) for(int y=0;y<GRID;y++) grid[x][y]=0;
        for(int i=0;i<N;i++) {
            if(ags[i].state==0) S++; else if(ags[i].state==1) E++;
            else if(ags[i].state==2) { I++; grid[ags[i].x][ags[i].y]++; }
            else R++;
        }
        fprintf(f, "%d,%d,%d,%d,%d\n", t, S, E, I, R);
        for(int i=0;i<N;i++) { ags[i].x=rand_mt()%GRID; ags[i].y=rand_mt()%GRID; }
        for(int i=0;i<N;i++) {
            ags[i].t++;
            if(ags[i].state==0) {
                int ni=0;
                for(int dx=-1;dx<=1;dx++) for(int dy=-1;dy<=1;dy++) 
                    ni += grid[(ags[i].x+dx+GRID)%GRID][(ags[i].y+dy+GRID)%GRID];
                if(ni>0 && rand_real() < (1.0-exp(-0.5*ni))) { ags[i].state=1; ags[i].t=0; }
            } else if(ags[i].state==1 && ags[i].t>ags[i].dE) { ags[i].state=2; ags[i].t=0; }
            else if(ags[i].state==2 && ags[i].t>ags[i].dI) { ags[i].state=3; ags[i].t=0; }
            else if(ags[i].state==3 && ags[i].t>ags[i].dR) { ags[i].state=0; ags[i].t=0; }
        }
    }
    fclose(f); free(ags);
    printf("C rep %d done.\n", id);
}
int main() { for(int i=1;i<=REPS;i++) run(i); return 0; }
EOF

# --- C++ Code ---
cat << 'EOF' > Partie2/seirs_abm.cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <iomanip>
#define N 20000
#define GRID 300
struct Ag { int st; int t; int x,y; int dE,dI,dR; };
int grid[GRID][GRID];
void run(int id) {
    std::mt19937 g(id*999+1);
    std::uniform_real_distribution<> u(0,1);
    std::uniform_int_distribution<> d(0,GRID-1);
    std::exponential_distribution<> eE(1.0/3.0), eI(1.0/7.0), eR(1.0/365.0);
    std::vector<Ag> ags(N);
    for(int i=0;i<N;i++) {
        ags[i].x=d(g); ags[i].y=d(g); ags[i].t=0;
        ags[i].dE=(int)eE(g); ags[i].dI=(int)eI(g); ags[i].dR=(int)eR(g);
        ags[i].st=(i<20)?2:0;
    }
    std::string fn="res_cpp_"+std::to_string(id)+".csv"; if(id<10) fn="res_cpp_0"+std::to_string(id)+".csv";
    std::ofstream f(fn); f<<"jour,S,E,I,R\n";
    for(int t=0; t<=730; t++) {
        for(int x=0;x<GRID;x++) for(int y=0;y<GRID;y++) grid[x][y]=0;
        int S=0,E=0,I=0,R=0;
        for(auto &a:ags) {
            if(a.st==0) S++; else if(a.st==1) E++;
            else if(a.st==2) { I++; grid[a.x][a.y]++; } else R++;
        }
        f<<t<<","<<S<<","<<E<<","<<I<<","<<R<<"\n";
        for(auto &a:ags) { a.x=d(g); a.y=d(g); }
        for(auto &a:ags) {
            a.t++;
            if(a.st==0) {
                int ni=0;
                for(int dx=-1;dx<=1;dx++) for(int dy=-1;dy<=1;dy++) 
                    ni+=grid[(a.x+dx+GRID)%GRID][(a.y+dy+GRID)%GRID];
                if(ni>0 && u(g)<(1.0-exp(-0.5*ni))) { a.st=1; a.t=0; }
            } else if(a.st==1 && a.t>a.dE) { a.st=2; a.t=0; }
            else if(a.st==2 && a.t>a.dI) { a.st=3; a.t=0; }
            else if(a.st==3 && a.t>a.dR) { a.st=0; a.t=0; }
        }
    }
    std::cout << "Cpp rep " << id << " done." << std::endl;
}
int main() { for(int i=1;i<=30;i++) run(i); return 0; }
EOF

# --- Python Code ---
cat << 'EOF' > Partie2/seirs_abm.py
import numpy as np
import pandas as pd
from scipy.signal import convolve2d
from tqdm import tqdm

def run(id):
    np.random.seed(id*555)
    st = np.zeros(20000, dtype=int)
    st[:20] = 2 # I
    tm = np.zeros(20000, dtype=int)
    dE = np.random.exponential(3.0, 20000).astype(int)
    dI = np.random.exponential(7.0, 20000).astype(int)
    dR = np.random.exponential(365.0, 20000).astype(int)
    
    res = []
    kern = np.ones((3,3),int)
    
    for t in range(731):
        # Stats
        u,c = np.unique(st, return_counts=True)
        d = {0:0, 1:0, 2:0, 3:0}
        for k,v in zip(u,c): d[k]=v
        res.append([t, d[0], d[1], d[2], d[3]])
        
        # Deplacement
        px = np.random.randint(0,300, 20000)
        py = np.random.randint(0,300, 20000)
        tm += 1
        
        # Infection
        gI = np.zeros((300,300),int)
        maskI = (st==2)
        np.add.at(gI, (px[maskI], py[maskI]), 1)
        ni_g = convolve2d(gI, kern, mode='same', boundary='wrap')
        
        maskS = (st==0)
        if np.any(maskS):
            ni = ni_g[px[maskS], py[maskS]]
            prob = 1.0 - np.exp(-0.5 * ni)
            inf = np.random.random(len(ni)) < prob
            idx = np.where(maskS)[0][inf]
            st[idx]=1; tm[idx]=0
            
        # Transitions
        mE = (st==1)&(tm>dE); st[mE]=2; tm[mE]=0
        mI = (st==2)&(tm>dI); st[mI]=3; tm[mI]=0
        mR = (st==3)&(tm>dR); st[mR]=0; tm[mR]=0
        
    pd.DataFrame(res, columns=["jour","S","E","I","R"]).to_csv(f"res_py_{id:02d}.csv", index=False)

if __name__=="__main__":
    print("Start Python...")
    for i in tqdm(range(1,31)): run(i)
EOF

# --- Notebook Partie 2 ---
cat << 'EOF' > Partie2/gen_nb2.py
import json
nb = {
 "cells": [
  {"cell_type": "markdown", "source": ["# Analyse Partie 2 : Multi-Agents (SMA)"], "metadata": {}},
  {"cell_type": "code", "source": [
    "import pandas as pd\nimport numpy as np\nimport matplotlib.pyplot as plt\nimport seaborn as sns\nimport glob\nplt.style.use('seaborn-v0_8-whitegrid')\nplt.rcParams['figure.figsize']=(12,8)\n\ndef load(p): \n    fs = sorted(glob.glob(p))\n    l = [pd.read_csv(f) for f in fs]\n    return l\n\ndc = load('res_c_*.csv')\ndcpp = load('res_cpp_*.csv')\ndpy = load('res_py_*.csv')\n\n# Fonction Convergence\ndef plot_conv(data, lbl, col):\n    df = pd.concat(data)\n    grp = df.groupby(df.index)\n    mu = grp['I'].mean()/20000\n    std = grp['I'].std()/20000\n    plt.plot(mu, label=lbl, color=col, linewidth=2)\n    plt.fill_between(mu.index, mu-std, mu+std, color=col, alpha=0.2)\n\nplt.figure()\nplot_conv(dc, 'C', '#1f77b4')\nplot_conv(dcpp, 'C++', '#2ca02c')\nplot_conv(dpy, 'Python', '#d62728')\nplt.title('Convergence : Moyenne des infectés ± Écart-type', fontsize=16)\nplt.xlabel('Jours')\nplt.ylabel('Ratio Infectés')\nplt.legend()\nplt.show()\n\n# Boxplot\np_c = [d['I'].max()/20000 for d in dc]\np_cpp = [d['I'].max()/20000 for d in dcpp]\np_py = [d['I'].max()/20000 for d in dpy]\n\nplt.figure()\ndf_box = pd.DataFrame({'C':p_c, 'C++':p_cpp, 'Python':p_py})\nsns.boxplot(data=df_box, palette=['#1f77b4', '#2ca02c', '#d62728'])\nplt.title('Distribution des Pics Infectieux (30 réplications)', fontsize=16)\nplt.ylabel('Pic (Ratio)')\nplt.show()"
  ], "metadata": {}, "outputs": [], "execution_count": None}
 ],
 "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}}, "nbformat": 4, "nbformat_minor": 4
}
with open('Partie2/partie2_analyse.ipynb', 'w') as f: json.dump(nb, f, indent=1)
EOF
python3 Partie2/gen_nb2.py
rm Partie2/gen_nb2.py

# --- Makefiles ---
cat << 'EOF' > Partie1/Makefile
all: seirs_c_euler seirs_c_rk4
seirs_c_euler: seirs_c_euler.c
	gcc seirs_c_euler.c -o seirs_c_euler -O3
seirs_c_rk4: seirs_c_rk4.c
	gcc seirs_c_rk4.c -o seirs_c_rk4 -O3
clean:
	rm -f seirs_c_euler seirs_c_rk4
EOF
sed -i 's/^    /\t/' Partie1/Makefile

cat << 'EOF' > Partie2/Makefile
all: seirs_abm_c seirs_abm_cpp
seirs_abm_c: seirs_abm.c
	gcc seirs_abm.c -o seirs_abm_c -O3 -lm
seirs_abm_cpp: seirs_abm.cpp
	g++ seirs_abm.cpp -o seirs_abm_cpp -O3
clean:
	rm -f seirs_abm_c seirs_abm_cpp
EOF
sed -i 's/^    /\t/' Partie2/Makefile

echo ">>> INSTALLATION TERMINÉE !"
echo ">>> Suis les instructions de l'étape 2 pour exécuter."
