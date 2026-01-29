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

// Fonction exponentielle pour les durées
int nExp(double m) { 
    return (int)(-m * log(1.0 - ((double)rand() / (double)RAND_MAX))); 
}

void run(int id) {
    // Initialisation du hasard avec une graine unique par réplication
    srand(id * 54321 + time(NULL));

    Ag *ags = malloc(N * sizeof(Ag));
    char fn[64]; sprintf(fn, "res_c_%02d.csv", id); 
    FILE *f = fopen(fn, "w");
    fprintf(f, "jour,S,E,I,R\n");

    // Initialisation Agents
    for(int i=0; i<N; i++) {
        ags[i].x = rand() % GRID; 
        ags[i].y = rand() % GRID; 
        ags[i].t = 0;
        ags[i].dE = nExp(3.0); 
        ags[i].dI = nExp(7.0); 
        ags[i].dR = nExp(365.0);
        // 20 premiers infectés
        ags[i].state = (i < 20) ? 2 : 0;
    }

    // Boucle temporelle
    for(int t=0; t<=TMAX; t++) {
        int S=0, E=0, I=0, R=0;
        
        // Reset grille et comptage
        for(int x=0; x<GRID; x++) for(int y=0; y<GRID; y++) grid[x][y] = 0;
        
        for(int i=0; i<N; i++) {
            if(ags[i].state == 0) S++;
            else if(ags[i].state == 1) E++;
            else if(ags[i].state == 2) { I++; grid[ags[i].x][ags[i].y]++; }
            else R++;
        }
        fprintf(f, "%d,%d,%d,%d,%d\n", t, S, E, I, R);

        // Déplacement
        for(int i=0; i<N; i++) { 
            ags[i].x = rand() % GRID; 
            ags[i].y = rand() % GRID; 
        }

        // Mises à jour états
        for(int i=0; i<N; i++) {
            ags[i].t++;
            if(ags[i].state == 0) { // S -> E
                int ni = 0;
                // Voisinage de Moore
                for(int dx=-1; dx<=1; dx++) {
                    for(int dy=-1; dy<=1; dy++) {
                         // Modulo pour tore
                         int nx = (ags[i].x + dx + GRID) % GRID;
                         int ny = (ags[i].y + dy + GRID) % GRID;
                         ni += grid[nx][ny];
                    }
                }
                if(ni > 0) {
                    double prob = 1.0 - exp(-0.5 * ni);
                    if( ((double)rand() / (double)RAND_MAX) < prob ) {
                        ags[i].state = 1; 
                        ags[i].t = 0;
                    }
                }
            } 
            else if(ags[i].state == 1 && ags[i].t > ags[i].dE) { ags[i].state = 2; ags[i].t = 0; } // E->I
            else if(ags[i].state == 2 && ags[i].t > ags[i].dI) { ags[i].state = 3; ags[i].t = 0; } // I->R
            else if(ags[i].state == 3 && ags[i].t > ags[i].dR) { ags[i].state = 0; ags[i].t = 0; } // R->S
        }
    }
    fclose(f); 
    free(ags);
    printf("Correction C rep %d done.\n", id);
}

int main() { 
    for(int i=1; i<=30; i++) run(i); 
    return 0; 
}
