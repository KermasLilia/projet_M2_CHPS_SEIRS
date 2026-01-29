#include <iostream>   // entrées / sorties (cout)
#include <vector>     // vecteurs dynamiques
#include <cmath>      // fonctions math (exp)
#include <random>     // générateurs aléatoires
#include <fstream>    // fichiers
#include <iomanip>    // formatage sortie

#define N 20000        // nombre total d'agents
#define GRID 300       // taille de la grille (300x300)

struct Ag {
    int st; // état de l’agent  état (0=S, 1=E, 2=I, 3=R)
    int t;  // temps passé dans l’état courant
    int x, y;  // position sur la grilleg
    int dE, dI, dR;  // durées des états E, I et R
};

int grid[GRID][GRID]; // grille qui compte les infectés par case

// Fonction qui lance une simulation (une répétition)
void run(int id) {

    // Générateur aléatoire initialisé avec l'id
    std::mt19937 g(id*999+1);

    // Lois de probabilité
    std::uniform_real_distribution<> u(0,1);      // réel entre 0 et 1
    std::uniform_int_distribution<> d(0,GRID-1);  // entier position grille
    std::exponential_distribution<> eE(1.0/3.0);   // durée état E
    std::exponential_distribution<> eI(1.0/7.0);   // durée état I
    std::exponential_distribution<> eR(1.0/365.0); // durée état R

    // Création des agents
    std::vector<Ag> ags(N);

    for(int i=0;i<N;i++) {
        ags[i].x = d(g);     // position x aléatoire
        ags[i].y = d(g);     // position y aléatoire
        ags[i].t = 0;        // temps initial à 0
        ags[i].dE = (int)eE(g);
        ags[i].dI = (int)eI(g);
        ags[i].dR = (int)eR(g);

        // Les 20 premiers agents sont infectés
        ags[i].st = (i < 20) ? 2 : 0;
    }

    // Nom du fichier de sortie
    std::string fn = "res_cpp_" + std::to_string(id) + ".csv";
    if(id < 10) fn = "res_cpp_0" + std::to_string(id) + ".csv";

    std::ofstream f(fn);
    f << "jour,S,E,I,R\n"; // en-tête du fichier CSV

    // Boucle temporelle (730 jours)
    for(int t=0; t<=730; t++) {

        // Réinitialisation de la grille
        for(int x=0;x<GRID;x++)
            for(int y=0;y<GRID;y++)
                grid[x][y] = 0;

        int S=0, E=0, I=0, R=0;

        // Comptage des états et placement des infectés sur la grille
        for(auto &a : ags) {
            if(a.st == 0) S++;
            else if(a.st == 1) E++;
            else if(a.st == 2) {
                I++;
                grid[a.x][a.y]++; // on compte les infectés
            }
            else R++;
        }

        // Sauvegarde des résultats du jour
        f << t << "," << S << "," << E << "," << I << "," << R << "\n";

        // Déplacement aléatoire des agents
        for(auto &a : ags) {
            a.x = d(g);
            a.y = d(g);
        }

        // Mise à jour des états
        for(auto &a : ags) {
            a.t++; // incrément du temps dans l’état

            // Agent susceptible
            if(a.st == 0) {
                int ni = 0;

                // Comptage des infectés dans le voisinage
                for(int dx=-1; dx<=1; dx++)
                    for(int dy=-1; dy<=1; dy++)
                        ni += grid[(a.x+dx+GRID)%GRID][(a.y+dy+GRID)%GRID];

                // Probabilité d’infection
                if(ni > 0 && u(g) < (1.0 - exp(-0.5 * ni))) {
                    a.st = 1; // devient exposé
                    a.t = 0;
                }
            }

            // Transition E -> I
            else if(a.st == 1 && a.t > a.dE) {
                a.st = 2;
                a.t = 0;
            }

            // Transition I -> R
            else if(a.st == 2 && a.t > a.dI) {
                a.st = 3;
                a.t = 0;
            }

            // Transition R -> S (perte d’immunité)
            else if(a.st == 3 && a.t > a.dR) {
                a.st = 0;
                a.t = 0;
            }
        }
    }

    std::cout << "Cpp rep " << id << " done." << std::endl;
}

// Programme principal : lance 30 simulations
int main() {
    for(int i=1; i<=30; i++)
        run(i);
    return 0;
}
