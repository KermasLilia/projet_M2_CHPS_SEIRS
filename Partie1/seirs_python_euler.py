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
