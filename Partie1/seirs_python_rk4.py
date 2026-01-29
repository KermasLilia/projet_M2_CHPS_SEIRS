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
