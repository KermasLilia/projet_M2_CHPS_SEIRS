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
