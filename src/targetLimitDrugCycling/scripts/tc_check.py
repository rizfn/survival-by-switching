import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

plt.rcParams.update({"font.size":20,"axes.labelsize":26,"xtick.labelsize":20,
"ytick.labelsize":20,"legend.fontsize":18,"lines.linewidth":2.8,"axes.linewidth":1.4,
"xtick.major.size":7,"ytick.major.size":7,"xtick.major.width":1.4,"ytick.major.width":1.4,
"figure.dpi":120,"savefig.bbox":"tight"})

# reduced target-cell model:  dT = lam - d T - a T I ;  dI = (a T - delta) I
l1,d1,a1,de1 = 4.0,1.0,1.0,5.0
l2,d2,a2,de2 = 0.2,1.0,5.0,1.25

def R0(l,d,a,de): return a*l/(d*de)
lb,db,ab,deb = (l1+l2)/2,(d1+d2)/2,(a1+a2)/2,(de1+de2)/2
print("R0_1 =",round(R0(l1,d1,a1,de1),3),"  R0_2 =",round(R0(l2,d2,a2,de2),3))
print("R0bar =",round(R0(lb,db,ab,deb),3))
print("T0_1=",l1/d1,"T*_1=",de1/a1,"   T0_2=",l2/d2,"T*_2=",de2/a2)
print("endemic I* (avg) =",round(lb/deb - db/ab,4))

def field(t,s,l,d,a,de):
    T,I=s; return [l-d*T-a*T*I,(a*T-de)*I]
def static(l,d,a,de,tend=60,s0=(2.0,0.05)):
    s=solve_ivp(field,(0,tend),s0,args=(l,d,a,de),t_eval=np.linspace(0,tend,3000),rtol=1e-10,atol=1e-13)
    return s.t,s.y
def switch(T,tend=60,s0=(2.0,0.05)):
    s=list(s0);t=0;ta=[0];Ia=[s0[1]];Ta=[s0[0]]
    while t<tend:
        for env in ((l1,d1,a1,de1),(l2,d2,a2,de2)):
            dt=T/2; tt=np.linspace(0,dt,max(40,int(300*dt)))
            sol=solve_ivp(field,(0,dt),s,args=env,t_eval=tt,rtol=1e-10,atol=1e-13)
            ta.extend(t+sol.t[1:]);Ia.extend(sol.y[1][1:]);Ta.extend(sol.y[0][1:])
            s=[sol.y[0][-1],sol.y[1][-1]];t+=dt
            if t>=tend:break
    return np.array(ta),np.array(Ta),np.array(Ia)

t1,Y1=static(l1,d1,a1,de1); t2,Y2=static(l2,d2,a2,de2)
ts,Ts,Is=switch(0.2)
print("I_final: env1=%.2e  env2=%.2e  switch=%.4f"%(Y1[1][-1],Y2[1][-1],Is[-1]))

fig,ax=plt.subplots(figsize=(9,6))
ax.semilogy(t1,Y1[1],color="#1b6ca8",label="drug 1 only  ($R_0{=}0.8$)")
ax.semilogy(t2,Y2[1],color="#c0392b",label="drug 2 only  ($R_0{=}0.8$)")
ax.semilogy(ts,Is,color="#2e7d32",label="fast switch  ($\\bar R_0{\\approx}2$)")
ax.set_xlabel("time"); ax.set_ylabel("infected cells  $I$")
ax.set_xlim(0,60); ax.set_ylim(1e-6,1)
ax.legend(frameon=False,loc="lower left")
tag="l1_4_d1_1_a1_1_de1_5_l2_0.2_d2_1_a2_5_de2_1.25_T_0.2"
fig.savefig(f"src/targetLimitDrugCycling/plots/{tag}.svg"); fig.savefig(f"src/targetLimitDrugCycling/plots/targetcell_switching_{tag}.pdf")
print("figure saved")