#!/usr/bin/env python3
"""Plot bifurcation/phase-line/potential and check explicit toggle and Hopf calculations."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.integrate import cumulative_trapezoid
OUT=Path(__file__).resolve().parents[1]/'figures'
plt.rcParams.update({'font.size':13,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight'})
a=.3
H=lambda x:x**4/(1+x**4)
Hp=lambda x:4*x**3/(1+x**4)**2
branch=lambda x:(x-a)/H(x)
fold1=brentq(lambda x:3*x-4*a-x**5,.4,.6)
fold2=brentq(lambda x:3*x-4*a-x**5,1,1.4)
b_up,b_down=branch(fold1),branch(fold2)
assert 0<b_down<b_up
b_ref=3.
f=lambda x:a+b_ref*H(x)-x
edges=np.linspace(0,5,3001)
roots=[brentq(f,l,r) for l,r in zip(edges[:-1],edges[1:]) if f(l)*f(r)<0]
assert len(roots)==3 and [b_ref*Hp(x)-1<0 for x in roots]==[True,False,True]
fig,axes=plt.subplots(1,3,figsize=(11.3,3.6))
for l,r,stable in [(a+.0001,fold1,True),(fold1,fold2,False),(fold2,8.3,True)]:
 xx=np.linspace(l,r,800);axes[0].plot(branch(xx),xx,'-' if stable else '--',color='C0' if stable else '.45',lw=1.6)
low_jump=brentq(lambda x:a+b_down*H(x)-x,a,fold1)
high_jump=brentq(lambda x:a+b_up*H(x)-x,fold2,8.3)
for b,l,r in [(b_up,fold1,high_jump),(b_down,fold2,low_jump)]:
 axes[0].annotate('',xy=(b,r),xytext=(b,l),arrowprops={'arrowstyle':'->','color':'C3','lw':1.3})
 axes[0].plot(b,fold1 if b==b_up else fold2,'ko',ms=3)
axes[0].axvline(b_ref,color='.7',lw=.8)
axes[0].set(xlabel='feedback strength b',ylabel='steady concentration x',xlim=(.8,7),ylim=(0,7.5),title='Stable branches and hysteresis')
axes[1].axhline(0,color='.6',lw=.8)
for x,stable in zip(roots,[True,False,True]):axes[1].plot(x,0,'ko',mfc='black' if stable else 'white',ms=7)
for l,r in zip([0]+roots,roots+[4]):
 m=(l+r)/2; direction=np.sign(f(m));span=min(.2,(r-l)*.3)
 axes[1].annotate('',xy=(m+direction*span,0),xytext=(m-direction*span,0),arrowprops={'arrowstyle':'->','color':'C0'})
axes[1].annotate('basin boundary',xy=(roots[1],.015),xytext=(1.8,.15),ha='center',fontsize=12,arrowprops={'arrowstyle':'->','color':'.4','lw':.8})
axes[1].set(xlabel='concentration x',xlim=(0,4),ylim=(-.2,.3),yticks=[],title='Two basins at fixed b=3')
xx=np.linspace(0,4,1600);V=-cumulative_trapezoid(f(xx),xx,initial=0)
axes[2].plot(xx,V,color='C0')
for x,stable in zip(roots,[True,False,True]):axes[2].plot(x,np.interp(x,xx,V),'ko',mfc='black' if stable else 'white',ms=5)
axes[2].set(xlabel='concentration x',ylabel='V(x)',title='Two wells and a barrier')
fig.tight_layout();fig.savefig(OUT/'bistability_branches.png',dpi=180);plt.close(fig)
# All positive n=3 toggle branches: each parameter above threshold has two stable asymmetric equilibria.
ac=3/2**(4/3)
for aa in [ac*.5,ac*.99,ac*1.01,2,3,5]:
 sym=brentq(lambda z:z*(1+z**3)-aa,0,aa)
 h=3*sym**3/(1+sym**3)
 assert (h<1)==(aa<ac)
 S=brentq(lambda z:z-z**-2-aa,.01,aa+2)
 if aa>ac:
  D=np.sqrt(S*S-4/S);x,y=(S+D)/2,(S-D)/2
  assert np.allclose([x*(1+y**3),y*(1+x**3)],[aa,aa])
  slopes=[-3*aa*y*y/(1+y**3)**2,-3*aa*x*x/(1+x**3)**2]
  J=np.array([[-1,slopes[0]],[slopes[1],-1]])
  assert np.all(np.linalg.eigvals(J)<0)
  assert np.allclose(slopes,[-3/(S**3-1)]*2)
# Hopf coefficient from the actual multilinear matrices, independent of its simplified scalar expression.
gamma=1.;sc=2**(1/3);omega=np.sqrt(3)*gamma
beta2=2*gamma*sc**2;beta3=-2*gamma*sc
P=np.array([[0,0,1],[1,0,0],[0,1,0]],dtype=complex)
A=-gamma*np.eye(3)-2*gamma*P
zeta=np.exp(-2j*np.pi/3);q=np.array([1,zeta.conjugate(),zeta]);p=q/3
B=lambda v,w:beta2*(P@v)*(P@w)
C=lambda v,w,z:beta3*(P@v)*(P@w)*(P@z)
assert np.allclose(A@q,1j*omega*q) and np.allclose(A.T@p,-1j*omega*p)
terms=C(q,q,q.conjugate())-2*B(q,np.linalg.solve(A,B(q,q.conjugate())))+B(q.conjugate(),np.linalg.solve(2j*omega*np.eye(3)-A,B(q,q)))
l1=np.vdot(p,terms).real/(2*omega)
assert np.isclose(l1,-5*sc/(6*np.sqrt(3))) and l1<0
# Independent explicit derivative at the toggle critical point confirms the cubic saturation.
st=2**(-1/3);at=3/2**(4/3)
f3=at*(-6/(1+st**3)**2+108*st**3/(1+st**3)**3-162*st**6/(1+st**3)**4)
assert np.isclose(f3,8*st)
print('Bistability branches, toggle thresholds/stability and negative Hopf coefficient verified')
