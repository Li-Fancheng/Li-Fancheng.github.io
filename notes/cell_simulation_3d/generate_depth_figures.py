#!/usr/bin/env python3
"""Companion calculations for the expanded notes; no experimental data.

Run generate_figures.py first. This script adds nine figure policies to its manifest.
"""
from pathlib import Path
import sys
sys.dont_write_bytecode = True
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import generate_figures as base

ROOT=Path(__file__).resolve().parent
POLICY=ROOT.parents[1]/'data/cell_simulation_3d_figures.json'
MANIFEST=json.loads(POLICY.read_text())
COLORS=base.COLORS


def save(fig,slug,name):
    folder=ROOT/slug/'figures';folder.mkdir(parents=True,exist_ok=True)
    fig.savefig(folder/f'{name}.svg',metadata={'Date':None})
    fig.savefig(folder/f'{name}.png',dpi=300)
    w,h=fig.get_size_inches()
    MANIFEST[f'{slug}/figures/{name}']={'width':820,'kind':'calculation',
        'intrinsic_width':round(w*96),'intrinsic_height':round(h*96)}
    plt.close(fig)


def contact_balance():
    fig,axs=base.pair(['A  Local force balance at a contact line','B  Adhesion work and contact angle'])
    ax=axs[0];theta=np.deg2rad(55);gamma=1.;contact=2*np.cos(theta)
    vecs=[(np.cos(theta),np.sin(theta)),(np.cos(theta),-np.sin(theta)),(-contact,0)]
    for vec,col in zip(vecs,COLORS):
        ax.annotate('',xy=vec,xytext=(0,0),arrowprops={'arrowstyle':'->','lw':2.6,'color':col})
    ax.plot([0,1.2],[0,0],':',color='#64748b',lw=1)
    ax.add_patch(Arc((0,0),.65,.65,theta1=0,theta2=55,color='#64748b'))
    ax.text(.36,.16,r'$\theta$');ax.text(.59,.86,r'$\gamma\mathbf{t}_1$')
    ax.text(.59,-.93,r'$\gamma\mathbf{t}_2$');ax.text(-1.18,.16,r'$\Gamma\mathbf{e}$')
    ax.set(xlim=(-1.4,1.35),ylim=(-1.1,1.1),aspect='equal');ax.axis('off')
    angle=np.linspace(0,90,300);work=2*(1-np.cos(np.deg2rad(angle)))
    axs[1].plot(angle,work);axs[1].set(xlabel=r'Angle $\theta$ (degrees)',ylabel=r'Adhesion ratio $W/\gamma$',xticks=[0,30,60,90])
    assert np.linalg.norm(np.sum(vecs,axis=0))<1e-14
    save(fig,'foundations','contact_balance')


def interface_motion():
    fig,axs=base.pair(['A  Keep the bare field mobility fixed','B  Match the interface mobility'])
    t=np.linspace(0,.1,400);epsilons=[.01,.02,.04];cw=1/6;bare=cw/.02
    for eps in epsilons:
        mobility=bare*eps/cw
        axs[0].plot(t,np.sqrt(1-4*mobility*t),label=fr'$\epsilon={eps}$')
    axs[1].plot(t,np.sqrt(1-4*t),lw=2.4,label=r'$\mathcal{M}=M\epsilon/c_W=1$')
    axs[1].text(.009,.64,'All three interface thicknesses coincide',fontsize=10)
    for ax in axs:
        ax.set(xlabel='Time $t$',ylabel='Sphere radius $R$',ylim=(.4,1.03));ax.legend(frameon=False)
    save(fig,'phase_field','interface_motion')


def surface_shear():
    fig,axs=base.pair(['A  Equal-area shear of a triangular patch','B  Same geometry, different material'])
    pts=[];idx={}
    for j in range(3):
        for i in range(3):
            idx[(i,j)]=len(pts);pts.append([i+.5*j,np.sqrt(3)*j/2])
    pts=np.array(pts);shear=np.array([[1,.45],[0,1]])
    for transform,color,alpha in [(np.eye(2),'#94a3b8',.5),(shear,COLORS[0],1)]:
        q=pts@transform.T
        for (i,j),a in idx.items():
            for step in [(1,0),(0,1),(-1,1)]:
                target=(i+step[0],j+step[1])
                if target in idx:
                    b=idx[target];axs[0].plot(q[[a,b],0],q[[a,b],1],color=color,alpha=alpha,lw=1.2)
        axs[0].scatter(*q.T,color=color,s=13,alpha=alpha)
    axs[0].set(aspect='equal',xlabel='$x$',ylabel='$y$');axs[0].grid(False)
    s=np.linspace(0,.9,400);angle=np.array([0,np.pi/3,2*np.pi/3])
    stretch=np.sqrt((np.cos(angle)[None,:]+s[:,None]*np.sin(angle)[None,:])**2+np.sin(angle)[None,:]**2)
    energy=np.sum((stretch-1)**2,axis=1)/np.sqrt(3)
    axs[1].plot(s,energy,label='Fixed-reference spring network')
    axs[1].plot(s,np.sqrt(3)*s*s/8,'--',label='Small-strain limit')
    axs[1].plot(s,np.zeros_like(s),color=COLORS[2],label='Area-tension membrane')
    axs[1].set(xlabel='Shear $s$',ylabel='Energy increase per area / $k$');axs[1].legend(frameon=False,fontsize=9)
    assert abs(energy[1]/s[1]**2-np.sqrt(3)/8)<.002
    assert abs(np.linalg.det(shear)-1)<1e-14
    save(fig,'surface_mesh','surface_shear')


def shape_index():
    fig=plt.figure(figsize=(10,4.25),layout='constrained')
    ax=fig.add_subplot(1,2,1,projection='3d')
    v,f=base.sphere_mesh(2);_,vol=base.geometry(v,f);v=v/vol**(1/3)
    ax.add_collection3d(Poly3DCollection((v+[-2,0,0])[f],facecolor=COLORS[0],edgecolor='none',alpha=.7))
    cube=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]],float)-.5
    cf=[[0,3,2,1],[4,5,6,7],[0,1,5,4],[1,2,6,5],[2,3,7,6],[3,0,4,7]]
    ax.add_collection3d(Poly3DCollection(cube[cf],facecolor=COLORS[1],edgecolor='#334155',linewidth=.6,alpha=.6))
    v,f=base.sphere_mesh(0);_,vol=base.geometry(v,f);v=v/vol**(1/3)
    ax.add_collection3d(Poly3DCollection((v+[2,0,0])[f],facecolor=COLORS[2],edgecolor='#334155',linewidth=.6,alpha=.65))
    ax.set(xlim=(-3,3),ylim=(-1.3,1.3),zlim=(-1,1),title='A  Three shapes at equal volume')
    ax.set_box_aspect((6,2.6,2));ax.view_init(22,-65);ax.set_axis_off()
    ax=fig.add_subplot(1,2,2)
    values=[(36*np.pi)**(1/3),6,2*np.sqrt(3)/(np.sqrt(2)/3)**(2/3)]
    ax.bar(['Sphere','Cube','Octahedron'],values,color=COLORS,width=.6)
    for i,value in enumerate(values):ax.text(i,value+.09,f'{value:.3f}',ha='center')
    ax.set(ylim=(0,6.7),ylabel=r'Shape index $A/V^{2/3}$',title='B  Unchanged by uniform rescaling')
    ax.grid(axis='y',alpha=.16)
    save(fig,'vertex_3d','shape_index')


def potts_volume():
    fig,axs=base.pair(['A  Match a fixed physical volume energy','B  Exact cost of adding or removing one voxel'])
    bulk=5.;target=1.
    for h in [.08,.12,.18]:
        n=np.arange(np.ceil(.7/h**3),np.floor(1.3/h**3)+1, max(1,round(.008/h**3)))
        physical=n*h**3;lam=bulk*h**6/(2*target);v0=target/h**3
        discrete=lam*(n-v0)**2
        axs[0].plot(physical,discrete,label=fr'$h={h}$')
        assert np.max(abs(discrete-bulk*(physical-target)**2/(2*target)))<1e-12
    axs[0].set(xlabel='Physical volume $V$',ylabel='Volume energy $E_V$');axs[0].legend(frameon=False)
    h=.12;vol=np.linspace(.75,1.25,400);lam=bulk*h**6/2;n=vol/h**3;v0=1/h**3
    add=lam*(2*(n-v0)+1)/(bulk*h**3)
    remove=lam*(-2*(n-v0)+1)/(bulk*h**3)
    axs[1].plot(vol,add,label='Add one voxel');axs[1].plot(vol,remove,label='Remove one voxel')
    axs[1].axhline(0,color='#64748b',lw=.7)
    axs[1].set(xlabel='Physical volume $V$',ylabel=r'$\Delta E_V/(K_Vh^3)$');axs[1].legend(frameon=False)
    save(fig,'lattice_potts','potts_volume')


def transport_motility():
    fig,axs=base.pair(['A  Steady nutrient profiles in a sphere','B  Equal long-time active diffusivity'])
    r=np.linspace(0,1,400)
    for radius in [.4,.7,.95]:
        axs[0].plot(r,1-radius**2*(1-r*r),label=fr'$R={radius}$')
    axs[0].axhline(.3,color='#64748b',ls=':',label='Example threshold')
    axs[0].set(xlabel='Relative radial position $r/R$',ylabel='Concentration $c$',ylim=(0,1.04));axs[0].legend(frameon=False)
    t=np.geomspace(.002,50,500)
    for tau in [.1,1,5]:
        msd=6*(t+tau*np.expm1(-t/tau))
        axs[1].loglog(t,msd,label=fr'$\tau_p={tau}$')
        assert np.all(msd>0)
    axs[1].loglog(t,6*t,':',color='#64748b',label=r'$6D_{act}t$')
    axs[1].set(xlabel='Time $t$',ylabel='Mean-square displacement');axs[1].legend(frameon=False,fontsize=9)
    save(fig,'particles','transport_motility')


def rheology_scales():
    fig,axs=base.pair(['A  Maxwell storage and loss moduli','B  Length-dependent versus local relaxation'])
    freq=np.geomspace(.01,100,400)
    axs[0].loglog(freq,freq**2/(1+freq**2),label=r"$G'/G$")
    axs[0].loglog(freq,freq/(1+freq**2),label=r"$G''/G$")
    axs[0].set(xlabel=r'Frequency $\omega\tau$',ylabel='Normalized modulus');axs[0].legend(frameon=False)
    length=np.geomspace(.1,10,400)
    axs[1].loglog(length,length**2,label=r'Poroelastic: $\tau\propto\ell^2$')
    axs[1].loglog(length,np.ones_like(length),label='Fixed local material time')
    axs[1].set(xlabel=r'Perturbation length $\ell/\ell_0$',ylabel=r'Relaxation time $\tau/\tau_0$');axs[1].legend(frameon=False,fontsize=9)
    save(fig,'continuum_fluids','rheology_scales')


def buckling():
    fig,axs=base.pair(['A  Compression softens an out-of-plane mode','B  Foundation stiffness selects a wavelength'])
    q=np.linspace(0,2,400)
    for compression in [1.,2.,2.6]:
        axs[0].plot(q,q**4-compression*q*q+1,label=fr'$N={compression}$')
    axs[0].axhline(0,color='#64748b',lw=.8)
    axs[0].set(xlabel='Wavenumber $q$',ylabel=r'Mode stiffness $\mathcal{K}(q)$',ylim=(-1,9));axs[0].legend(frameon=False)
    k=np.geomspace(.01,100,400)
    axs[1].loglog(k,2*np.pi/k**.25)
    axs[1].set(xlabel='Foundation stiffness $K_f$',ylabel=r'Critical wavelength $2\pi/q_c$')
    assert abs(1-2+1)<1e-14
    save(fig,'morphogenesis','buckling')


def curvature_noise():
    fig,axs=base.pair(['A  Noise and smoothing set a useful scale','B  Bias and noise contributions'])
    h=np.geomspace(.015,.5,400)
    for sigma in [1e-4,5e-4,2e-3]:
        rms=np.sqrt(h**4+sigma**2/h**4)
        axs[0].loglog(h,rms,label=fr'$\sigma_x={sigma:g}$')
        star=sigma**.25
        axs[0].plot(star,np.sqrt(star**4+sigma**2/star**4),'o',color=axs[0].lines[-1].get_color())
    axs[0].set(xlabel='Smoothing length $h$',ylabel='RMS curvature-error proxy');axs[0].legend(frameon=False,fontsize=9)
    sigma=5e-4
    axs[1].loglog(h,h**2,label='Bias amplitude $h^2$')
    axs[1].loglog(h,sigma/h**2,label=r'Noise amplitude $\sigma_x/h^2$')
    axs[1].loglog(h,np.sqrt(h**4+sigma**2/h**4),'--',label='Total RMS proxy')
    axs[1].set(xlabel='Smoothing length $h$',ylabel='Error amplitude');axs[1].legend(frameon=False,fontsize=9)
    save(fig,'validation','curvature_noise')


def main():
    for draw in [contact_balance,interface_motion,surface_shear,shape_index,
                 potts_volume,transport_motility,rheology_scales,buckling,curvature_noise]:
        draw()
    POLICY.write_text(json.dumps(MANIFEST,ensure_ascii=False,indent=2)+'\n')
    print('Generated 9 additional SVG/PNG figure pairs; calculation checks passed.')


if __name__=='__main__':main()
