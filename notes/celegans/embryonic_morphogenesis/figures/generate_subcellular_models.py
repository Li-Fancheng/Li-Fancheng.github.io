"""Dimensionless subcellular examples for note 07; no embryo data are fitted.

Dependencies: NumPy, SciPy, Matplotlib. Run from any working directory.
The five outputs solve stated field equations or exact kinematic examples.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Ellipse

OUT=Path(__file__).resolve().parent
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,
                     'svg.fonttype':'none','font.family':'DejaVu Sans','savefig.facecolor':'white'})
TEAL='#27858b';ROSE='#b65376';INK='#304857';GOLD='#bb8a23'
report={}
def save(fig,name):
    fig.savefig(OUT/(name+'.svg'),bbox_inches='tight')
    fig.savefig(OUT/(name+'.png'),dpi=170,bbox_inches='tight')
    plt.close(fig)


def screened_flow():
    """Solve (Gamma-eta Laplacian)v=grad(tau) on a periodic square."""
    n=128;L=12.;dx=L/n;x=np.arange(n)*dx-L/2;X,Y=np.meshgrid(x,x)
    tau=np.exp(-(X*X+Y*Y)/(2*.8**2))
    k=2*np.pi*np.fft.fftfreq(n,d=dx);KX,KY=np.meshgrid(k,k);K2=KX*KX+KY*KY
    F=np.fft.fft2(tau);eta=1.;fig,axs=plt.subplots(1,3,figsize=(13,4.4),layout='constrained')
    cs=axs[0].contourf(X,Y,tau,levels=16,cmap='magma');fig.colorbar(cs,ax=axs[0],shrink=.7,label='Active tensile stress')
    axs[0].set_title('Localized contractility, fixed cell outline')
    vals={};flows=[]
    for gamma in [4.,.16]:
        vx=np.fft.ifft2(1j*KX*F/(gamma+eta*K2)).real
        vy=np.fft.ifft2(1j*KY*F/(gamma+eta*K2)).real
        flows.append((vx,vy))
    vmax=max(np.hypot(vx,vy).max() for vx,vy in flows)
    for ax,gamma,(vx,vy) in zip(axs[1:],[4.,.16],flows):
        residual=np.fft.fft2(vx)*(gamma+eta*K2)-1j*KX*F
        error=float(np.max(np.abs(residual))/np.max(np.abs(KX*F)))
        assert error<1e-10
        speed=np.hypot(vx,vy);cf=ax.contourf(X,Y,speed,levels=np.linspace(0,vmax,16),cmap='Blues')
        ax.streamplot(x,x,vx,vy,density=.7,color=TEAL,linewidth=.9,arrowsize=.9)
        ax.set_title(f'Friction {gamma:g}, screening length {np.sqrt(eta/gamma):g}')
        vals[str(gamma)]={'screening_length':float(np.sqrt(eta/gamma)),'spectral_relative_residual':error,'peak_speed':float(speed.max())}
    for ax in axs:ax.set_aspect('equal');ax.set_xlim(-5,5);ax.set_ylim(-5,5);ax.set_xlabel('Surface coordinate x');ax.set_ylabel('Surface coordinate y')
    fig.colorbar(cf,ax=axs[1:].tolist(),orientation='horizontal',shrink=.7,pad=.13,label='Speed (common scale)')
    fig.suptitle('Surface flow can change while the cell center stays fixed',color=INK)
    save(fig,'32_cortical_flow_screening');report['screening']={'L':L,'n':n,'eta':eta,'patch_width':.8,'cases':vals}


def transport_exchange():
    """Exact constant-speed advection-diffusion-decay of a labelled Gaussian."""
    x=np.linspace(-5,7,361);t=np.linspace(0,8,201);X,T=np.meshgrid(x,t)
    sig=.5;D=.025;v=.6;koff=.25
    var=sig**2+2*D*T
    mobile=sig/np.sqrt(var)*np.exp(-(X+2-v*T)**2/(2*var))*np.exp(-koff*T)
    # A moving recruitment wave creates the same total concentration pattern,
    # while each molecule stays in place until it unbinds. S may have either sign.
    ct0=np.exp(-(x+2)**2/(2*sig**2))
    dx=x[1]-x[0]
    mass=np.trapezoid(mobile,x,axis=1)
    exact=np.sqrt(2*np.pi)*sig*np.exp(-koff*t)
    error=float(np.max(np.abs(mass-exact))/exact[0]);assert error<.001
    fig,axs=plt.subplots(2,2,figsize=(12,7),layout='constrained')
    for ax,title in zip(axs[0],['Material transport: molecules move','Recruitment wave: molecules exchange locally']):
        im=ax.pcolormesh(x,t,mobile,shading='auto',cmap='magma',vmin=0,vmax=1)
        ax.set_title(title);ax.set_xlabel('Position x');ax.set_ylabel('Time t')
    for x0 in np.linspace(-2.8,-1.2,7):axs[0,0].plot(x0+v*t[t<=4],t[t<=4],color='cyan',lw=.8,alpha=.75)
    for x0 in [-2,0,2]:
        t0=max(0,(x0+2)/v);tt=np.linspace(t0,min(t0+1.4,8),30);axs[0,1].plot(np.full_like(tt,x0),tt,color='cyan',lw=1.3)
    for ti in [0,2,4,8]:
        j=np.argmin(abs(t-ti));axs[1,0].plot(x,mobile[j],label=f't={ti}')
    axs[1,0].set(title='Identical intensity profiles in the two constructions',xlabel='Position x',ylabel='Labelled intensity');axs[1,0].legend(ncol=2)
    axs[1,1].axis('off');axs[1,1].text(.02,.9,'A travelling bright pattern does not identify\nthe trajectories of its constituent molecules.\n\nCyan lines: illustrative single-molecule paths.\nLeft: x(t) = x0 + vt.\nRight: fixed x during a membrane residence.\n\nSingle-particle tracking and local assembly rates\nseparate these two possibilities.',va='top',fontsize=12,color=INK,linespacing=1.5)
    fig.colorbar(im,ax=axs[0].tolist(),shrink=.8,label='Fluorescence pattern')
    save(fig,'33_transport_and_exchange');report['transport']={'v':v,'D':D,'koff':koff,'sigma0':sig,'relative_integrated_mass_error':error,'note':'Right panel is a prescribed recruitment-wave construction, not an independently fitted reaction model.'}


def excitable_strip():
    """FitzHugh-Nagumo normal form, no-flux 1D strip; u,h are scaled variables."""
    a=.7;b=.8;eps=.08;D=.6
    ur=brentq(lambda u:u-u**3/3-(u+a)/b,-2,-1);hr=(ur+a)/b
    x=np.linspace(-20,20,401);dx=x[1]-x[0];dt=.004;steps=15000
    u=np.full_like(x,ur);h=np.full_like(x,hr);frames=[];times=[]
    for n in range(steps+1):
        time=n*dt
        if n%50==0:frames.append((u.copy(),h.copy()));times.append(time)
        lap=(np.r_[u[1],u[:-1]]-2*u+np.r_[u[1:],u[-2]])/dx**2
        I=(1.5*np.exp(-((x+12)/1.3)**2)) if 1<=time<2 else 0
        du=u-u**3/3-h+D*lap+I;dh=eps*(u+a-b*h)
        u=u+dt*du;h=h+dt*dh
    frames=np.array(frames);times=np.array(times);assert np.isfinite(frames).all()
    assert frames[:,0,:].max()>1.4
    eigen=np.linalg.eigvals([[1-ur*ur,-1],[eps,-eps*b]])
    assert np.max(eigen.real)<0
    fig,axs=plt.subplots(2,2,figsize=(12,7.2),layout='constrained')
    for j,title,cmap in [(0,'Fast excitation u','magma'),(1,'Delayed recovery h','viridis')]:
        im=axs[0,j].pcolormesh(x,times,frames[:,j,:],shading='auto',cmap=cmap)
        axs[0,j].set(xlabel='Position x',ylabel='Time t',title=title);fig.colorbar(im,ax=axs[0,j],shrink=.8)
    for ti in [0,5,12,22]:
        jj=np.argmin(abs(times-ti));axs[1,0].plot(x,frames[jj,0],label=f't={ti}')
    axs[1,0].set(xlabel='Position x',ylabel='u',title='Pulse propagation after a local transient input');axs[1,0].legend(ncol=2)
    j=np.argmin(abs(x+5));axs[1,1].plot(times,frames[:,0,j],color=ROSE,label='u');axs[1,1].plot(times,frames[:,1,j],color=TEAL,label='h')
    axs[1,1].set(xlabel='Time t',title='Local excitation precedes recovery (x=-5)');axs[1,1].legend()
    save(fig,'34_excitable_cortex_strip');report['excitation']={'a':a,'b':b,'epsilon':eps,'D':D,'dt':dt,'dx':dx,'rest':[ur,hr],'rest_eigenvalues_real':eigen.real.tolist(),'input':'1.5 exp(-((x+12)/1.3)^2), 1<=t<2','note':'Normal form; scaled u,h are not molecular concentrations.'}


def affine_alignment():
    """Exact push-forward of initially isotropic material line segments."""
    rng=np.random.default_rng(9);n=180
    ctr=rng.uniform([-1.8,-1.4],[1.8,1.4],(n,2));angles=rng.uniform(0,np.pi,n);p=np.c_[np.cos(angles),np.sin(angles)]
    fig,axs=plt.subplots(1,4,figsize=(13.8,4.4),layout='constrained')
    for ax,strain in zip(axs[:3],[0,.8,1.6]):
        F=np.diag([np.exp(-strain),1.]);c=ctr@F.T;v=p@F.T
        seg=np.stack([c-.16*v,c+.16*v],axis=1);ax.add_collection(LineCollection(seg,colors=TEAL,lw=1.3))
        ax.add_patch(plt.Rectangle((-2*np.exp(-strain),-1.6),4*np.exp(-strain),3.2,fill=False,color=INK,lw=1))
        ax.set_title(f'Compression strain {strain:g}');ax.set_xlim(-2.2,2.2);ax.set_ylim(-1.8,1.8);ax.set_aspect('equal');ax.axis('off')
    strain=np.linspace(0,2,100);theta=np.linspace(0,np.pi,20000,endpoint=False)
    measured=[]
    for s in strain:
        vx=np.exp(-s)*np.cos(theta);vy=np.sin(theta);measured.append(np.mean((vy*vy-vx*vx)/(vx*vx+vy*vy)))
    exact=np.tanh(strain/2);err=float(np.max(np.abs(np.array(measured)-exact)));assert err<1e-9
    axs[3].plot(strain,exact,color=ROSE);axs[3].set(xlabel='Compression strain',ylabel=r'$S_y=\langle p_y^2-p_x^2\rangle$',title='Alignment follows deformation',ylim=(0,.85))
    fig.suptitle('The same initial rods align along y when compressed along x',color=INK)
    save(fig,'35_flow_alignment_geometry');report['alignment']={'F':'diag(exp(-strain),1)','order':'Sy=tanh(strain/2)','angular_quadrature_max_error':err,'note':'Exact affine kinematics; excludes turnover, active alignment, and bending.'}


def bulk_recirculation():
    """Exact 2D Stokes solution in an ellipse, tangent moving boundary."""
    a=2.;b=1.;C=1.;eta=1.
    x=np.linspace(-a,a,241);y=np.linspace(-b,b,161);X,Y=np.meshgrid(x,y);inside=X*X/a**2+Y*Y/b**2<=1
    psi=C*Y*(1-X*X/a**2-Y*Y/b**2)
    U=C*(1-X*X/a**2-3*Y*Y/b**2);V=2*C*X*Y/a**2
    pressure=-2*eta*C*(1/a**2+3/b**2)*X
    theta=np.linspace(0,2*np.pi,500);xb=a*np.cos(theta);yb=b*np.sin(theta)
    ub=C*(1-xb*xb/a**2-3*yb*yb/b**2);vb=2*C*xb*yb/a**2
    normal_flux=ub*xb/a**2+vb*yb/b**2;assert np.max(abs(normal_flux))<1e-12
    fig,axs=plt.subplots(1,3,figsize=(13,4.4),layout='constrained')
    u=np.ma.masked_where(~inside,U);v=np.ma.masked_where(~inside,V)
    im=axs[0].contourf(X,Y,u,levels=np.linspace(-2,1,16),cmap='RdBu_r');axs[0].streamplot(x,y,u,v,color=INK,density=.85,linewidth=.8)
    axs[0].set_title('Peripheral drive, central return');fig.colorbar(im,ax=axs[0],orientation='horizontal',shrink=.7,pad=.18,label=r'$u_x$')
    axs[1].contour(X,Y,np.ma.masked_where(~inside,psi),levels=np.linspace(-.38,.38,16),colors=TEAL);axs[1].set_title('Closed streamlines inside one cell')
    axs[2].plot(C*(1-3*y*y/b**2),y,color=ROSE);axs[2].axvline(0,color='#abb5bd',lw=1);axs[2].set(xlabel=r'$u_x(0,y)$',ylabel='Transverse position y',title='Center and periphery move oppositely')
    for ax in axs[:2]:ax.add_patch(Ellipse((0,0),2*a,2*b,fill=False,color=INK));ax.set_aspect('equal');ax.set_xlabel('x');ax.set_ylabel('y')
    save(fig,'36_cytoplasmic_return_flow');report['bulk_flow']={'a':a,'b':b,'C':C,'eta':eta,'streamfunction':'C*y*(1-x^2/a^2-y^2/b^2)','pressure':'-2*eta*C*(1/a^2+3/b^2)*x','max_boundary_normal_flux':float(np.max(abs(normal_flux))),'divergence':'-2Cx/a^2+2Cx/a^2=0','note':'2D analytic Stokes example; not the fitted 3D Niwayama et al. model.'}

if __name__=='__main__':
    screened_flow();transport_exchange();excitable_strip();affine_alignment();bulk_recirculation()
    (OUT/'subcellular_model_parameters.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
