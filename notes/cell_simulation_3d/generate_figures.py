#!/usr/bin/env python3
"""Deterministic teaching figures, not full 3-D tissue simulations.

Run from any directory. Equations and parameter choices are described in body.md.
Only NumPy and Matplotlib are required; no experimental data are used.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
COLORS = ['#247b91', '#d68a42', '#8868ad']
plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'mathtext.fontset': 'cm', 'font.size': 11, 'axes.titlesize': 12,
    'axes.labelsize': 12, 'axes.spines.top': False, 'axes.spines.right': False,
    'axes.prop_cycle': matplotlib.cycler(color=COLORS),
    'svg.hashsalt': 'cell-simulation-3d', 'figure.facecolor': 'white',
    'savefig.facecolor': 'white',
})
POLICY = REPO / 'data/cell_simulation_3d_figures.json'
# Regenerating teaching plots must not remove the original-paper image policies.
MANIFEST = {
    key: value for key, value in (json.loads(POLICY.read_text()) if POLICY.exists() else {}).items()
    if value.get('kind') == 'paper'
}


def save(fig, slug, name, width=820):
    folder = ROOT / slug / 'figures'
    folder.mkdir(parents=True, exist_ok=True)
    fig.savefig(folder / f'{name}.svg', metadata={'Date': None})
    fig.savefig(folder / f'{name}.png', dpi=300)
    w, h = fig.get_size_inches()
    MANIFEST[f'{slug}/figures/{name}'] = {
        'width': width, 'kind': 'calculation',
        'intrinsic_width': round(w * 96), 'intrinsic_height': round(h * 96),
    }
    plt.close(fig)


def pair(titles):
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.25), layout='constrained')
    for ax, title in zip(axs, titles):
        ax.set_title(title, loc='left', pad=12)
        ax.grid(alpha=.16)
    return fig, axs


def sphere_mesh(level=0):
    vertices = np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]], float)
    faces = np.array([[0,2,4],[2,1,4],[1,3,4],[3,0,4],
                      [2,0,5],[1,2,5],[3,1,5],[0,3,5]], int)
    for _ in range(level):
        points = list(vertices)
        cache = {}
        def midpoint(a, b):
            key = tuple(sorted((int(a), int(b))))
            if key not in cache:
                p = (vertices[a] + vertices[b]) / 2
                cache[key] = len(points)
                points.append(p / np.linalg.norm(p))
            return cache[key]
        refined = []
        for a,b,c in faces:
            ab, bc, ca = midpoint(a,b), midpoint(b,c), midpoint(c,a)
            refined.extend([[a,ab,ca],[ab,b,bc],[ca,bc,c],[ab,bc,ca]])
        vertices, faces = np.array(points), np.array(refined)
    return vertices, faces


def geometry(vertices, faces):
    tri = vertices[faces]
    area = np.linalg.norm(np.cross(tri[:,1]-tri[:,0], tri[:,2]-tri[:,0]), axis=1).sum()/2
    volume = np.einsum('ij,ij->i', tri[:,0], np.cross(tri[:,1],tri[:,2])).sum()/6
    return float(area), float(volume)


def setup3d(ax, title):
    ax.set_title(title, loc='left', pad=0)
    ax.set(xlim=(-1.8,1.8), ylim=(-1.5,1.5), zlim=(-1,1.4))
    ax.set_box_aspect((3.6,3,2.4))
    ax.view_init(24,-58)
    ax.set_axis_off()


def representations():
    fig = plt.figure(figsize=(10, 6.5), layout='constrained')
    centers = np.array([[-.8,-.35,0],[.8,-.35,0],[0,1,.05]])
    v, f = sphere_mesh(1)
    for n,title in enumerate(['A  Cell centres','B  Subcellular elements','C  Independent surface meshes']):
        ax=fig.add_subplot(2,3,n+1,projection='3d')
        setup3d(ax,title)
        for i,c in enumerate(centers):
            if n == 0:
                ax.scatter(*c, s=95, color=COLORS[i], depthshade=False)
                # A coarse sphere indicates the assumed radius, not resolved shape.
                ax.add_collection3d(Poly3DCollection((v*.66+c)[f], facecolor=COLORS[i], alpha=.06, linewidth=.2))
            elif n == 1:
                points=np.vstack([np.zeros((1,3)), v*.54])+c
                ax.scatter(*points.T, s=13, color=COLORS[i], depthshade=False)
            else:
                ax.add_collection3d(Poly3DCollection((v*.72+c)[f], facecolor=COLORS[i], edgecolor='#334155', alpha=.32, linewidth=.45))
    ax=fig.add_subplot(2,3,4,projection='3d')
    setup3d(ax,'D  Shared polyhedral faces')
    # Cubes share complete faces; this is deliberately a schematic tessellation.
    corners=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]],float)
    cube_faces=[[0,3,2,1],[4,5,6,7],[0,1,5,4],[1,2,6,5],[2,3,7,6],[3,0,4,7]]
    for i,c in enumerate([[-1,-.7,-.4],[0,-.7,-.4],[-1,.3,-.4]]):
        ax.add_collection3d(Poly3DCollection((corners+np.array(c))[cube_faces], facecolor=COLORS[i], edgecolor='#334155', alpha=.22, linewidth=.7))
    grid=np.linspace(-1.8,1.8,240)
    x,y=np.meshgrid(grid,grid)
    fields=[]
    for c in centers:
        r=np.sqrt((x-c[0])**2+(y-c[1])**2)
        fields.append(1/(1+np.exp((r-.73)/.065)))
    fields=np.array(fields)
    ax=fig.add_subplot(2,3,5)
    rgb=np.ones((*x.shape,3))
    from matplotlib.colors import to_rgb
    for col,phi in zip(COLORS,fields):
        rgb-=phi[...,None]*(1-np.array(to_rgb(col)))
    ax.imshow(np.clip(rgb,0,1),origin='lower',extent=(-1.8,1.8,-1.8,1.8))
    ax.contour(x,y,fields.max(axis=0),levels=[.5],colors='#334155',linewidths=.65)
    ax.set(title='E  Phase fields: central slice',aspect='equal')
    ax.axis('off')
    ax=fig.add_subplot(2,3,6)
    g=np.linspace(-1.8,1.8,33)
    gx,gy=np.meshgrid(g,g)
    distances=np.array([(gx-c[0])**2+(gy-c[1])**2 for c in centers])
    labels=distances.argmin(axis=0)+1
    labels[distances.min(axis=0)>.73**2]=0
    from matplotlib.colors import ListedColormap
    ax.pcolormesh(gx,gy,labels,cmap=ListedColormap(['white']+COLORS),vmin=0,vmax=3,shading='nearest',edgecolors='#e2e8f0',linewidth=.2)
    ax.set(title='F  Voxel labels: central slice',aspect='equal')
    ax.axis('off')
    save(fig,'foundations','representations')


def interface():
    fig,axs=pair(['A  Exact planar interface','B  Interface-energy quadrature'])
    x=np.linspace(-8,8,600)
    for eps in [.5,1,1.8]:
        axs[0].plot(x,1/(1+np.exp(x/eps)),label=fr'$\epsilon={eps}$')
    axs[0].set(xlabel='$x$',ylabel=r'$\phi(x)$',ylim=(-.03,1.03))
    axs[0].legend(frameon=False)
    hs=np.geomspace(.08,3,50)
    for offset in [0,.5]:
        errors=[]
        for h in hs:
            xs=np.arange(-12,12+h/10,h)+offset*h
            p=1/(1+np.exp(xs))
            density=p*p*(1-p)**2 # epsilon=gamma=1; c_W=1/6
            energy=np.trapezoid(density,xs)*6
            errors.append(abs(energy-1))
        axs[1].loglog(hs,errors,label=fr'grid offset $={offset}h$')
    axs[1].set(xlabel=r'$h/\epsilon$',ylabel='Relative surface-energy error')
    axs[1].legend(frameon=False,loc='upper left')
    # Tail truncation produces an error floor at fine resolution.
    xx=np.linspace(-16,16,5000); pp=1/(1+np.exp(xx))
    assert abs(np.trapezoid(6*pp**2*(1-pp)**2,xx)-1)<1e-10
    save(fig,'phase_field','interface')


def mesh_convergence():
    fig=plt.figure(figsize=(10,4.25),layout='constrained')
    ax=fig.add_subplot(1,2,1,projection='3d')
    v,f=sphere_mesh(2)
    ax.add_collection3d(Poly3DCollection(v[f],facecolor=COLORS[0],edgecolor='#334155',alpha=.45,linewidth=.45))
    ax.set(xlim=(-1,1),ylim=(-1,1),zlim=(-1,1),title='A  Refined octahedral sphere')
    ax.set_box_aspect((1,1,1)); ax.view_init(22,35); ax.set_axis_off()
    ax=fig.add_subplot(1,2,2)
    ns,ea,ev=[],[],[]
    for level in range(6):
        v,f=sphere_mesh(level); a,volume=geometry(v,f)
        ns.append(len(f)); ea.append(abs(a/(4*np.pi)-1)); ev.append(abs(volume/(4*np.pi/3)-1))
    assert all(np.diff(ea)<0) and all(np.diff(ev)<0)
    ax.loglog(ns,ea,'o-',label='Surface area')
    ax.loglog(ns,ev,'s-',label='Enclosed volume')
    ax.set(title='B  Error against analytic unit sphere',xlabel='Number of triangles',ylabel='Relative geometric error')
    ax.legend(frameon=False);ax.grid(alpha=.16)
    save(fig,'surface_mesh','mesh_convergence')


def prism():
    fig,axs=pair(['A  Fixed-volume hexagonal prism','B  Equilibrium cell height'])
    c=np.sqrt(8*np.sqrt(3)); g=2.; volume=1.
    a=np.linspace(.15,3.7,500)
    for tension in [.35,.8,1.6]:
        energy=g*a+tension*c*volume/np.sqrt(a)
        star=(tension*c*volume/(2*g))**(2/3)
        axs[0].plot(a,energy,label=fr'$\Gamma_l={tension}$')
        axs[0].plot(star,g*star+tension*c/np.sqrt(star),'o',color=axs[0].lines[-1].get_color())
    axs[0].set(xlabel='Cross-section area $a$',ylabel='Energy $E$')
    axs[0].legend(frameon=False)
    t=np.linspace(.2,2,300)
    axs[1].plot(t,volume/(t*c*volume/(2*g))**(2/3))
    axs[1].set(xlabel=r'Lateral tension $\Gamma_l$',ylabel='Optimal height $h_*$',ylim=(0,4))
    save(fig,'vertex_3d','prism')


def potts():
    fig,axs=pair(['A  Cubic nearest-neighbour anisotropy','B  Uphill copy acceptance'])
    theta=np.linspace(0,90,300);rad=np.deg2rad(theta)
    axs[0].plot(theta,np.abs(np.cos(rad))+np.abs(np.sin(rad)))
    axs[0].axhline(1,color='#64748b',ls=':',label='Axial interface')
    axs[0].set(xlabel='Interface normal angle (degrees)',ylabel='Relative crossing-bond density',xticks=[0,30,45,60,90])
    axs[0].legend(frameon=False)
    delta=np.linspace(0,4,400)
    for temp in [.25,.7,1.5]:
        axs[1].plot(delta,np.exp(-delta/temp),label=fr'$T_{{eff}}={temp}$')
    axs[1].set(xlabel=r'Energy increase $\Delta\mathcal{H}$',ylabel='Acceptance probability',ylim=(0,1.04))
    axs[1].legend(frameon=False)
    save(fig,'lattice_potts','potts')


def particle_potential():
    fig,axs=pair(['A  Morse interaction energy','B  Radial force: positive is repulsive'])
    r=np.linspace(.65,2.2,500);a=3.
    u=np.exp(-2*a*(r-1))-2*np.exp(-a*(r-1))
    force=2*a*(np.exp(-2*a*(r-1))-np.exp(-a*(r-1)))
    axs[0].plot(r,u);axs[1].plot(r,force)
    for ax in axs:
        ax.axvline(1,color='#64748b',ls=':',lw=1)
        ax.axhline(0,color='#64748b',lw=.7)
        ax.set_xlabel('Distance $r/r_0$')
    axs[0].set_ylabel('Energy $U/D$');axs[1].set_ylabel('Force $F r_0/D$')
    # Analytic derivative and central difference must agree away from endpoints.
    assert np.max(abs(force[1:-1]+np.gradient(u,r)[1:-1]))<.01
    save(fig,'particles','particle_potential')


def material_response():
    fig,axs=pair(['A  Isochoric uniaxial deformation','B  Maxwell stress relaxation'])
    stretch=np.linspace(.5,2,400)
    axs[0].plot(stretch,(stretch**2+2/stretch-3)/2)
    axs[0].set(xlabel=r'Axial stretch $\lambda$',ylabel=r'Energy density $\Psi/\mu$')
    t=np.linspace(0,5,400)
    for tau in [.4,1,2]:
        axs[1].plot(t,np.exp(-t/tau),label=fr'$\tau={tau}$')
    axs[1].set(xlabel='Time $t$',ylabel='Normalized modulus $G(t)/G_0$',ylim=(0,1.03))
    axs[1].legend(frameon=False)
    save(fig,'continuum_fluids','material_response')


def lumen_trajectory(r0,dt=.001,end=5):
    t=np.arange(round(end/dt)+1)*dt
    radius=np.zeros_like(t);radius[0]=r0
    def rhs(r): return r**-3-r**-1
    for i in range(len(t)-1):
        r=radius[i]
        k1=rhs(r);k2=rhs(r+dt*k1/2);k3=rhs(r+dt*k2/2);k4=rhs(r+dt*k3)
        radius[i+1]=r+dt*(k1+2*k2+2*k3+k4)/6
    return t,radius


def lumen():
    fig,axs=pair(['A  Fixed-solute spherical lumen','B  Radius relaxes toward $r_*=1$'])
    r=np.linspace(.55,2,400)
    axs[0].plot(r,r**-3,label='Osmotic drive $r^{-3}$')
    axs[0].plot(r,r**-1,label='Laplace resistance $r^{-1}$')
    axs[0].set(xlabel='Radius $r$',ylabel='Scaled pressure',ylim=(0,6.2))
    axs[0].legend(frameon=False)
    for r0 in [.6,1.4,1.8]:
        t,radius=lumen_trajectory(r0)
        axs[1].plot(t,radius,label=fr'$r(0)={r0}$')
        assert np.isfinite(radius).all() and np.all(radius>0)
        assert abs(radius[-1]-1)<abs(r0-1)
        t2,r2=lumen_trajectory(r0,dt=.0005)
        assert abs(radius[-1]-r2[-1])<1e-7
    axs[1].axhline(1,color='#64748b',ls=':',lw=1)
    axs[1].set(xlabel=r'Scaled time $\tau$',ylabel='Scaled radius $r$',ylim=(.5,1.9))
    axs[1].legend(frameon=False)
    save(fig,'morphogenesis','lumen')


def inference():
    fig,axs=pair(['A  Shape constrains a ratio','B  Add an independent tension measurement'])
    g,w=np.meshgrid(np.linspace(.35,2,350),np.linspace(.1,1.9,350))
    ratio_loss=((w/g-.8)/.08)**2
    for ax,loss in zip(axs,[ratio_loss,ratio_loss+((g-1)/.15)**2]):
        contours=ax.contourf(g,w,loss,levels=[0,1,4,9,16,25],cmap='Blues_r',extend='max')
        ax.plot([.35,2],[.28,1.6],color=COLORS[1],ls='--',lw=1.2,label=r'$W/\gamma=0.8$')
        ax.plot(1,.8,'o',color=COLORS[1],ms=5)
        ax.set(xlabel=r'Surface tension $\gamma$',ylabel='Adhesion $W$')
        ax.grid(False)
    fig.colorbar(contours,ax=axs,label='Synthetic weighted loss',shrink=.85)
    save(fig,'validation','inference')


def main():
    for draw in [representations,interface,mesh_convergence,prism,potts,
                 particle_potential,material_response,lumen,inference]:
        draw()
    target=REPO/'data/cell_simulation_3d_figures.json'
    target.write_text(json.dumps(MANIFEST,ensure_ascii=False,indent=2)+'\n')
    count = sum(value['kind'] == 'calculation' for value in MANIFEST.values())
    print(f'Generated {count} teaching SVG/PNG pairs; paper policies preserved; numerical checks passed.')


if __name__=='__main__':
    main()
