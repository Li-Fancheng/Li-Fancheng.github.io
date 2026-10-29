"""Original, reproducible mechanism comparisons for the introductory course.

Use requirements.txt. All simulations here use dimensionless coordinates and
fresh renewal preparations unless explicitly stated. No sample MSD is computed
for the stable Levy flight, whose theoretical second moment is infinite.
"""
import json
import sys
sys.dont_write_bytecode = True
import numpy as np
from scipy.integrate import quad, simpson
from scipy.linalg import cholesky, toeplitz
from scipy.special import gamma, gammaincc
from scipy.stats import levy_stable
import generate_figures as style
from matplotlib import pyplot as plt
from matplotlib.patches import FancyArrowPatch

ROOT = style.ROOT
BLUE, ORANGE, GREEN, INK = style.BLUE, style.ORANGE, style.GREEN, style.INK


def brownian_overview():
    """The same free dynamics, with an additional endpoint condition in C."""
    rng = np.random.default_rng(2026102914)
    diffusion, duration, endpoint, steps = .3, 2., 1.2, 400
    t = np.linspace(0, duration, steps + 1)
    w = np.c_[np.zeros(5), np.cumsum(
        rng.normal(scale=np.sqrt(duration / steps), size=(5, steps)), axis=1)]
    free = np.sqrt(2 * diffusion) * w
    bridges = free - t[None, :] / duration * free[:, -1:] + endpoint * t / duration
    fig, ax = plt.subplots(1, 3, figsize=(12.6, 4.2), layout="constrained")
    for path in free:
        ax[0].plot(t, path, color=BLUE, alpha=.55, lw=.85)
    ax[0].axhline(0, color=INK, ls="--", lw=1.2)
    ax[0].set(title="A  Repeated free trajectories", xlabel="time t", ylabel="position x")
    x = np.linspace(-3.8, 3.8, 500)
    for time, color in zip([.15, .6, 2.], [GREEN, ORANGE, BLUE]):
        density = style.gaussian(x, 0., 2 * diffusion * time)
        ax[1].plot(x, density, color=color, label=f"t = {time:g}")
    ax[1].set(title="B  Unconditioned densities", xlabel="position x", ylabel="density p(x,t)", ylim=(0, 1.4))
    ax[1].legend(fontsize=10)
    for path in bridges:
        ax[2].plot(t, path, color=ORANGE, alpha=.55, lw=.85)
    ax[2].plot(t, endpoint * t / duration, color=INK, ls="--", lw=1.2)
    ax[2].scatter([0, duration], [0, endpoint], color=INK, s=25, zorder=5)
    ax[2].set(title="C  Conditioned on the endpoint", xlabel="time t", ylabel="position x")
    # Identical trajectory axes make the additional endpoint constraint visible.
    extent = max(np.max(np.abs(free)), np.max(np.abs(bridges))) + .15
    for panel in [ax[0], ax[2]]:
        panel.set(xlim=(0, duration), ylim=(-extent, extent))
    style.save(fig, "14_brownian_views")


def inertial_crossover():
    u = np.geomspace(1e-3, 100, 800)
    exact = 2 * (u + np.expm1(-u))
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.9), layout="constrained")
    ax[0].loglog(u, exact, color=BLUE, label="exact MSD")
    short, long = u < 0.4, u > 3
    ax[0].loglog(u[short], u[short]**2, color=ORANGE, ls="--", label="ballistic: slope 2")
    ax[0].loglog(u[long], 2*u[long], color=GREEN, ls="--", label="diffusive: slope 1")
    ax[0].set(xlabel=r"$t/\tau_m$", ylabel=r"$\langle\Delta X^2\rangle/(D\tau_m)$", title="A  One particle, two time regimes")
    ax[0].legend(fontsize=9)
    lag = np.linspace(0, 6, 300)
    ax[1].plot(lag, np.exp(-lag), color=BLUE)
    ax[1].axvline(1, color=ORANGE, ls=":", label="velocity memory time")
    ax[1].set(xlabel=r"$t/\tau_m$", ylabel=r"$C_v(t)/C_v(0)$", title="B  Velocity direction loses its memory", ylim=(0,1.05))
    ax[1].legend(fontsize=9)
    style.save(fig, "09_inertial_crossover")


def ring_current():
    length, speed, diffusion = 6., 0.8, 0.4
    fig = plt.figure(figsize=(10,4.2), layout="constrained")
    gs = fig.add_gridspec(2,2)
    circle = fig.add_subplot(gs[:,0]); density = fig.add_subplot(gs[0,1]); current = fig.add_subplot(gs[1,1])
    theta = np.linspace(0,2*np.pi,400)
    circle.plot(np.cos(theta),np.sin(theta),color=BLUE,alpha=.5)
    for angle in np.linspace(0,2*np.pi,8,endpoint=False):
        circle.add_patch(FancyArrowPatch((np.cos(angle),np.sin(angle)),(np.cos(angle+.4),np.sin(angle+.4)),connectionstyle="arc3,rad=0.15",arrowstyle="-|>",mutation_scale=15,color=ORANGE,lw=2))
    circle.text(0,.1,r"$p_{\rm st}=1/L$",ha="center")
    circle.text(0,-.2,r"$J_{\rm st}=v/L>0$",ha="center",color=ORANGE)
    circle.set(title="A  Constant occupancy, circulating current",aspect="equal",xlim=(-1.4,1.4),ylim=(-1.25,1.25));circle.axis("off")
    x=np.linspace(0,length,100)
    density.plot(x,np.full_like(x,1/length),color=BLUE)
    density.set(xlim=(0,length),ylim=(0,.23),ylabel="density",title="B  Unchanged density and steady flux")
    current.plot(x,np.full_like(x,speed/length),color=ORANGE)
    current.set(xlim=(0,length),ylim=(0,.2),xlabel="arc length x",ylabel="probability / time")
    style.save(fig,"10_ring_current")
    return {"normalization_error":abs(length*(1/length)-1),"stationary_current":speed/length,"diffusion":diffusion}


def levy_comparison():
    rng=np.random.default_rng(2026102911)
    alpha, duration, count, speed, scale=1.3,20.,1000,1.,.3
    times=np.linspace(0,duration,count+1);dt=duration/count
    increments=levy_stable.rvs(alpha,0,scale=dt**(1/alpha),size=count,random_state=rng)
    flight=np.r_[0,np.cumsum(increments)]
    wt,wx=[0.],[0.]
    while wt[-1]<duration:
        run=scale*rng.pareto(alpha)
        end=min(duration,wt[-1]+run)
        wx.append(wx[-1]+rng.choice([-1,1])*speed*(end-wt[-1]));wt.append(end)
    wt,wx=np.array(wt),np.array(wx)
    fig,ax=plt.subplots(1,2,figsize=(10,4.1),layout="constrained")
    ax[0].step(times,flight,where="post",color=BLUE,lw=1.1)
    ax[0].set(title="A  Stable flight: instantaneous jumps",xlabel="time",ylabel="position")
    ax[1].plot(wt,wx,color=ORANGE,lw=1.5)
    ax[1].plot([0,duration],[0,speed*duration],color=INK,ls=":",lw=1)
    ax[1].plot([0,duration],[0,-speed*duration],color=INK,ls=":",lw=1)
    ax[1].set(title="B  Constant-speed walk: long flights",xlabel="time",ylabel="position")
    style.save(fig,"11_levy_flight_walk")
    front_error=float(np.max(np.maximum(np.abs(wx)-speed*wt,0)))
    speed_error=float(np.max(np.abs(np.abs(np.diff(wx)/np.diff(wt))-speed)))
    # Validate the stable scale with a bounded characteristic function, not a moment.
    sample=levy_stable.rvs(alpha,0,scale=.2**(1/alpha),size=60000,random_state=rng)
    cf_error=float(abs(np.mean(np.exp(.7j*sample))-np.exp(-.2*.7**alpha)))
    assert front_error<1e-10 and speed_error<1e-9 and cf_error<.02
    return {"walk_front_max_excess":front_error,"walk_speed_max_error":speed_error,"stable_characteristic_function_error":cf_error}


def ctrw_fbm():
    rng=np.random.default_rng(2026102912)
    beta,hurst,tau0,duration,steps,samples=.6,.3,.05,20.,600,1024
    t=np.linspace(0,duration,steps+1);dt=duration/steps
    lag=np.arange(steps,dtype=float)
    rho=(np.abs(lag+1)**(2*hurst)+np.abs(lag-1)**(2*hurst)-2*lag**(2*hurst))/2
    covariance=toeplitz(rho)*dt**(2*hurst)
    factor=np.tril(cholesky(covariance,lower=True))  # only the lower factor is used
    fgn=rng.standard_normal((samples,steps))@factor.T
    fbm=np.c_[np.zeros(samples),np.cumsum(fgn,axis=1)]
    paths=np.zeros((samples,steps+1));counts=np.zeros_like(paths)
    for sample in range(samples):
        event_times=[];positions=[0.];now=0.
        while True:
            now+=tau0*rng.pareto(beta)
            if now>duration:break
            event_times.append(now);positions.append(positions[-1]+rng.choice([-1.,1.]))
        n=np.searchsorted(event_times,t,side="right")
        counts[sample]=n;paths[sample]=np.asarray(positions)[n]
    ctrw_inc=np.diff(paths,axis=1)
    msd_c=np.mean(paths**2,axis=0);msd_f=np.mean(fbm**2,axis=0)
    lags=np.arange(1,9)
    corr_c=np.array([np.mean(ctrw_inc[:,:-k]*ctrw_inc[:,k:])/np.mean(ctrw_inc**2) for k in lags])
    corr_f=np.array([np.mean(fgn[:,:-k]*fgn[:,k:])/np.mean(fgn**2) for k in lags])
    theory=rho[lags]
    fig,ax=plt.subplots(2,2,figsize=(10,7.4),layout="constrained")
    ax[0,0].step(t,paths[0],where="post",color=BLUE,lw=1.1)
    ax[0,0].set(title=r"A  CTRW: waits and jumps, $\beta=0.6$",xlabel="time",ylabel="position")
    ax[0,1].plot(t,fbm[0],color=ORANGE,lw=.9)
    ax[0,1].set(title=r"B  fBm: correlated Gaussian steps, $H=0.3$",xlabel="time",ylabel="position")
    selected=t>=.2
    ax[1,0].loglog(t[selected],msd_c[selected]/msd_c[-1],color=BLUE,label="CTRW")
    ax[1,0].loglog(t[selected],msd_f[selected]/msd_f[-1],color=ORANGE,label="fBm")
    ax[1,0].loglog(t[selected],(t[selected]/duration)**beta,color=INK,ls="--",label=r"$t^{0.6}$ guide")
    ax[1,0].set(title="C  Similar growth can hide different clocks",xlabel="time",ylabel="MSD / terminal MSD");ax[1,0].legend(fontsize=9)
    ax[1,1].axhline(0,color=INK,lw=1)
    ax[1,1].plot(lags,corr_c,"o-",color=BLUE,label="CTRW, sampled")
    ax[1,1].plot(lags,corr_f,"o",color=ORANGE,label="fBm, sampled")
    ax[1,1].plot(lags,theory,color=ORANGE,ls="--",label="fBm, exact")
    ax[1,1].set(title="D  Anti-persistence separates the mechanisms",xlabel="increment lag",ylabel="normalized correlation");ax[1,1].legend(fontsize=9)
    style.save(fig,"12_ctrw_fbm")
    delay=6
    def time_average(x,end):return float(np.mean((x[:,delay:end+1]-x[:,:end+1-delay])**2))
    results={"fbm_covariance_factor_relative_error":float(np.max(abs(factor@factor.T-covariance))/covariance[0,0]),
             "fbm_lag_one_correlation":float(corr_f[0]),"fbm_lag_one_exact":float(theory[0]),
             "ctrw_lag_one_correlation":float(corr_c[0]),
             "ctrw_mean_steps_first_5":float(np.mean(counts[:,150])),
             "ctrw_mean_steps_age_10_to_15":float(np.mean(counts[:,450]-counts[:,300])),
             "fbm_msd_first_5":float(np.mean(fbm[:,150]**2)),
             "fbm_msd_age_10_to_15":float(np.mean((fbm[:,450]-fbm[:,300])**2)),
             "time_average_delay":delay*dt,
             "ctrw_time_average_T10":time_average(paths,300),"ctrw_time_average_T20":time_average(paths,600),
             "fbm_time_average_T10":time_average(fbm,300),"fbm_time_average_T20":time_average(fbm,600),
             "fbm_time_average_exact":(delay*dt)**beta}
    assert results["fbm_covariance_factor_relative_error"]<1e-12
    assert abs(corr_f[0]-theory[0])<.02
    assert abs(results["fbm_time_average_T20"]/results["fbm_time_average_exact"]-1)<.05
    return results


def history_observables():
    rng=np.random.default_rng(2026102913)
    diffusion,duration,steps,samples=.2,2.,400,4000
    t=np.linspace(0,duration,steps+1);dt=duration/steps
    w=np.c_[np.zeros(samples),np.cumsum(rng.normal(scale=np.sqrt(dt),size=(samples,steps)),axis=1)]
    bridge=np.sqrt(2*diffusion)*(w-t[None,:]/duration*w[:,-1:])
    exact=2*diffusion*t*(1-t/duration)
    variance=np.var(bridge,axis=0)
    fig,ax=plt.subplots(1,2,figsize=(10,4),layout="constrained")
    ax[0].plot(t,variance,color=BLUE,label="4000 sampled bridges")
    ax[0].plot(t,exact,color=ORANGE,ls="--",label=r"$2Dt(1-t/T)$")
    ax[0].set(title="A  Endpoints suppress bridge fluctuations",xlabel="time",ylabel="conditional variance");ax[0].legend(fontsize=9)
    sample_index=next(i for i, row in enumerate(w) if row.min() < -.1 and row.max() > .1)
    x=np.sqrt(2*diffusion)*w[sample_index]
    ax[1].plot(t,x,color=INK,lw=.9)
    ax[1].fill_between(t,x,0,where=x>=0,color=BLUE,alpha=.3)
    ax[1].fill_between(t,x,0,where=x<0,color=ORANGE,alpha=.3)
    ax[1].axhline(0,color=INK,lw=.6)
    ax[1].set(title="B  A signed area records the whole history",xlabel="time",ylabel="free Brownian position")
    style.save(fig,"13_history_observables")
    grid_area_variance=2*diffusion*dt*np.sum((duration-(np.arange(steps)+.5)*dt)**2)
    areas=np.trapezoid(np.sqrt(2*diffusion)*w,x=t,axis=1)
    return {"bridge_sample_variance_max_absolute_error":float(np.max(abs(variance-exact))),
            "bridge_endpoint_max_error":float(np.max(abs(bridge[:,[0,-1]]))),
            "brownian_area_grid_variance":float(grid_area_variance),
            "brownian_area_continuous_variance":2*diffusion*duration**3/3,
            "brownian_area_sample_variance":float(np.var(areas))}


def analytical_checks():
    errors=[];run_errors=[]
    for time in [.001,.1,1.,10.]:
        numeric=2*quad(lambda s:(time-s)*np.exp(-s),0,time,epsabs=1e-12)[0]
        errors.append(abs(numeric-2*(time+np.expm1(-time))))
        speed,rate=1.7,.6
        numeric=2*speed**2*quad(lambda s:(time-s)*np.exp(-2*rate*s),0,time,epsabs=1e-12)[0]
        prediction=speed**2*time/rate-speed**2/(2*rate**2)*(-np.expm1(-2*rate*time))
        run_errors.append(abs(numeric-prediction))
    # Direct density-weight difference under reversal, using an arbitrary record.
    path=np.array([.3,.1,.8,.2,1.7]);epsilon=.2;drift=.8;diffusion=.4
    dx=np.diff(path)
    log_ratio=np.sum(((-dx-drift*epsilon)**2-(dx-drift*epsilon)**2)/(4*diffusion*epsilon))
    beta,tau0,s=.6,.05,1e-6
    z=s*tau0
    one_minus_psi=np.exp(z)*z**beta*gamma(1-beta)*gammaincc(1-beta,z)
    asymptotic=gamma(1-beta)*tau0**beta*s**beta
    results={"inertial_correlation_integral_max_error":max(errors),"persistent_correlation_integral_max_error":max(run_errors),
             "ring_time_reversal_log_ratio_error":float(abs(log_ratio-drift/diffusion*(path[-1]-path[0]))),
             "lomax_waiting_laplace_small_s_relative_error":float(abs(one_minus_psi/asymptotic-1))}
    assert max(errors)<1e-10 and max(run_errors)<1e-10
    assert results["ring_time_reversal_log_ratio_error"]<1e-12
    assert results["lomax_waiting_laplace_small_s_relative_error"]<.01
    return results


if __name__=="__main__":
    brownian_overview()
    inertial_crossover()
    checks={"ring":ring_current(),"levy":levy_comparison(),"ctrw_fbm":ctrw_fbm(),
            "histories":history_observables(),"analytical":analytical_checks()}
    (ROOT/"intro_numerical_checks.json").write_text(json.dumps(checks,indent=2)+"\n")
    print(json.dumps(checks,indent=2))
