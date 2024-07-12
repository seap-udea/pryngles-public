import pryngles as pr
from pryngles import Consts
extension = 'cpixx'
print(f"Pryngles version: {pr.version}")

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
import multiprocessing as mp
import math
import time
import cv2 as cv
import pickle
from tqdm import tqdm
from joblib import Parallel,delayed,cpu_count

tmp_dir="tmp/"
fig_dir="figures/"

def save_figure(fig, fig_name, fig_dir='/tmp', dpi=600, formats=['png','svg','eps'],**options):
    for form in formats:
        fig_file = fig_dir + fig_name + '.' + form
        print(f"Saving {fig_file}...")
        fig.savefig(fig_file,format=form,dpi=dpi,bbox_inches='tight',**options)

def pool_handler(loc_num, params, func, multiplier,max_process=None):
    
    if max_process is None:
        max_process = math.floor(mp.cpu_count()*multiplier)
    
    # Max allowed thread count on the server is 14
    if max_process>14: max_process=14
    
    if loc_num > max_process:
        run_num = 0
        while loc_num > 0:
            lbound = run_num*max_process
            if loc_num > max_process: ubound = (run_num+1)*max_process
            else: ubound = lbound + loc_num
            processes = [mp.Process(target=func, args=(params[i]))
                        for i in range(lbound, ubound)]
            for p in processes:
                time.sleep(0.1)
                p.start()
            for p in processes:
                time.sleep(0.1)
                p.join()
            loc_num = loc_num - max_process
            run_num+=1
    else:
        processes = [mp.Process(target=func, args=(params[i]))
                    for i in range(0, loc_num)]
        for p in processes:
            time.sleep(0.1)
            p.start()
        for p in processes:
            time.sleep(0.1)
            p.join()

def simulation(iorb: float):
    
    # Constant parameters
    name = f"Orbit_incl80_fix"
    value = iorb

    #Parameters
    fou_file_ring = tmp_dir + "fou_ring_0_4_0_8.dat"
    fou_file_planet = tmp_dir + "fou_bmsca10.0_asurf0.5.dat"

    # Though ring is transparent we must provide options
    tau_ring = 0.0 # Ring is transparent
    ring_ri = 2.0
    ring_re = 2.25
    gamma = 80 # Ring inclination
    alphar = 0 # Ring roll

    # Scattering options
    interp = "spline"
    reference_plane = "Detector"

    # Number of spangles
    Ns = 30
    Nb = 0
    Np = 12000
    Nr = 100

    # Other parameters
    a = 1
    r_s = Consts.rsun
    r_p = Consts.rsaturn
    e = 0.0
    theta_end = 360
    n_theta = 361
    allow_non_uni = False
    normalize = True
    lambq_offset = 0.0
    limb_coeffs = [0.65]
    interpr = "spline"
    
    # Make sure all output is printed
    pr.Verbose.VERBOSITY=pr.VERB_NONE

    # Calculate starting position of observer and star
    ir, beta_obs, lamb_obs, lamb_star = pr.Util.calcStartingPosition(iorb,gamma,alphar)
    print(f"Initial conditions in the planet-centered system are:")
    print(f"\tRing inclination w.r.t. ecliptic: {ir*Consts.rad:.1f} deg")
    print(f"\tInclination of the observer w.r.t ecliptic: {beta_obs*Consts.rad:.1f} deg")
    print(f"\tInitial longitude of the observer: {lamb_obs*Consts.rad:.1f} deg")
    print(f"\tInitial longitude of the star: {lamb_star*Consts.rad:.1f} deg")

    # Intialize system
    sys = pr.System()
    s=sys.add(kind="Star",physics=dict(radius=r_s/sys.ul),optics=dict(limb_coeffs=limb_coeffs))
    p=sys.add(kind="Planet", primary=s,
                   radius=r_p/sys.ul,
                   orbit=dict(a=a, e=e),
                   physics=dict(radius=r_p/sys.ul),
                   optics=dict(nspangles=Np))
    r=sys.add(kind="Ring", primary=p,
                   physics=dict(fi=ring_ri, fe=ring_re, i=ir),
                   optics=dict(nspangles=Nr))
    RP=sys.ensamble_system(extension=extension,
                           fname_planet=fou_file_planet,
                           fname_ring=fou_file_ring)
    RP.behavior["interp_method_ring"] = interpr
    RP.behavior["allow_non_uni"] = allow_non_uni
    RP.reference_plane = reference_plane

    # Now we will prepare the interval of longitudes and see the configuration:
    lamb_initial = lamb_star
    lamb_final = lamb_initial + theta_end*Consts.deg
    lambs = np.linspace(lamb_initial,lamb_final,n_theta)

    # Determine starting position in eccentric orbit, default is apocenter
    if e >= 0.05:
        lq = np.linspace(0,2*np.pi,2000)
        d_max = 0
        lambq_max = 0
        for ii,l in enumerate(lq):
            RP.lambq = l
            RP.changeStellarPosition(lamb_initial)
            d = RP.rstar
            if d >= d_max:
                d_max = d
                lambq_max = l
        RP.lambq = lambq_max + lambq_offset*np.pi/180     

    # Initialise the starting position
    RP.changeObserver([lamb_obs,beta_obs])
    RP.changeStellarPosition(lamb_initial)
    RP.updateSpangles()

    # Save images showing the starting position of planet, ring and star
    ecl_fig,obs_fig,star_fig = RP.plotRingedPlanet(showstar=True,showfig=False)
    for rf in "ecl","obs","star":
        fig=eval(f"{rf}_fig")
        save_figure(fig,f"fig_with_{name}_{value}_and_oi_{iorb}_rl_{alphar}_ri_{gamma}_rin_{ring_ri}_rout_{ring_re}_{rf}",
                    fig_dir=fig_dir, dpi=300, formats=["png"])

    # Simulate
    Stot  = []
    Sp    = []
    Sr    = []
    Ptot  = []
    Pp    = []
    Pr    = []
    alpha = []
    vorbit= []
    rstar = []
    rstaro= []

    print("\n########################################")
    start_msg = f" Starting orbit simulation with: \n fou_file_planet = {RP.fname_planet} ,"+\
                f"\n fou_file_ring = {RP.fname_ring} ,"+\
                f"\n orbit inclination = {iorb} ,"+\
                f"\n ring roll = {alphar} ,"+\
                f"\n ring inclination = {gamma} ,"+\
                f"\n ring ri = {RP.fi} ,"+\
                f"\n ring re = {RP.fe} ,"+\
                f"\n observer i = {RP.eobs_ecl[1]/Consts.deg} ,"+\
                f"\n observer longitude = {RP.eobs_ecl[0]/Consts.deg} ,"+\
                f"\n ring inclination wrt ecl = {RP.i/Consts.deg} ,"+\
                f"\n number of ring spangles = {RP.Nrt} ,"+\
                f"\n number of planet spangles = {RP.Np} ,"+\
                f"\n ring opacity = {tau_ring} , "+\
                f"\n reference plane = {RP.reference_plane} "
    print(start_msg)
    print("########################################\n")

    # Start the orbit
    for lamb in lambs:
        RP.changeStellarPosition(lamb)
        #print("True anomaly: ", (lamb-lamb_initial)/Consts.deg)
        RP.updateSpangles()
        RP.updateReflection(taur=tau_ring,normalize=normalize)
        #print("used ring opacity: ", RP.taur)

        # Save the relevant data
        Stot  += [RP.Stot]
        Sp    += [RP.Stotp]
        Sr    += [RP.Stotr]
        Ptot  += [RP.Ptot]
        Pp    += [RP.Ptotp]
        Pr    += [RP.Ptotr]
        alpha += [np.arccos(RP.alphaps)/Consts.deg]
        vorbit+= [RP.vorbit]
        rstar += [RP.rstar]
        rstaro+= [RP.rstar_obs]

    true_anomaly = list((lambs-lamb_initial)/Consts.deg)

    # Save the output to a pickle file for revovery after running 
    parameters = {"Fourier planet": RP.fname_planet, "Fourier ring": RP.fname_ring,
                  "Orbit inclinatie": iorb, "Ring roll": alphar,
                  "Ring inclination": ir, "Ring ri": RP.fi, "Ring re": RP.fe,
                  "Observer inclination": RP.eobs_ecl[1]/Consts.deg, 
                  "Observer longitude": RP.eobs_ecl[0]/Consts.deg,
                  "Ring i_ecl": RP.i/Consts.deg,
                  "Nr": RP.Nrt,
                  "Np": RP.Np,
                  "Ring opacity": tau_ring,
                  "Ref plane": RP.reference_plane,
                  "Planet radius": RP.Rp}
    save_dict = {"lambda": true_anomaly, "alpha": alpha, "Stot": Stot,
                 "Sp": Sp, "Sr": Sr, "Ptot": Ptot, "Pp": Pp, "Pr": Pr,
                 "vorbit": vorbit, "rstar": rstar, "rstar_obs": rstaro,
                 "Param": parameters}

    # Pickle the data, if file already exists it will be overwritten
    with open(tmp_dir + f"data_with_{name}_{value}.pkl", "wb") as f:
        pickle.dump(save_dict, f)

#Options to Parallel: https://joblib.readthedocs.io/en/stable/generated/joblib.Parallel.html
iorbs = np.array([0,10,20,30,40,50,60,70,80,90])
backend='loky'
#backend='sequential'
tini=time.time()
Parallel(n_jobs=cpu_count()-2,backend=backend,verbose=50)(delayed(simulation)(iorb) for iorb in iorbs)
tend=time.time()
print(f"Parallel execution time = {tend-tini}")