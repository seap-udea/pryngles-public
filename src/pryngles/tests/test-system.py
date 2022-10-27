##################################################################
#                                                                #
#.#####...#####...##..##..##..##...####...##......######...####..#
#.##..##..##..##...####...###.##..##......##......##......##.....#
#.#####...#####.....##....##.###..##.###..##......####.....####..#
#.##......##..##....##....##..##..##..##..##......##..........##.#
#.##......##..##....##....##..##...####...######..######...####..#
#................................................................#

# PlanetaRY spanGLES                                             #
#                                                                #
##################################################################
# License http://github.com/seap-udea/pryngles-public            #
##################################################################
# Main contributors:                                             #
#   Jorge I. Zuluaga, Mario Sucerquia, Jaime A. Alvarado         #
##################################################################
import unittest
from pryngles import *
class Test(unittest.TestCase):
	def test_system_init(self):
	
	    Verbose.VERBOSITY=VERB_ALL
	
	    sys=System(resetable=True)
	    print("Nbodies = ",sys.nbodies)
	    print("G constant = ",sys.G)
	    print("G constant = ",sys.units)
	    print("Canonical units = ",sys.ul,sys.um,sys.ut)
	
	    sys=System(units=['m','kg','s'])
	    print("Nbodies = ",sys.nbodies)
	    print("G constant = ",sys.G)
	    print("G constant = ",sys.units)
	    print("Canonical units = ",sys.ul,sys.um,sys.ut)
	    print(sys)
	
	    sys.save_to("/tmp/system.pkl")
	    print(sys.status())
	    sys2=System("/tmp/system.pkl")
	    print(sys2.status())
	
	    Verbose.VERBOSITY=VERB_NONE
	
	def test_system_add(self):
	    global sys
	
	    Verbose.VERBOSITY=VERB_ALL
	
	    sys=System()
	    S=sys.add(m=8,radius=4,x=5,vy=2)
	    P=sys.add("Planet",parent=S,radius=2,x=1)
	    M=sys.add("Planet",name="Moon",parent=P,radius=2,x=2)
	    R=sys.add("Ring",parent=P,fi=1.3,fe=2.3)
	    print(sys)
	    print(sys.root)
	
	    """
	    for particle in sys.sim.particles:
	        print(particle)
	
	    print(len(sys.bodies),len(sys.sim.particles))
	    """
	
	    Verbose.VERBOSITY=VERB_NONE
	
	def test_system_remove(self):
	
	    Verbose.VERBOSITY=VERB_ALL
	
	    sys=System()
	    S=sys.add(name="Star",m=8,radius=4,x=5,vy=2)
	    P=sys.add("Planet",parent=S,name="Planet",radius=2,x=2)
	    M=sys.add("Planet",parent=P,name="Moon",radius=2,x=2)
	    R=sys.add("Ring",parent=P,name="Ring",fi=1.3,fe=2.3)
	    print(sys.bodies)
	    sys.remove("Ring")
	    print(sys.bodies)
	    sys.remove("Planet")
	    print(sys.bodies)
	    sys.remove("Star")
	    print(sys.bodies)
	
	    Verbose.VERBOSITY=VERB_NONE
	
	def test_sim(self):
	
	    global sys
	    plt.close("all")
	
	    Verbose.VERBOSITY=VERB_ALL
	
	    #Create system
	    sys=System(units=['au','msun','yr'])
	    S=sys.add(m=8,radius=4)
	    P1=sys.add("Planet",name="Planet1",parent=S,radius=2,a=1,M=90*Consts.deg,inc=70*Consts.deg)
	    M1P1=sys.add("Planet",name="Moon1P1",parent=P1,radius=2,a=0.1)
	    R=sys.add("Ring",name="Ring",parent=P1,radius=2)
	    P2=sys.add("Planet",name="Planet2",parent=S,radius=2,a=2,M=0*Consts.deg,inc=0*Consts.deg)
	    S.show_tree()
	    
	    #Initialize
	    orbit=sys.initialize_simulation(orbital_tree=[[S,[P1,M1P1]],P2])
	    sys.sim.status()
	    
	    #Check save to disk
	    sys.save_to("/tmp/system.pkl")
	    sys=System()
	    sys.load_from("/tmp/system.pkl")
	    sys.sim.status()
	    
	    #Animate
	    Plot.animate_rebound(sys.sim)
	    
	    Verbose.VERBOSITY=VERB_NONE
	
	def test_obs(self):
	
	    global sys
	    plt.close("all")
	
	    Verbose.VERBOSITY=VERB_NONE
	
	    nspangles=100
	
	    #Define system
	    sys=System(resetable=True)
	
	    #Add objects
	    S=sys.add(nspangles=nspangles,m=8,radius=1)
	    P=sys.add("Planet",parent=S,nspangles=nspangles,m=1,radius=0.2,a=5)
	
	    #Test setting observer without spangling
	    self.assertRaises(AssertionError,lambda:sys._set_observer(nvec=[1,0,0]))
	
	    #Spangle system
	    sys.initialize_simulation()
	    sys.spangle_system()
	
	    sys._set_observer(nvec=[-1,0,0])
	    sys.sg.plot3d()
	
	    sys._set_observer(nvec=[0,0,1])
	    sys.sg.plot3d()
	
	    Verbose.VERBOSITY=VERB_NONE
	
	def test_setluz(self):
	
	    global sys
	    plt.close("all")
	
	    Verbose.VERBOSITY=VERB_NONE
	    nspangles=500
	    sys=System()
	    S=sys.add("Star",nspangles=nspangles,m=1,radius=1)
	    D=sys.add("Ring",name="Disk",parent=S,nspangles=nspangles,fi=20,fe=30,i=0*Consts.deg)
	    P=sys.add("Planet",parent=S,nspangles=nspangles,radius=0.2,m=1e-3,a=10)
	    R=sys.add("Ring",parent=P,nspangles=nspangles,fi=1.5,fe=2.0,i=-20*Consts.deg)
	    M=sys.add("Planet",parent=P,name="Moon",nspangles=nspangles,radius=0.1,m=1e-6,a=1,M=30*Consts.deg)
	    K=sys.add("Ring",name="Cronoring",parent=M,nspangles=nspangles,fi=1.1,fe=1.5,i=20*Consts.deg)
	    
	    sys.initialize_simulation()
	    sys.spangle_system()
	    
	    #sys.sg.plot3d(center_at="Ring",not_plot=["Disk"])
	    #sys.sg.plot3d(center_at="Ring")
	    sys.sg.plot3d()
	    cond=(sys.sg.data.name=="Moon")&(sys.sg.data.hidden_by_luz!="")
	    print_df(sys.sg.data.loc[cond,["hidden_by_luz"]].head(10))
	
	    Verbose.VERBOSITY=VERB_NONE
	
	def test_spangle(self):
	
	    global sys
	
	    Verbose.VERBOSITY=VERB_NONE
	
	    nspangles=100
	    sys=System(resetable=False)
	    S=sys.add(name="Star",nspangles=nspangles,m=9,radius=1)
	    P=sys.add("Planet",parent=S,name="Planet",nspangles=nspangles,radius=0.2,a=2)
	    M=sys.add("Planet",parent=P,name="Moon",nspangles=nspangles,radius=0.1,a=1)
	    R=sys.add("Ring",parent=P,name="Ring",nspangles=nspangles,fi=1.3,fe=2.3,i=90*Consts.deg)
	
	    sys.initialize_simulation()
	    sys.spangle_system()
	
	    #Check addition columns
	    print(sys.source)
	    print(sys.sg.data.columns)
	
	    #Check save
	    sys.save_to("/tmp/system.pkl")
	
	    #Check plot
	    #sys.sp.plot3d(center_at="Ring",not_plot=["Star1","Star2"])
	    sys.sg.plot3d()
	
	    Verbose.VERBOSITY=VERB_NONE
	
	def test_update(self):
	
	    Verbose.VERBOSITY=VERB_NONE
	
	    nspangles=100
	    sys=System(resetable=True)
	    S=sys.add("Star",name="Star",nspangles=nspangles,m=8,radius=1)
	    P=sys.add("Planet",parent=S,name="Planet",nspangles=nspangles,radius=0.2,a=2)
	    M=sys.add("Planet",parent=P,name="Moon",nspangles=nspangles,radius=0.1,a=1)
	    R=sys.add("Ring",parent=P,name="Ring",nspangles=nspangles,fi=1.3,fe=2.3,i=90*Consts.deg)
	    print(P.radius)
	    sys.update_body(P,radius=0.5)
	    print(P.radius)
	    sys.update_body("Ring",fe=3.0)
	    print(R.radius)
	    sys.initialize_simulation()
	    sys.spangle_system()
	    self.assertRaises(AssertionError,lambda:sys.update_body("Ring",fe=3.0))
	
	    Verbose.VERBOSITY=VERB_NONE
	
	def test_reset(self):
	
	    global sys
	
	    Verbose.VERBOSITY=VERB_NONE
	
	    nspangles=100
	    sys=System(resetable=True)
	    S=sys.add("Star",name="Star",nspangles=nspangles,m=8,radius=1)
	    P=sys.add("Planet",parent=S,name="Planet",nspangles=nspangles,radius=0.2,a=2)
	    M=sys.add("Planet",parent=P,name="Moon",nspangles=nspangles,radius=0.1,a=1)
	    R=sys.add("Ring",parent=P,name="Ring",nspangles=nspangles,fi=1.3,fe=2.3,i=90*Consts.deg)
	    
	    sys.initialize_simulation()
	    sys.spangle_system()
	
	    #All transformations from here are not stored
	    sys.sg.plot3d()
	    sys._set_observer(nvec=[0,0,-1])
	    sys.sg.plot3d()
	
	    #All transformations from here are not stored
	    sys.reset()
	    sys.sg.plot3d()
	
	    Verbose.VERBOSITY=VERB_NONE
	
	def test_int(self):
	
	    global sys
	
	    Verbose.VERBOSITY=VERB_NONE
	
	    nspangles=100
	    sys=System()
	    S=sys.add("Star",name="Star",nspangles=nspangles,m=1,radius=1)
	    M=sys.add("Planet",parent=S,name="Moon",nspangles=nspangles,radius=0.1,m=1e-6,a=3)
	    P=sys.add("Planet",parent=S,name="Planet",nspangles=nspangles,radius=0.2,m=1e-3,a=5)
	    R=sys.add("Ring",parent=P,name="Ring",nspangles=nspangles,fi=1.3,fe=2.3,i=20*Consts.deg)
	    
	    sys.initialize_simulation()
	    sys.spangle_system()
	
	    sys.integrate(1)
	
	    sys._set_observer([0,0,1])
	    sys._set_luz()
	
	    sys.integrate(1)
	
	    sys._set_observer([0,0,1])
	    sys._set_luz()
	
	    sys.sg.plot3d()
	    sys.sg.plot3d(center_at="Ring")
	
	    Verbose.VERBOSITY=VERB_NONE
	
	def test_legacy(self):
	
	    global sys
	
	    Verbose.VERBOSITY=VERB_NONE
	    #"""
	    sys=System()
	    S=sys.add(kind="Star",
	              physics=dict(radius=Consts.rsun/sys.ul),
	              optics=dict(limb_coeffs=[0.65])
	             )
	    P=sys.add(kind="Planet",primary=S,
	              orbit=dict(a=0.2,e=0.0),
	              physics=dict(radius=Consts.rsaturn/sys.ul)
	             )
	    R=sys.add(kind="Ring",primary=P,
	              physics=dict(fi=1.5,fe=2.5,i=30*Consts.deg)
	             )
	    O=sys.add(kind="Observer",
	              optics=dict(lamb=90*Consts.deg,beta=90*Consts.deg)
	             )
	    RP=sys.ensamble_system()
	    ecliptic,observer,star=RP.plotRingedPlanet(showfig=1)
	    
	    Verbose.VERBOSITY=VERB_NONE
	
if __name__=="__main__":
    unittest.main(argv=['first-arg-is-ignored'],exit=False)
    