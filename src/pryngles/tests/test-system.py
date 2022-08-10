##################################################################
#                                                                #
#.#####...#####...##..##..##..##...####...##......######...####..#
#.##..##..##..##...####...###.##..##......##......##......##.....#
#.#####...#####.....##....##.###..##.###..##......####.....####..#
#.##......##..##....##....##..##..##..##..##......##..........##.#
#.##......##..##....##....##..##...####...######..######...####..#
#................................................................#
#                                                                #
# PlanetaRY spanGLES                                             #
# The bright-side of the light-curve of (ringed) exoplanets      #
#                                                                #
##################################################################
# Jorge I. Zuluaga, Mario Sucerquia, Jaime A. Alvarado (C) 2022  #
##################################################################
import unittest
from pryngles import *
class Test(unittest.TestCase):
    def test_system_init(self):
        
        sys=System()
        print(sys.nbodies)
        print(sys._sim.G)
        print(sys.ul,sys.um,sys.ut)
        
        sys=System(units=['m','kg','s'])
        print(sys.nbodies)
        print(sys._sim.G)
        print(sys.ul,sys.um,sys.ut)
        
        S=Star()
        sys=System(stars=S)
        print(sys.stars)
        print(sys.hashes)
        print(sys.nbodies,sys.nstars,sys.nplanets,sys.nrings,sys.nobservers)
        
        S2=Star()
        P=Planet(primary=S2)
        #Check that when planet does not use a star of the system an exception is raised
        self.assertRaises(AssertionError,lambda:System(stars=S,planets=P))
        
        P=Planet(primary=S)
        sys=System(stars=S,planets=P)
        print(sys.hashes)
        print(sys.nbodies,sys.nstars,sys.nplanets,sys.nrings,sys.nobservers)
        
        R=Planet(primary=P)
        #Check that planet cannot be initialized as a planet
        self.assertRaises(AssertionError,lambda:System(stars=S,planets=P,rings=R))

        R=Ring(primary=P)
        sys=System(stars=S,planets=P,rings=R)
        print(sys.hashes)
        print(sys.nbodies,sys.nstars,sys.nplanets,sys.nrings,sys.nobservers)
        
        
        P1=Planet(primary=S)
        P2=Planet(primary=S)
        P3=Planet(primary=S)
        sys=System(stars=S,planets=[P1,P2,P3])
        print(sys.hashes)
        print(sys.nbodies,sys.nstars,sys.nplanets,sys.nrings,sys.nobservers)
        
        """
        self.assertEqual(np.isclose([P.physics.wrot],
                                    [2*np.pi/PlanetDefaults.physics["prot"]],
                                    rtol=1e-7),
                         [True]*1)
        #Check exception: primary could not be different from None or Body
        self.assertRaises(AssertionError,lambda:Observer(primary="Nada"))
        """
        
    def test_system_add(self):
        sys=System()
        S=sys.add(kind="Star",orbit=dict(m=2))
        print(sys.hashes)
        print(sys.nbodies,sys.nstars,sys.nplanets,sys.nrings,sys.nobservers)
        print(sys.stars[0].orbit)
        print(S.orbit)
        
        S.update_body(orbit=dict(m=3))
        print(sys.stars[0].orbit)
        
        """
        self.assertEqual(np.isclose([P.physics.wrot],
                                    [2*np.pi/PlanetDefaults.physics["prot"]],
                                    rtol=1e-7),
                         [True]*1)
        #Check exception: primary could not be different from None or Body
        self.assertRaises(AssertionError,lambda:Observer(primary="Nada"))
        """
        
    def test_system_remove(self):
        sys=System()
        S=sys.add(kind="Star",orbit=dict(m=2))
        print(sys.hashes)
        print(sys.nbodies,sys.nstars,sys.nplanets,sys.nrings,sys.nobservers)

        sys.remove(body_hash=S.hash)
        print(sys.hashes)
        print(sys.nbodies,sys.nstars,sys.nplanets,sys.nrings,sys.nobservers)
        
        sys=System()
        S=sys.add(kind="Star")
        P=sys.add(kind="Planet",primary=S)
        R=sys.add(kind="Ring",primary=P)
        print(sys.hashes)
        print(sys.nbodies,sys.nstars,sys.nplanets,sys.nrings,sys.nobservers)

        sys.remove(body_hash=S.hash)
        print(sys.hashes)
        print(sys.nbodies,sys.nstars,sys.nplanets,sys.nrings,sys.nobservers)

        sys=System()
        S=sys.add(kind="Star")
        P=sys.add(kind="Planet",primary=S)
        R=sys.add(kind="Ring",primary=P)
        print(sys.hashes)
        print(sys.nbodies,sys.nstars,sys.nplanets,sys.nrings,sys.nobservers)

        sys.remove(body_hash=P.hash)
        print(sys.hashes)
        print(sys.nbodies,sys.nstars,sys.nplanets,sys.nrings,sys.nobservers)
        
        sys=System()
        S=sys.add(kind="Star")
        P=sys.add(kind="Planet",primary=S)
        R=sys.add(kind="Ring",primary=P)
        print(sys.hashes)
        print(sys.nbodies,sys.nstars,sys.nplanets,sys.nrings,sys.nobservers)

        sys.remove(body_hash=R.hash)
        print(sys.hashes)
        print(sys.nbodies,sys.nstars,sys.nplanets,sys.nrings,sys.nobservers)
        
        """
        self.assertEqual(np.isclose([P.physics.wrot],
                                    [2*np.pi/PlanetDefaults.physics["prot"]],
                                    rtol=1e-7),
                         [True]*1)
        #Check exception: primary could not be different from None or Body
        self.assertRaises(AssertionError,lambda:Observer(primary="Nada"))
        """
        

if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
