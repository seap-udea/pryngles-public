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
    def test_planet(self):
        S=Star()

        #Check exception: primary is mandatory for planets
        self.assertRaises(ValueError,lambda:Planet())

        P=Planet(primary=S)
        
        print(P.physics)
        print(P.hash)
        
        #Check derived properties
        self.assertEqual(np.isclose([P.physics.wrot],
                                    [2*np.pi/PlanetDefaults.physics["prot"]],
                                    rtol=1e-7),
                         [True]*1)
        
        P.update_body(orbit=dict(a=5),physics=dict(rho=0.2))
        print(P.orbit,P.physics)
        
        #Check exception: primary could not be different from None or Body
        self.assertRaises(AssertionError,lambda:Planet(primary="Nada"))
        
    def test_sp(self):
        S=Star()
        P=Planet(primary=S)
        P.spangle_body()
        print_df(P.sp.data.tail())
        P.sp.plot3d()
        

if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
