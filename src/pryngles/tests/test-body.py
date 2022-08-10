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
    def test_fun(self):
        B=Body(BodyDefaults,"Body",None,dict(x=0),dict(mass=2.3),dict())
        print(B)
        print(B.orbit)
        print(B.physics)
        print(B.optics)
        
        B.update_body(orbit=dict(m=2))
        print(B.orbit)
        
        
        C=Body(BodyDefaults,"Body",B,dict(),dict(),dict())
        print(C)
        print(B)
        """
        self.assertEqual(np.isclose([P.physics.wrot],
                                    [2*np.pi/PlanetDefaults.physics["prot"]],
                                    rtol=1e-7),
                         [True]*1)
        self.assertRaises(AssertionError,lambda:Observer(primary="Nada"))
        """
        

if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
